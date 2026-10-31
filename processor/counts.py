"""
Gene counts from any of the supported layouts, converted to one genes x samples
integer matrix (the shape of rnaseq-quant's gene_counts.csv).

Tried in this order; the first layout found is used and every other count file
is listed as ignored:

  rnaseq-quant   gene_counts.csv: gene_id, gene_name, one column per sample
  featureCounts  a table whose header starts Geneid, Chr, Start, End, Strand, Length
  matrix         any other CSV/TSV/XLSX with >= 100 rows and >= 2 numeric columns:
                 genes in rows, the first column the gene id, an optional text column
                 with gene names, one numeric column per sample (a title row above the
                 header is skipped)
  STAR           <sample>ReadsPerGene.out.tab, one per sample (column chosen by strandedness)
  RSEM           <sample>.genes.results, one per sample (expected_count)
  two-column     one file per sample, gene and count, as htseq-count writes it
                 (also the authors' *_Genes_ReadCount.txt of Discover 427)

Per-sample files are named by the nearest sam-/pool- folder above them (SPARC
layout), else by the file name without its tool suffix. Non-integer counts
(RSEM, salmon) are rounded, as tximport does. Summary rows (htseq's __*, STAR's
N_*) are dropped.
"""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from processor import tables as T

log = logging.getLogger("bulk-de")

FEATURECOUNTS_HEAD = ["geneid", "chr", "start", "end", "strand", "length"]
NOT_SAMPLES = {"length", "effective_length", "chr", "chrom", "start", "end", "strand", "gene_length"}
NAME_HINT = re.compile(r"name|symbol|^gene$", re.I)
SAMPLE_SUFFIX = re.compile(r"(ReadsPerGene\.out\.tab|\.genes\.results|_Genes_ReadCount|\.htseq|\.counts?|"
                           r"[._-]count(s)?|[._-]readcounts?)$", re.I)
MIN_MATRIX_ROWS = 100


@dataclass
class Counts:
    matrix: pd.DataFrame                    # genes x samples, int64; index gene_id
    gene_names: pd.Series                   # gene_id -> gene name ("" when unknown)
    layout: str
    files: list[Path]
    notes: list[str] = field(default_factory=list)
    skipped: dict[str, str] = field(default_factory=dict)


class NoCounts(Exception):
    pass


def _ints(df: pd.DataFrame, what: str, notes: list[str]) -> pd.DataFrame:
    vals = df.apply(lambda c: pd.to_numeric(c.astype(str).str.strip(), errors="coerce"))
    if vals.isna().any().any():
        bad = int(vals.isna().any(axis=1).sum())
        notes.append(f"{what}: {bad} rows with a missing or non-numeric count were dropped")
        vals = vals.dropna()
    if (vals < 0).any().any():
        raise ValueError("negative values: not counts")
    if not np.allclose(vals.to_numpy(), np.round(vals.to_numpy())):
        notes.append(f"{what}: non-integer counts (estimated counts) were rounded")
    return vals.round().astype(np.int64)


def sample_name(path: Path, root: Path) -> str:
    for parent in path.relative_to(root).parents:
        if parent.name.lower().startswith(("sam-", "pool-")):
            return parent.name
    name = path.name
    for ext in (".gz", ".txt", ".tsv", ".tab", ".csv"):
        if name.lower().endswith(ext) and not name.endswith("ReadsPerGene.out.tab"):
            name = name[: -len(ext)]
    return SAMPLE_SUFFIX.sub("", name) or path.stem


# ---- matrix layouts ----------------------------------------------------------------------------

def _matrix(df: pd.DataFrame, what: str, notes: list[str]) -> tuple[pd.DataFrame, pd.Series]:
    """df with a header: first column = gene id, numeric columns = samples, a text column = names."""
    dup = df.columns.duplicated()
    if dup.any():
        notes.append(f"{what}: columns named twice ({', '.join(map(str, df.columns[dup]))}); the second copy was ignored")
        df = df.loc[:, ~dup]
    id_col = df.columns[0]
    rest = [c for c in df.columns[1:] if T.snake(c) not in NOT_SAMPLES]
    sample_cols = [c for c in rest if T.numeric_share(df[c]) >= 0.95]
    text_cols = [c for c in rest if c not in sample_cols]
    name_col = next((c for c in text_cols if NAME_HINT.search(str(c))), text_cols[0] if text_cols else None)
    if len(sample_cols) < 2:
        raise ValueError("fewer than two numeric sample columns")
    ids = df[id_col].astype(str).str.strip()
    keep = (ids != "") & ~ids.str.startswith(("__", "N_"))
    df = df[keep]
    names = df[name_col].astype(str).str.strip() if name_col is not None else df[id_col].astype(str).str.strip()
    counts = _ints(df[sample_cols], what, notes)
    names = names.loc[counts.index]
    counts.index = df.loc[counts.index, id_col].astype(str).str.strip().to_numpy()
    names.index = counts.index
    if counts.index.duplicated().any():
        n = int(counts.index.duplicated().sum())
        notes.append(f"{what}: {n} duplicated gene ids; their counts were summed")
        counts = counts.groupby(level=0, sort=False).sum()
    names = names[~names.index.duplicated()].reindex(counts.index).fillna("")
    counts.columns = [str(c).strip() for c in counts.columns]
    return counts, names


def read_quant(path: Path, notes: list[str]):
    df = T.with_header(T.read_raw(path))
    if list(df.columns[:2]) != ["gene_id", "gene_name"]:
        raise ValueError("gene_counts.csv without gene_id, gene_name columns")
    return _matrix(df, path.name, notes)


def is_featurecounts(raw: pd.DataFrame) -> bool:
    return raw.shape[1] > 6 and [str(v).strip().lower() for v in raw.iloc[0, :6]] == FEATURECOUNTS_HEAD


def read_featurecounts(raw: pd.DataFrame, what: str, notes: list[str]):
    df = raw.iloc[1:].reset_index(drop=True)
    df.columns = [str(v).strip() for v in raw.iloc[0]]
    df = df.drop(columns=list(df.columns[1:6]))
    df.columns = [df.columns[0]] + [_bam_name(c) for c in df.columns[1:]]
    return _matrix(df, what, notes)


def _bam_name(c: str) -> str:
    c = Path(c).name
    return re.sub(r"([._]?(Aligned\.sortedByCoord\.out|Aligned\.out|sorted|dedup))?\.bam$", "", c) or c


# ---- per-sample layouts ------------------------------------------------------------------------

def read_star(path: Path, notes: list[str]) -> tuple[pd.Series, str]:
    """STAR --quantMode GeneCounts: col 2 unstranded, 3 forward, 4 reverse. A stranded library puts
    most reads in one of 3/4; we take that column when it holds >= 3x the other, else column 2."""
    raw = T.read_raw(path)
    raw = raw[~raw[0].astype(str).str.startswith("N_")]
    vals = raw.iloc[:, 1:4].apply(pd.to_numeric, errors="raise").astype(np.int64)
    s3, s4 = int(vals.iloc[:, 1].sum()), int(vals.iloc[:, 2].sum())
    if max(s3, s4) >= 3 * max(min(s3, s4), 1):
        col, how = (1, "forward-stranded (column 3)") if s3 > s4 else (2, "reverse-stranded (column 4)")
    else:
        col, how = 0, "unstranded (column 2)"
    return pd.Series(vals.iloc[:, col].to_numpy(), index=raw[0].astype(str).to_numpy()), how


def read_rsem(path: Path, notes: list[str]) -> tuple[pd.Series, str]:
    df = T.with_header(T.read_raw(path))
    if "gene_id" not in df.columns or "expected_count" not in df.columns:
        raise ValueError("no gene_id / expected_count columns")
    v = pd.to_numeric(df["expected_count"], errors="raise")
    return pd.Series(np.round(v.to_numpy()).astype(np.int64), index=df["gene_id"].astype(str).to_numpy()), \
        "expected_count, rounded"


def read_two_column(raw: pd.DataFrame) -> pd.Series:
    if T.numeric_share(raw[1].iloc[:1]) == 0:      # a header row
        raw = raw.iloc[1:]
    ids = raw[0].astype(str).str.strip()
    raw = raw[~ids.str.startswith("__")]
    v = pd.to_numeric(raw[1], errors="raise")
    if (v < 0).any():
        raise ValueError("negative values: not counts")
    return pd.Series(np.round(v.to_numpy()).astype(np.int64), index=raw[0].astype(str).str.strip().to_numpy())


# ---- detection ---------------------------------------------------------------------------------

def load(files: list[Path], root: Path) -> Counts:
    tabs = [p for p in files if T.is_table(p) or p.name.endswith(("ReadsPerGene.out.tab", ".genes.results"))]
    notes: list[str] = []
    skipped: dict[str, str] = {}

    quant = [p for p in tabs if p.name == "gene_counts.csv"]
    if quant:
        m, names = read_quant(quant[0], notes)
        return _done(m, names, "rnaseq-quant", [quant[0]], notes, skipped)

    star = [p for p in tabs if p.name.endswith("ReadsPerGene.out.tab")]
    rsem = [p for p in tabs if p.name.endswith(".genes.results")]
    featurecounts, matrices, two_col = [], [], []
    for p in tabs:
        if p in star or p in rsem or p.name.lower() in ("samples.csv", "gene_tpm.csv"):
            continue
        try:
            raw = T.read_raw(p)
        except Exception as exc:  # an unreadable file is not fatal
            skipped[str(p)] = f"{type(exc).__name__}: {exc}"
            continue
        if raw.empty:
            continue
        if is_featurecounts(raw):
            featurecounts.append((p, raw))
        elif raw.shape[1] == 2 and len(raw) >= MIN_MATRIX_ROWS and T.numeric_share(raw[1].iloc[1:200]) >= 0.95:
            two_col.append((p, raw))
        elif len(raw) >= MIN_MATRIX_ROWS and raw.shape[1] >= 3:
            matrices.append((p, raw))

    if featurecounts:
        p, raw = featurecounts[0]
        m, names = read_featurecounts(raw, p.name, notes)
        return _done(m, names, "featureCounts", [p], notes, skipped)
    best = None
    for p, raw in matrices:
        try:
            m, names = _matrix(T.with_header(raw), p.name, notes_m := [])
        except ValueError as exc:
            skipped[str(p)] = str(exc)
            continue
        integer = not any("rounded" in n for n in notes_m)
        key = (integer, len(m))
        if best is None or key > best[0]:
            best = (key, p, m, names, notes_m)
    if best is not None:
        _, p, m, names, notes_m = best
        return _done(m, names, "matrix", [p], notes + notes_m, skipped)

    for layout, items, reader in (("STAR", star, read_star), ("RSEM", rsem, read_rsem),
                                  ("two-column", two_col, None)):
        if not items:
            continue
        cols: dict[str, pd.Series] = {}
        used: list[Path] = []
        hows: set[str] = set()
        for item in items:
            p = item[0] if isinstance(item, tuple) else item
            sid = sample_name(p, root)
            if sid in cols:
                skipped[str(p)] = f"a second file for sample {sid}"
                continue
            try:
                if reader is None:
                    s = read_two_column(item[1])
                else:
                    s, how = reader(p, notes)
                    hows.add(how)
            except Exception as exc:  # one bad sample is skipped
                skipped[str(p)] = f"{type(exc).__name__}: {exc}"
                log.error("%s: unreadable, skipped: %s", p, exc)
                continue
            s = s.groupby(level=0, sort=False).sum()
            cols[sid] = s
            used.append(p)
            log.info("%s: %s, %d genes, %d counts", sid, layout, len(s), int(s.sum()))
        if not cols:
            continue
        if hows:
            notes.append(f"{layout}: " + "; ".join(sorted(hows)))
        common = set.intersection(*(set(s.index) for s in cols.values()))
        union = set.union(*(set(s.index) for s in cols.values()))
        order = [g for g in next(iter(cols.values())).index if g in common]
        if len(common) < len(union):
            notes.append(f"{len(union) - len(common)} genes not in every sample's file were dropped")
        m = pd.DataFrame({k: s.reindex(order).to_numpy() for k, s in cols.items()}, index=order).astype(np.int64)
        return _done(m, pd.Series("", index=m.index), layout, used, notes, skipped)
    raise NoCounts("no count table found (gene_counts.csv, featureCounts, a genes x samples table, "
                   "STAR ReadsPerGene.out.tab, RSEM *.genes.results or per-sample two-column files)")


def _done(m, names, layout, used, notes, skipped) -> Counts:
    if m.shape[1] < 2:
        raise NoCounts(f"{layout}: fewer than two samples")
    m.index = m.index.astype(str)
    m.index.name = "gene_id"
    names = names.astype(str).str.strip()
    # Files that name genes by symbol only (htseq, the authors' tables) have no separate name column.
    symbol = ~m.index.str.match(r"^(ENS[A-Z]*G\d|\d+$)")
    names = names.where((names != "") | ~symbol, pd.Series(m.index, index=m.index))
    return Counts(m, names, layout, list(used), notes, skipped)
