"""
Sample metadata for the design.

Any CSV/TSV/XLSX sheet that is not the count table can describe the samples: a
user sample sheet, rnaseq-quant's samples.csv, or SPARC samples.xlsx /
subjects.xlsx. Column names become snake_case ('pool id' -> pool_id). The sheet
and column whose values name the most count samples is used (values compared
without case and without a sam-/sub-/pool- prefix, so 'pool-FN1', 'FN1' and
'fn1' match); SAMPLE_SHEET restricts the choice to files of that name.

When that sheet has a subject_id column and another sheet (SPARC subjects.xlsx)
is keyed by subject_id, the subject columns are joined on (prefixed subject_
when a name clashes). Several rows can match one library (SPARC samples.xlsx
lists the animals of a pooled library): a column is kept for that library when
all its rows agree, and left empty otherwise.

Only the columns the design names are ever written out, so clinical fields in
a sheet never reach the outputs.
"""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from processor import tables as T

log = logging.getLogger("bulk-de")
PREFIX = re.compile(r"^(sam|sub|pool)-", re.I)
SPECIES = {"homo sapiens": "human", "human": "human", "mus musculus": "mouse", "mouse": "mouse",
           "rattus norvegicus": "rat", "rat": "rat"}


def key(v) -> str:
    return PREFIX.sub("", str(v).strip()).lower()


@dataclass
class Metadata:
    table: pd.DataFrame                       # index: count sample ids; string columns ("" = missing)
    sheet: str | None = None
    match_column: str | None = None
    subjects_sheet: str | None = None
    n_matched: int = 0
    notes: list[str] = field(default_factory=list)
    skipped: dict[str, str] = field(default_factory=dict)


def _sheets(files: list[Path], exclude: set[Path], wanted: str) -> tuple[list[tuple[Path, pd.DataFrame]], dict]:
    out, skipped = [], {}
    for p in files:
        if p in exclude or not T.is_table(p) or p.name.endswith((".genes.results", "ReadsPerGene.out.tab")):
            continue
        if wanted and p.name != wanted:
            continue
        if p.stat().st_size > T.MAX_SHEET_BYTES:
            continue
        try:
            df = T.read_sheet(p)
        except Exception as exc:  # an unreadable sheet is not fatal
            skipped[str(p)] = f"{type(exc).__name__}: {exc}"
            continue
        if not df.empty and df.shape[1] >= 2:
            out.append((p, df))
    return out, skipped


def find(files: list[Path], exclude: set[Path], samples: list[str], wanted: str) -> Metadata:
    sheets, skipped = _sheets(files, exclude, wanted)
    wanted_keys = {key(s): s for s in samples}
    best = None
    for p, df in sheets:
        for col in df.columns:
            hits = set(df[col].map(key)) & set(wanted_keys)
            score = (len(hits), p.name.lower().startswith("samples"), col == "sample_id")
            if hits and (best is None or score > best[0]):
                best = (score, p, df, col)
    empty = pd.DataFrame(index=pd.Index(samples, name="sample_id"))
    if best is None:
        return Metadata(empty, notes=["no sheet names the count samples; no metadata"], skipped=skipped)
    (n, _, _), path, df, col = best
    meta = Metadata(empty, sheet=path.name, match_column=col, n_matched=n, skipped=skipped)

    if "subject_id" in df.columns:
        subj = next(((p, s) for p, s in sheets if p != path and "subject_id" in s.columns
                     and "sample_id" not in s.columns), None)
        if subj is not None:
            sp, sdf = subj
            sdf = sdf.drop_duplicates("subject_id")
            sdf = sdf.rename(columns={c: f"subject_{c}" for c in sdf.columns if c in df.columns and c != "subject_id"})
            df = df.merge(sdf, on="subject_id", how="left")
            meta.subjects_sheet = sp.name

    df = df.assign(_key=df[col].map(key))
    df = df[df["_key"].isin(wanted_keys)]
    rows, mixed = {}, set()
    for k, g in df.groupby("_key", sort=False):
        row = {}
        for c in g.columns.drop("_key"):
            vals = set(g[c]) - {""}
            if len(vals) == 1:
                row[c] = vals.pop()
            elif len(vals) > 1:
                row[c] = ""
                mixed.add(c)
        rows[wanted_keys[k]] = row
    table = pd.DataFrame.from_dict(rows, orient="index").reindex(samples).fillna("")
    table.index.name = "sample_id"
    meta.table = table.drop(columns=[c for c in ("sample_id",) if c in table.columns])
    if mixed:
        meta.notes.append(f"{len(mixed)} columns differ between rows of the same library and were left empty "
                          f"for it (pooled animals, several samples of one subject)")
    log.info("metadata: %s, column %s matches %d of %d samples%s", path.name, col, n, len(samples),
             f", subjects from {meta.subjects_sheet}" if meta.subjects_sheet else "")
    return meta


def organism(table: pd.DataFrame) -> str | None:
    for c in ("organism", "species", "subject_species"):
        if c in table.columns:
            vals = {SPECIES.get(str(v).strip().lower()) for v in table[c]} - {None}
            if len(vals) == 1:
                return vals.pop()
    return None
