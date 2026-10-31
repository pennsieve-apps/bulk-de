"""
Finding and reading the input files, untrusted.

Files are only ever opened from a walk of INPUT_DIR (no path is taken from file
contents). The walk follows symlinks, since INPUT_DIR entries are symlinks, but a
link found below an entry must stay inside that entry's own tree (resolved), and
a folder is visited once (no loops). Tables are capped in size before reading;
XLSX goes through openpyxl with defusedxml installed (no entity expansion), and
formulas are never evaluated.
"""
from __future__ import annotations

import csv
import gzip
import io
import logging
import os
import re
from pathlib import Path

import pandas as pd

log = logging.getLogger("bulk-de")

MAX_TABLE_BYTES = 1 << 30          # a count table read into memory (1 GiB)
MAX_SHEET_BYTES = 50 << 20         # a metadata sheet (50 MiB)
TABLE_SUFFIXES = (".csv", ".tsv", ".txt", ".tab", ".counts", ".xlsx", ".csv.gz", ".tsv.gz", ".txt.gz", ".tab.gz")


def _inside(path: str, root: str) -> bool:
    return path == root or path.startswith(root.rstrip(os.sep) + os.sep)


def walk(root: Path) -> list[Path]:
    """Every regular file under root, sorted, following symlinks that stay inside their top-level entry."""
    found: list[Path] = []
    root_real = os.path.realpath(root)
    for entry in sorted(root.iterdir()) if root.is_dir() else []:
        if entry.name.startswith("."):
            continue
        base = os.path.realpath(entry)
        if entry.is_file():
            found.append(entry)
            continue
        seen: set[str] = set()
        for dirpath, dirnames, filenames in os.walk(entry, followlinks=True):
            real = os.path.realpath(dirpath)
            if real in seen or not (_inside(real, base) or _inside(real, root_real)):
                dirnames[:] = []
                continue
            seen.add(real)
            dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
            for name in sorted(filenames):
                p = Path(dirpath) / name
                rp = os.path.realpath(p)
                if name.startswith(".") or not os.path.isfile(rp):
                    continue
                if not (_inside(rp, base) or _inside(rp, root_real)):
                    log.warning("%s: a link that leaves its input folder; ignored", p)
                    continue
                found.append(p)
    return found


def is_table(p: Path) -> bool:
    return p.name.lower().endswith(TABLE_SUFFIXES)


def snake(name) -> str:
    """'Subject experimental group' -> 'subject_experimental_group'."""
    s = re.sub(r"[^0-9a-zA-Z]+", "_", str(name).strip()).strip("_").lower()
    return s or "column"


def _text(path: Path, limit: int) -> str:
    if path.stat().st_size > limit:
        raise ValueError(f"larger than {limit >> 20} MiB")
    raw = path.read_bytes()
    if path.name.lower().endswith(".gz"):
        with gzip.GzipFile(fileobj=io.BytesIO(raw)) as fh:
            raw = fh.read(limit + 1)
        if len(raw) > limit:
            raise ValueError(f"larger than {limit >> 20} MiB uncompressed")
    return raw.decode("utf-8-sig", errors="replace")


def read_raw(path: Path, limit: int = MAX_TABLE_BYTES) -> pd.DataFrame:
    """Every cell as a string (no header interpretation), comment lines ('#') dropped.
    XLSX: the first sheet; CSV/TSV/TXT: comma or tab, whichever the first data line has."""
    if path.name.lower().endswith(".xlsx"):
        if path.stat().st_size > limit:
            raise ValueError(f"larger than {limit >> 20} MiB")
        return pd.read_excel(path, sheet_name=0, header=None, dtype=str, engine="openpyxl")
    text = _text(path, limit)
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if not lines:
        return pd.DataFrame()
    sep = "\t" if "\t" in lines[min(len(lines) - 1, 5)] else ","
    rows = list(csv.reader(lines, delimiter=sep))          # ragged rows (a title line) are fine
    width = max(len(r) for r in rows)
    return pd.DataFrame([r + [""] * (width - len(r)) for r in rows], dtype=str).replace("", None)


def numeric_share(col: pd.Series) -> float:
    s = col.dropna().astype(str).str.strip()
    s = s[s != ""]
    if s.empty:
        return 0.0
    return float(pd.to_numeric(s, errors="coerce").notna().mean())


def with_header(raw: pd.DataFrame) -> pd.DataFrame:
    """Use the last non-numeric row before the data as the header: deposited spreadsheets often have a
    title or study row above the column names. Data starts at the first row that is mostly numbers."""
    raw = raw.dropna(axis=1, how="all").dropna(axis=0, how="all").reset_index(drop=True)
    start = 0
    for i in range(min(len(raw), 20)):
        row = raw.iloc[i]
        vals = row.dropna().astype(str).str.strip()
        vals = vals[vals != ""]
        if len(vals) and pd.to_numeric(vals, errors="coerce").notna().mean() >= 0.5:
            start = i
            break
    header_row = max(start - 1, 0)
    header = [str(v).strip() if isinstance(v, str) and v.strip() else f"column_{j}"
              for j, v in enumerate(raw.iloc[header_row])]
    df = raw.iloc[header_row + 1:].reset_index(drop=True) if start > 0 else raw.iloc[1:].reset_index(drop=True)
    df.columns = header
    return df


def read_sheet(path: Path) -> pd.DataFrame:
    """A metadata sheet with its header row; column names in snake_case, cells stripped strings."""
    raw = read_raw(path, MAX_SHEET_BYTES).dropna(axis=1, how="all").dropna(axis=0, how="all")
    if raw.empty:
        return raw
    df = raw.iloc[1:].reset_index(drop=True)
    df.columns = _unique([snake(c) for c in raw.iloc[0].fillna("")])
    df = df.map(lambda v: "" if v is None or (isinstance(v, float) and pd.isna(v)) else str(v).strip())
    df = df.loc[:, [c for c in df.columns if (df[c] != "").any()]]
    return df


def _unique(names: list[str]) -> list[str]:
    out, seen = [], {}
    for n in names:
        k = seen.get(n, 0)
        out.append(n if k == 0 else f"{n}_{k + 1}")
        seen[n] = k + 1
    return out
