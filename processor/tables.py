"""
Reading the input tables, untrusted.

Files come from processor/safe_inputs.py's InputDir (the shared walker: links that leave
the input are skipped) and are read only through it, capped: text tables up to
MAX_TABLE_BYTES (uncompressed for .gz), metadata sheets up to MAX_SHEET_BYTES, and an
.xlsx only when the declared size of its members is below MAX_XLSX_UNCOMPRESSED, parsed
by openpyxl with defusedxml (no entity expansion; formulas are never evaluated).
"""
from __future__ import annotations

import csv
import gzip
import io
import re
from pathlib import Path

import pandas as pd

from processor.safe_inputs import InputDir, UnsafeInput

MAX_TABLE_BYTES = 200 << 20             # a count table read into memory
MAX_SHEET_BYTES = 50 << 20              # a metadata sheet
MAX_XLSX_UNCOMPRESSED = 400 << 20       # the declared size of an .xlsx's members
TABLE_SUFFIXES = (".csv", ".tsv", ".txt", ".tab", ".counts", ".xlsx", ".csv.gz", ".tsv.gz", ".txt.gz", ".tab.gz")


def is_table(p: Path) -> bool:
    return p.name.lower().endswith(TABLE_SUFFIXES)


def snake(name) -> str:
    """'Subject experimental group' -> 'subject_experimental_group'."""
    s = re.sub(r"[^0-9a-zA-Z]+", "_", str(name).strip()).strip("_").lower()
    return s or "column"


def _text(inp: InputDir, path: Path, limit: int) -> str:
    raw = inp.read_bytes(path, limit)
    if path.name.lower().endswith(".gz"):
        with gzip.GzipFile(fileobj=io.BytesIO(raw)) as fh:
            raw = fh.read(limit + 1)
        if len(raw) > limit:
            raise UnsafeInput(f"{inp.rel(path)} is larger than {limit >> 20} MiB uncompressed")
    return raw.decode("utf-8-sig", errors="replace")


def read_raw(inp: InputDir, path: Path, limit: int = MAX_TABLE_BYTES) -> pd.DataFrame:
    """Every cell as a string or None (no header interpretation), comment lines ('#') dropped, ragged
    rows padded. XLSX: the first sheet; CSV/TSV/TXT: comma or tab, whichever the data lines have."""
    if path.name.lower().endswith(".xlsx"):
        if inp.xlsx_uncompressed_size(path) > MAX_XLSX_UNCOMPRESSED:
            raise UnsafeInput(f"{inp.rel(path)} unpacks to more than {MAX_XLSX_UNCOMPRESSED >> 20} MiB")
        data = inp.read_bytes(path, limit)
        return pd.read_excel(io.BytesIO(data), sheet_name=0, header=None, dtype=str, engine="openpyxl")
    lines = [ln for ln in _text(inp, path, limit).splitlines() if ln.strip() and not ln.startswith("#")]
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


def _blank(v) -> bool:
    return v is None or (isinstance(v, float) and pd.isna(v)) or str(v).strip() == ""


def with_header(raw: pd.DataFrame) -> pd.DataFrame:
    """Use the last non-numeric row before the data as the header: deposited spreadsheets often have a
    title or study row above the column names. Data starts at the first row that is mostly numbers.
    A header one cell shorter than the data (R's write.table with row names: no name above the gene
    ids) is shifted right by one, and the gene-id column is called gene_id."""
    raw = raw.dropna(axis=1, how="all").dropna(axis=0, how="all").reset_index(drop=True)
    start = 0
    for i in range(min(len(raw), 20)):
        vals = raw.iloc[i].dropna().astype(str).str.strip()
        vals = vals[vals != ""]
        if len(vals) and pd.to_numeric(vals, errors="coerce").notna().mean() >= 0.5:
            start = i
            break
    header_row = max(start - 1, 0)
    cells = list(raw.iloc[header_row])
    data = raw.iloc[header_row + 1:].reset_index(drop=True)
    filled = [i for i, v in enumerate(cells) if not _blank(v)]
    data_width = max((sum(not _blank(v) for v in r) for r in data.head(50).itertuples(index=False)), default=0)
    if filled and filled == list(range(len(filled))) and len(filled) == data_width - 1:
        cells = ["gene_id"] + cells[:-1]
    data.columns = [str(v).strip() if not _blank(v) else f"column_{j}" for j, v in enumerate(cells)]
    return data


def read_sheet(inp: InputDir, path: Path) -> pd.DataFrame:
    """A metadata sheet with its header row; column names in snake_case, cells stripped strings."""
    raw = read_raw(inp, path, MAX_SHEET_BYTES).dropna(axis=1, how="all").dropna(axis=0, how="all")
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
