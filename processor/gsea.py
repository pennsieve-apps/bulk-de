"""
Preranked gene set enrichment (GSEA, Subramanian et al. 2005) with gseapy, on
genes ranked by the Wald statistic of each contrast.

Gene sets are baked into the image (GENESETS_DIR, see genesets/build.py), all
from MSigDB 2026.2 (CC BY 4.0): the Hallmark collection (H for human, MH for
mouse) and GO biological process (C5 GO:BP / M5 GO:BP, from GO-basic 2026-07-26).
Genes are mapped to MSigDB symbols with MSigDB's own chip files: Ensembl gene ids
(version dropped) through the Ensembl chip, other ids and gene names through the
symbol-remapping chip (aliases and outdated symbols included). Rat genes map to
their mouse orthologs (MSigDB's Alliance of Genome Resources mapping) and use the
mouse sets. When several genes map to one symbol, the one with the largest |stat|
is kept.
"""
from __future__ import annotations

import gzip
import json
import logging
import os
import re
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger("bulk-de")

COLLECTIONS = ("hallmark", "go-bp")
SET_SPECIES = {"human": "human", "mouse": "mouse", "rat": "mouse"}
ENSEMBL = re.compile(r"^(ENS[A-Z]*G\d{11})(\.\d+)?$")
MIN_SIZE, MAX_SIZE, PERMUTATIONS = 15, 500, 1000


def genesets_dir() -> Path:
    return Path(os.environ.get("GENESETS_DIR", "/opt/genesets"))


@lru_cache(maxsize=1)
def release() -> dict:
    p = genesets_dir() / "genesets.json"
    return json.loads(p.read_text()) if p.is_file() else {}


@lru_cache(maxsize=8)
def gene_sets(collection: str, species: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    with gzip.open(genesets_dir() / f"{collection}.{species}.gmt.gz", "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                out[parts[0]] = parts[2:]
    return out


@lru_cache(maxsize=4)
def id_map(organism: str) -> dict[str, str]:
    """lower-case Ensembl id or symbol -> MSigDB symbol (of the gene-set species)."""
    out: dict[str, str] = {}
    with gzip.open(genesets_dir() / f"map.{organism}.tsv.gz", "rt") as fh:
        for line in fh:
            k, _, v = line.rstrip("\n").partition("\t")
            out.setdefault(k, v)
    return out


def guess_organism(gene_ids: pd.Index, gene_names: pd.Series) -> tuple[str, str]:
    ids = pd.Series(gene_ids.astype(str))
    for prefix, org in (("ENSMUSG", "mouse"), ("ENSRNOG", "rat"), ("ENSG", "human")):
        if ids.str.startswith(prefix).mean() > 0.5:
            return org, f"Ensembl ids ({prefix})"
    names = gene_names.where(gene_names != "", pd.Series(gene_ids, index=gene_names.index)).astype(str)
    letters = names[names.str.contains("[A-Za-z]{2}")]
    if len(letters) and (letters == letters.str.upper()).mean() > 0.8:
        return "human", "upper-case gene symbols"
    return "mouse", "capitalised gene symbols (mouse or rat; set ORGANISM=rat for rat)"


def map_genes(gene_ids: pd.Index, gene_names: pd.Series, organism: str) -> pd.Series:
    """gene_id -> MSigDB symbol ('' when unmapped)."""
    m = id_map(organism)
    out = []
    for gid, name in zip(gene_ids.astype(str), gene_names.reindex(gene_ids).fillna("").astype(str)):
        e = ENSEMBL.match(gid)
        sym = m.get(e.group(1).lower(), "") if e else ""
        if not sym:
            sym = m.get(name.lower(), "") if name else ""
        if not sym and not e and re.search("[A-Za-z]", gid):   # numeric ids (Entrez) are not symbols
            sym = m.get(gid.lower(), "")
        out.append(sym)
    return pd.Series(out, index=gene_ids)


def ranking(res: pd.DataFrame, symbols: pd.Series) -> pd.Series:
    r = pd.DataFrame({"symbol": symbols.reindex(res["gene_id"]).to_numpy(), "stat": res["stat"].to_numpy()})
    r = r[(r["symbol"] != "") & r["stat"].notna()]
    r = r.loc[r["stat"].abs().groupby(r["symbol"]).idxmax()]
    return r.set_index("symbol")["stat"].sort_values(ascending=False)


def run(rnk: pd.Series, collections: list[str], species: str, seed: int, threads: int) -> pd.DataFrame:
    """GSEApy's Rust core (gseapy.gse.prerank_rs, as gseapy.prerank calls it), read set by set: gseapy.prerank
    keeps every set's running-sum vector as Python lists, > 3 GB for GO:BP on 13 k genes."""
    from gseapy.gse import prerank_rs

    genes, values = rnk.index.tolist(), rnk.to_numpy(dtype=float).tolist()
    rows = []
    for coll in collections:
        res = prerank_rs(genes, values, gene_sets(coll, species), 1.0, MIN_SIZE, MAX_SIZE, PERMUTATIONS,
                         max(threads, 1), seed)
        if not res.summaries:
            log.warning("GSEA %s: no gene set has %d-%d of the ranked genes", coll, MIN_SIZE, MAX_SIZE)
        for gs in res.summaries:
            running = np.asarray(gs.run_es)
            hits = list(gs.hits)
            if gs.es >= 0:
                lead = [h for h in hits if h <= int(running.argmax())]
            else:
                lead = [h for h in hits if h >= int(running.argmin())][::-1]
            rows.append({"collection": coll, "term": gs.term, "size": len(hits), "es": gs.es, "nes": gs.nes,
                         "pvalue": gs.pval, "fdr": gs.fdr, "n_leading_edge": len(lead),
                         "leading_edge": ";".join(genes[h] for h in lead)})
    out = pd.DataFrame(rows, columns=["collection", "term", "size", "es", "nes", "pvalue", "fdr", "n_leading_edge",
                                      "leading_edge"]).astype({"es": float, "nes": float, "pvalue": float, "fdr": float})
    if out.empty:
        return out
    out = out.assign(_abs=-out["nes"].abs())
    return out.sort_values(["collection", "pvalue", "_abs"]).drop(columns="_abs").reset_index(drop=True)


def names_from_ensembl(gene_ids: pd.Index, gene_names: pd.Series, organism: str) -> pd.Series:
    """Fill empty gene names of Ensembl ids from the chip (human and mouse only: the rat chip gives mouse symbols)."""
    if organism not in ("human", "mouse") or not (genesets_dir() / f"map.{organism}.tsv.gz").is_file():
        return gene_names
    m = id_map(organism)
    filled = [n or (m.get(e.group(1).lower(), "") if (e := ENSEMBL.match(g)) else "")
              for g, n in zip(gene_ids.astype(str), gene_names.reindex(gene_ids).fillna("").astype(str))]
    return pd.Series(filled, index=gene_ids)
