"""
bulk-de: differential expression between groups of bulk RNA-seq samples.

Reads gene counts in any supported layout (processor/counts.py) and a sample
sheet (processor/metadata.py) from INPUT_DIR, fits DESeq2's model with PyDESeq2
for DESIGN, tests each contrast, ranks genes for preranked GSEA, and writes,
under OUTPUT_DIR/<LABEL>/:

  design.json              formula, contrasts, reference levels, samples used, gene filter
  samples.csv              per sample: design values, library size, genes detected, size factor, used
  pca.csv, pca.png         PCA of the VST of the 500 most variable genes
  <contrast>.csv           per gene: DESeq2 Wald results (+ shrunk log2 fold change)
  <contrast>_volcano.png, <contrast>_ma.png
  <contrast>_gsea.csv      preranked GSEA per gene-set collection
  de.json                  what was found in the input (layout, files, samples, genes, organism, mapping)
  summary.json             per contrast: genes up/down at FDR, top genes, top gene sets
  provenance.json          versions, parameters, inputs, timings, memory, citations
"""
from __future__ import annotations

import json
import logging
import math
import os
import platform
import re
import sys
import time
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("bulk-de")
log.setLevel(logging.INFO)   # the Lambda runtime sets the root logger to WARNING first
for noisy in ("matplotlib", "PIL", "fontTools", "gseapy"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "bulk-de"
REPO_ROOT = Path(__file__).resolve().parent.parent
LABEL_RE = re.compile(r"^[a-z0-9-]+$")
ORGANISMS = ("auto", "human", "mouse", "rat")
CONDITION_FALLBACK = ("subject_experimental_group", "sample_experimental_group")

DEFAULTS = {
    "LABEL": "",
    "DESIGN": "~ condition",
    "CONTRASTS": "",
    "REFERENCE_LEVELS": "",
    "SAMPLE_SHEET": "",
    "MIN_COUNT": "10",
    "MIN_SAMPLES": "0",
    "FDR": "0.05",
    "LFC_SHRINK": "true",
    "GENE_SETS": "hallmark,go-bp",
    "ORGANISM": "auto",
    "SEED": "42",
    "THREADS": "0",
}
PARAMETERS = tuple(DEFAULTS)


class ParameterError(ValueError):
    pass


def read_parameters(env=os.environ) -> dict:
    from processor import de as D
    from processor.gsea import COLLECTIONS

    raw = {k: (str(env.get(k)).strip() if env.get(k) not in (None, "") else v) for k, v in DEFAULTS.items()}
    p: dict = {k: raw[k] for k in ("DESIGN", "CONTRASTS", "REFERENCE_LEVELS", "SAMPLE_SHEET")}
    try:
        p["terms"] = D.parse_design(raw["DESIGN"])
        p["contrasts"] = D.parse_contrasts(raw["CONTRASTS"])
        p["references"] = D.parse_pairs(raw["REFERENCE_LEVELS"], "REFERENCE_LEVELS")
    except D.DesignError as exc:
        raise ParameterError(str(exc)) from exc
    if "/" in p["SAMPLE_SHEET"] or p["SAMPLE_SHEET"].startswith("."):
        raise ParameterError("SAMPLE_SHEET is a file name, not a path")
    label = raw["LABEL"] or re.sub(r"[^a-z0-9-]+", "-", p["terms"][-1].lower().replace("_", "-")).strip("-")
    if not LABEL_RE.match(label):
        raise ParameterError(f"LABEL must match [a-z0-9-]+, got {label!r}")
    p["LABEL"] = label
    try:
        for k in ("MIN_COUNT", "MIN_SAMPLES", "SEED", "THREADS"):
            p[k] = int(raw[k])
        p["FDR"] = float(raw["FDR"])
    except ValueError as exc:
        raise ParameterError(f"bad numeric parameter: {exc}") from exc
    if min(p["MIN_COUNT"], p["MIN_SAMPLES"], p["THREADS"]) < 0:
        raise ParameterError("MIN_COUNT, MIN_SAMPLES and THREADS must be >= 0")
    if not 0 < p["FDR"] < 1:
        raise ParameterError("FDR must be in (0, 1)")
    shrink = raw["LFC_SHRINK"].lower()
    if shrink not in ("true", "false"):
        raise ParameterError("LFC_SHRINK must be true or false")
    p["LFC_SHRINK"] = shrink == "true"
    sets = [s.strip().lower() for s in raw["GENE_SETS"].split(",") if s.strip()]
    if sets == ["none"]:
        sets = []
    bad = [s for s in sets if s not in COLLECTIONS]
    if bad:
        raise ParameterError(f"GENE_SETS: unknown {', '.join(bad)} (use hallmark, go-bp or none)")
    p["GENE_SETS"] = sets
    p["ORGANISM"] = raw["ORGANISM"].lower()
    if p["ORGANISM"] not in ORGANISMS:
        raise ParameterError(f"ORGANISM must be one of {', '.join(ORGANISMS)}")
    return p


# ---- helpers -----------------------------------------------------------------------------------

def _rel(path, root: Path):
    try:
        return str(Path(path).relative_to(root))
    except ValueError:
        return str(path)


def _clean_json(obj):
    import numpy as np

    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        obj = float(obj)
        return None if math.isnan(obj) or math.isinf(obj) else obj
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, dict):
        return {str(k): _clean_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean_json(v) for v in obj]
    return obj


def write_json(path: Path, data) -> None:
    with open(path, "w") as fh:
        json.dump(_clean_json(data), fh, indent=2)
        fh.write("\n")


def write_csv(df, path: Path) -> None:
    df.to_csv(path, index=False, float_format="%.6g", lineterminator="\n")


def tool_versions() -> dict:
    out = {"python": platform.python_version()}
    for pkg in ("pydeseq2", "gseapy", "numpy", "scipy", "pandas", "scikit-learn", "formulaic", "matplotlib",
                "openpyxl"):
        try:
            out[pkg] = version(pkg)
        except PackageNotFoundError:
            out[pkg] = None
    return out


def peak_rss_mb():
    try:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return round(peak / (1024 * 1024 if sys.platform == "darwin" else 1024), 1)
    except (ImportError, OSError):
        return None


def cgroup_peak_mb():
    try:
        return round(int(open("/sys/fs/cgroup/memory.peak").read()) / 1024 / 1024, 1)
    except (OSError, ValueError):
        return None


def app_version() -> str:
    try:
        return (REPO_ROOT / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def citations(keys: list[str]) -> list[dict]:
    refs = json.loads((Path(__file__).parent / "references.json").read_text(encoding="utf-8"))
    out = []
    for k in keys:
        r = refs[k]
        where = str(r["volume"]) + (f"({r['issue']})" if r.get("issue") else "") + (f":{r['pages']}" if r.get("pages") else "")
        c = {"doi": r["doi"], "text": f"{r['short_authors']}. {r['title']}. {r['journal_short']} {r['year']}"
                                      + (f";{where}" if where else "") + "."}
        if r.get("pmid"):
            c["pmid"] = r["pmid"]
        out.append(c)
    return out


def _top(res, fdr: float, up: bool, n: int = 10) -> list[dict]:
    sig = res[(res["padj"] < fdr) & ((res["log2_fold_change"] > 0) if up else (res["log2_fold_change"] < 0))]
    return [{"gene_id": r["gene_id"], "gene_name": r["gene_name"], "log2_fold_change": r["log2_fold_change"],
             "padj": r["padj"]} for _, r in sig.nsmallest(n, "padj").iterrows()]


# ---- main --------------------------------------------------------------------------------------

def main() -> int:
    started = datetime.now(timezone.utc)
    t0 = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 2
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    try:
        params = read_parameters()
    except ParameterError as exc:
        log.error("%s", exc)
        return 2
    n_cpus = int(os.environ.get("BULK_DE_THREADS") or params["THREADS"] or 1)

    import numpy as np
    import pandas as pd

    from processor import counts as C
    from processor import de as D
    from processor import figures as F
    from processor import gsea as G
    from processor import metadata as M
    from processor import tables as T

    np.random.seed(params["SEED"])
    files = T.walk(input_dir)
    port = input_dir / "counts"
    count_files = [f for f in files if port.is_dir() and T._inside(str(f), str(port))] or files
    try:
        counts = C.load(count_files, input_dir)
    except C.NoCounts as exc:
        log.error("%s under %s", exc, input_dir)
        return 1
    samples = list(counts.matrix.columns)
    log.info("counts: %s layout, %d genes x %d samples from %d file(s)", counts.layout, counts.matrix.shape[0],
             len(samples), len(counts.files))
    for n in counts.notes:
        log.info("counts: %s", n)
    if params["SAMPLE_SHEET"] and not any(f.name == params["SAMPLE_SHEET"] for f in files):
        log.error("SAMPLE_SHEET %s not found under INPUT_DIR", params["SAMPLE_SHEET"])
        return 2
    meta = M.find(files, set(counts.files), samples, params["SAMPLE_SHEET"])
    table = meta.table
    notes = counts.notes + meta.notes
    terms = params["terms"]
    if "condition" in terms and "condition" not in table.columns:
        alias = next((c for c in CONDITION_FALLBACK if c in table.columns), None)
        if alias:
            table = table.assign(condition=table[alias])
            notes.append(f"condition taken from {alias}")
    missing = [t for t in terms if t not in table.columns]
    if missing:
        cols = ", ".join(table.columns) or "none"
        log.error("DESIGN names %s, not a column of the sample sheet %s (columns: %s)", ", ".join(missing),
                  meta.sheet or "(none found)", cols)
        return 2 if meta.sheet else 1

    design_vals = table[terms].astype(str).apply(lambda c: c.str.strip())
    na = design_vals.isin(["", "NA", "NaN", "nan", "None"]).any(axis=1)
    statuses = {s: {"sample_id": s, "used": not na[s]} for s in samples}
    for s in samples:
        if na[s]:
            statuses[s]["reason"] = "no value for " + ", ".join(t for t in terms if design_vals.at[s, t] in
                                                                ("", "NA", "NaN", "nan", "None"))
            log.warning("%s: left out, %s", s, statuses[s]["reason"])
    used = [s for s in samples if not na[s]]
    if len(used) < 3:
        log.error("%d samples have every DESIGN value; need at least 3", len(used))
        return 1
    md = design_vals.loc[used]
    try:
        design = D.build(terms, md, params["contrasts"], params["references"])
    except D.DesignError as exc:
        log.error("%s", exc)
        return 2

    # organism, gene names
    organism, how = params["ORGANISM"], "ORGANISM"
    if organism == "auto":
        organism = M.organism(table.loc[used])
        how = "sample sheet species"
        if organism is None:
            organism, how = G.guess_organism(counts.matrix.index, counts.gene_names)
    gene_names = G.names_from_ensembl(counts.matrix.index, counts.gene_names, organism)
    log.info("organism %s (%s)", organism, how)

    out = output_dir / params["LABEL"]
    out.mkdir(parents=True, exist_ok=True)

    # filter and fit
    t_fit = time.monotonic()
    mat = counts.matrix[used].T                                    # samples x genes
    last = terms[-1]
    min_samples = params["MIN_SAMPLES"] or int(md[last].value_counts().min())
    filt = D.filter_genes(mat, params["MIN_COUNT"], min_samples)
    log.info("%d of %d genes have >= %d counts in >= %d samples", filt.shape[1], mat.shape[1], params["MIN_COUNT"],
             min_samples)
    if filt.shape[1] < 10:
        log.error("%d genes pass MIN_COUNT=%d in MIN_SAMPLES=%d samples", filt.shape[1], params["MIN_COUNT"],
                  min_samples)
        return 1
    fit = D.fit(filt, md, design, params["FDR"], params["LFC_SHRINK"], n_cpus)
    t_fit = time.monotonic() - t_fit
    log.info("PyDESeq2: %s, %d contrasts, %.1fs", design.formula, len(design.contrasts), t_fit)

    # samples.csv, pca
    lib = counts.matrix.sum(axis=0)
    detected = (counts.matrix > 0).sum(axis=0)
    srows = []
    for s in samples:
        srows.append({"sample_id": s, **{t: design_vals.at[s, t] for t in terms}, "library_size": int(lib[s]),
                      "n_genes_detected": int(detected[s]),
                      "size_factor": float(fit.size_factors[s]) if s in fit.size_factors.index else None,
                      "used": statuses[s]["used"]})
        statuses[s].update(library_size=int(lib[s]), n_genes_detected=int(detected[s]))
        log.info("%s: %d counts, %d genes detected, %s%s", s, lib[s], detected[s],
                 "used" if statuses[s]["used"] else "left out",
                 f", size factor {fit.size_factors[s]:.3f}" if s in fit.size_factors.index else "")
    write_csv(pd.DataFrame(srows), out / "samples.csv")
    scores, var = D.pca(fit.vst)
    write_csv(pd.concat([md, scores], axis=1).reset_index().rename(columns={"index": "sample_id"}), out / "pca.csv")
    figures = []
    try:
        figures.append(F.pca(scores, var, md, last, terms[0] if len(terms) > 1 else None,
                             f"PCA of VST counts, top {min(D.PCA_TOP_GENES, fit.vst.shape[1])} genes · "
                             f"{params['LABEL']}", out / "pca.png"))
    except Exception as exc:  # figures are a convenience; the tables stand without them
        log.exception("PCA figure failed: %s", exc)

    # contrasts
    symbols = None
    gsea_note = None
    if params["GENE_SETS"]:
        if not G.genesets_dir().is_dir():
            gsea_note = f"gene sets not found in {G.genesets_dir()}; GSEA skipped"
        else:
            symbols = G.map_genes(filt.columns, gene_names, organism)
            n_mapped = int((symbols != "").sum())
            log.info("GSEA: %d of %d tested genes map to MSigDB %s symbols", n_mapped, len(symbols),
                     G.SET_SPECIES[organism])
            if n_mapped < 0.3 * len(symbols):
                notes.append(f"only {n_mapped} of {len(symbols)} tested genes map to MSigDB symbols; is ORGANISM "
                             f"({organism}) right?")
    if gsea_note:
        notes.append(gsea_note)
        log.warning("%s", gsea_note)
    t_gsea = 0.0
    per_contrast = []
    for c in design.contrasts:
        r = fit.results[c.name]
        res = pd.DataFrame({"gene_id": r.index.astype(str), "gene_name": gene_names.reindex(r.index).fillna("").to_numpy()})
        for col in ("base_mean", "log2_fold_change", "lfc_se", "stat", "pvalue", "padj", "log2_fold_change_shrunk"):
            if col in r:
                res[col] = r[col].to_numpy()
        res = res.sort_values("pvalue", na_position="last").reset_index(drop=True)
        write_csv(res, out / f"{c.name}.csv")
        sig = res[res["padj"] < params["FDR"]]
        n_up, n_down = int((sig["log2_fold_change"] > 0).sum()), int((sig["log2_fold_change"] < 0).sum())
        title = f"{c.column}: {c.level} vs {c.reference}"
        try:
            figures.append(F.volcano(res, params["FDR"], f"{title} ({n_up} up, {n_down} down at FDR {params['FDR']:g})",
                                     out / f"{c.name}_volcano.png"))
            figures.append(F.ma(res, params["FDR"], title, out / f"{c.name}_ma.png"))
        except Exception as exc:  # figures are a convenience
            log.exception("figures for %s failed: %s", c.name, exc)
        entry = {"contrast": c.text, "name": c.name, "column": c.column, "level": c.level, "reference": c.reference,
                 "n_samples_level": int((md[c.column] == c.level).sum()),
                 "n_samples_reference": int((md[c.column] == c.reference).sum()),
                 "n_tested": int(res["padj"].notna().sum()), "n_significant": len(sig), "n_up": n_up,
                 "n_down": n_down, "lfc_shrunk": fit.shrunk[c.name],
                 "top_up": _top(res, params["FDR"], True), "top_down": _top(res, params["FDR"], False)}
        if symbols is not None:
            t = time.monotonic()
            rnk = G.ranking(res, symbols)
            gs = G.run(rnk, params["GENE_SETS"], G.SET_SPECIES[organism], params["SEED"], n_cpus)
            write_csv(gs, out / f"{c.name}_gsea.csv")
            t_gsea += time.monotonic() - t
            entry["gsea"] = {
                coll: {"n_sets": int((gs["collection"] == coll).sum()),
                       "n_significant": int(((gs["collection"] == coll) & (gs["fdr"] < params["FDR"])).sum()),
                       "top": [{"term": x["term"], "nes": x["nes"], "fdr": x["fdr"]}
                               for _, x in gs[(gs["collection"] == coll)].nsmallest(10, "fdr").iterrows()]}
                for coll in params["GENE_SETS"]}
            entry["gsea_ranked_genes"] = len(rnk)
            log.info("%s: GSEA on %d ranked genes, %s, %.1fs", c.name, len(rnk),
                     ", ".join(f"{k} {v['n_significant']}/{v['n_sets']} sets at FDR" for k, v in entry["gsea"].items()),
                     time.monotonic() - t)
        per_contrast.append(entry)
        log.info("%s: %d tested, %d up, %d down at FDR %g%s", c.name, entry["n_tested"], n_up, n_down, params["FDR"],
                 "" if fit.shrunk[c.name] or not params["LFC_SHRINK"] else
                 " (no shrunk LFC: the reference is not the factor's base level)")

    # json
    write_json(out / "design.json", {
        "formula": "~ " + " + ".join(terms), "terms": terms, "contrasts": [c.text for c in design.contrasts],
        "reference_levels": design.base, "levels": design.levels,
        "samples": [{"sample_id": s, **{t: md.at[s, t] for t in terms}} for s in used],
        "samples_left_out": [{"sample_id": s, "reason": statuses[s]["reason"]} for s in samples if not statuses[s]["used"]],
        "genes": {"input": int(mat.shape[1]), "tested": int(filt.shape[1]), "min_count": params["MIN_COUNT"],
                  "min_samples": min_samples},
    })
    write_json(out / "de.json", {
        "layout": counts.layout, "files": [_rel(f, input_dir) for f in counts.files],
        "n_genes": int(counts.matrix.shape[0]), "n_samples": len(samples), "samples": samples,
        "sample_sheet": meta.sheet, "sample_sheet_column": meta.match_column, "subjects_sheet": meta.subjects_sheet,
        "n_samples_matched": meta.n_matched, "organism": organism, "organism_from": how,
        "gene_sets": {"collections": params["GENE_SETS"], "release": G.release().get("source"),
                      "species": G.SET_SPECIES[organism],
                      "n_genes_mapped": int((symbols != "").sum()) if symbols is not None else None},
        "skipped": {**{_rel(k, input_dir): v for k, v in counts.skipped.items()},
                    **{_rel(k, input_dir): v for k, v in meta.skipped.items()}},
        "notes": notes,
    })
    write_json(out / "summary.json", {
        "label": params["LABEL"], "formula": "~ " + " + ".join(terms), "fdr": params["FDR"],
        "n_samples": len(used), "n_samples_left_out": len(samples) - len(used),
        "n_genes_tested": int(filt.shape[1]), "organism": organism, "layout": counts.layout,
        "pca": {"variance_explained": var[:5], "n_genes": min(D.PCA_TOP_GENES, int(filt.shape[1]))},
        "contrasts": per_contrast, "notes": notes, "figures": [p.name for p in figures],
    })
    keys = ["pydeseq2", "deseq2", "vst", "bh"] + (["apeglm"] if any(fit.shrunk.values()) else [])
    if symbols is not None:
        keys += ["gsea", "gseapy"]
        if "hallmark" in params["GENE_SETS"]:
            keys.append("msigdb_hallmark")
        if G.SET_SPECIES[organism] == "mouse":
            keys.append("msigdb_mouse")
        if "go-bp" in params["GENE_SETS"]:
            keys += ["go", "go_2023"]
    prov_params = {k: v for k, v in params.items() if k in PARAMETERS}
    prov_params["GENE_SETS"] = ",".join(params["GENE_SETS"]) or "none"
    write_json(out / "provenance.json", {
        "app": APP_NAME, "version": app_version(), "tools": tool_versions(), "parameters": prov_params,
        "resolved": {"design": design.formula, "tokens": design.tokens, "contrasts": [c.text for c in design.contrasts],
                     "organism": organism, "min_samples": min_samples,
                     "threads": {"n": n_cpus, "source": os.environ.get("BULK_DE_THREADS_SOURCE", "THREADS")}},
        "gene_sets": G.release() if symbols is not None else None,
        "inputs": {"counts": [_rel(f, input_dir) for f in counts.files], "sample_sheet": meta.sheet,
                   "subjects_sheet": meta.subjects_sheet},
        "started": started.isoformat(), "finished": datetime.now(timezone.utc).isoformat(),
        "timings_seconds": {"deseq2": round(t_fit, 2), "gsea": round(t_gsea, 2),
                            "total": round(time.monotonic() - t0, 2)},
        "peak_rss_mb": peak_rss_mb(), "cgroup_memory_peak_mb": cgroup_peak_mb(),
        "samples": list(statuses.values()), "notes": notes,
        "citations": citations(keys),
        "outputs": sorted(p.name for p in out.iterdir()) + ["provenance.json"],
    })
    log.info("done: %d samples, %d genes tested, %s, %.1fs, peak RSS %s MB -> %s", len(used), filt.shape[1],
             "; ".join(f"{e['name']} {e['n_up']} up / {e['n_down']} down" for e in per_contrast),
             time.monotonic() - t0, peak_rss_mb(), out)
    for n in notes:
        log.info("note: %s", n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
