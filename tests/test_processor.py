"""End-to-end tests on the synthetic fixtures (tests/make_fixture.py), no Docker."""
from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from processor import counts as C
from tests.make_fixture import DOWN, GRADE, N_GENES, SEX, UP, human_names, make_fixture, mouse_names

DOI = {k: v["doi"] for k, v in json.loads((Path(__file__).parent.parent / "processor" / "references.json")
                                         .read_text(encoding="utf-8")).items()}
PARAMS = ("LABEL", "DESIGN", "CONTRASTS", "REFERENCE_LEVELS", "SAMPLE_SHEET", "MIN_COUNT", "MIN_SAMPLES", "FDR",
          "LFC_SHRINK", "GENE_SETS", "ORGANISM", "SEED", "THREADS")


@pytest.fixture(scope="session")
def fx(tmp_path_factory):
    return make_fixture(tmp_path_factory.mktemp("fixture"))


def run(input_dir: Path, output_dir: Path, monkeypatch, genesets: Path, **params) -> int:
    for k in PARAMS:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    monkeypatch.setenv("BULK_DE_THREADS", "1")
    monkeypatch.setenv("GENESETS_DIR", str(genesets))
    for k, v in params.items():
        monkeypatch.setenv(k, str(v))
    from processor import gsea
    from processor.main import main

    for f in (gsea.release, gsea.gene_sets, gsea.id_map):
        f.cache_clear()
    return main()


def test_rejoin_tables(fx, tmp_path, monkeypatch):
    """Deposited genes x samples CSV + separate sample xlsx, design with a donor block."""
    assert run(fx / "rejoin", tmp_path, monkeypatch, fx / "genesets", DESIGN="~ donor + tissue") == 0
    out = tmp_path / "tissue"
    name = "tissue-meniscus-vs-cartilage"
    for f in ("design.json", "samples.csv", "pca.csv", "pca.png", f"{name}.csv", f"{name}_volcano.png",
              f"{name}_ma.png", f"{name}_gsea.csv", "de.json", "summary.json", "provenance.json"):
        assert (out / f).is_file(), f
    res = pd.read_csv(out / f"{name}.csv", keep_default_na=False, na_values=[""])
    assert list(res.columns) == ["gene_id", "gene_name", "base_mean", "log2_fold_change", "lfc_se", "stat", "pvalue",
                                 "padj", "log2_fold_change_shrunk"]
    sig = res[res["padj"] < 0.05].set_index("gene_name")
    names = human_names()
    assert {names[i] for i in UP} <= set(sig.index) and {names[i] for i in DOWN} <= set(sig.index)
    assert (sig.loc[[names[i] for i in UP], "log2_fold_change"] > 1).all()
    assert len(set(sig.index) - {names[i] for i in UP + DOWN}) <= 3
    # shrinkage pulls estimates towards zero, never past it
    both = res.dropna(subset=["log2_fold_change_shrunk"])
    assert (both["log2_fold_change_shrunk"].abs() <= both["log2_fold_change"].abs() + 1e-6).all()
    assert res["gene_id"].iloc[0].startswith("ENSG")

    design = json.loads((out / "design.json").read_text())
    assert design["formula"] == "~ donor + tissue" and design["contrasts"] == ["tissue:meniscus:cartilage"]
    assert len(design["samples"]) == 12 and set(design["samples"][0]) == {"sample_id", "donor", "tissue"}
    de = json.loads((out / "de.json").read_text())
    assert de["layout"] == "matrix" and de["sample_sheet"] == "sample_info.xlsx" and de["organism"] == "human"
    assert de["n_genes"] == N_GENES and de["gene_sets"]["n_genes_mapped"] > 0
    gs = pd.read_csv(out / f"{name}_gsea.csv")
    top = gs[gs["collection"] == "hallmark"].iloc[0]
    assert top["term"] == "SET_UP" and top["nes"] > 0 and top["fdr"] < 0.05
    pca = pd.read_csv(out / "pca.csv")
    assert {"sample_id", "donor", "tissue", "PC1", "PC2"} <= set(pca.columns)
    # PC1 separates the tissues
    assert pca.groupby("tissue")["PC1"].mean().diff().abs().iloc[-1] > pca["PC1"].std() * 0.5

    summary = json.loads((out / "summary.json").read_text())
    c = summary["contrasts"][0]
    assert c["n_up"] >= len(UP) and c["n_down"] >= len(DOWN) and c["lfc_shrunk"]
    assert c["top_up"][0]["gene_name"] in {names[i] for i in UP}
    prov = json.loads((out / "provenance.json").read_text())
    dois = {x["doi"] for x in prov["citations"]}
    assert {DOI[k] for k in ("pydeseq2", "deseq2", "apeglm", "gsea", "gseapy", "msigdb_hallmark", "go")} <= dois
    assert DOI["msigdb_mouse"] not in dois
    assert all(set(x) <= {"doi", "pmid", "text"} and x["doi"] and x["text"] for x in prov["citations"])
    for k in ("app", "version", "tools", "parameters", "inputs", "started", "finished", "samples", "peak_rss_mb"):
        assert k in prov
    assert prov["app"] == "bulk-de" and prov["version"] == "0.1.0" and prov["tools"]["pydeseq2"]
    # clinical columns of the sheet (Age) never reach the outputs
    assert "age" not in pd.read_csv(out / "samples.csv").columns


def test_rejoin_second_factor_and_reference(fx, tmp_path, monkeypatch):
    assert run(fx / "rejoin", tmp_path, monkeypatch, fx / "genesets", DESIGN="~ tissue + oa_grade",
               CONTRASTS="oa_grade:high:low,tissue:cartilage:meniscus", LABEL="grade", GENE_SETS="none") == 0
    out = tmp_path / "grade"
    g = pd.read_csv(out / "oa_grade-high-vs-low.csv", keep_default_na=False, na_values=[""])
    sig = set(g.loc[g["padj"] < 0.05, "gene_name"])
    assert {human_names()[i] for i in GRADE} <= sig
    t = pd.read_csv(out / "tissue-cartilage-vs-meniscus.csv", keep_default_na=False, na_values=[""])
    up = t.set_index("gene_name").loc[[human_names()[i] for i in UP]]
    assert (up["log2_fold_change"] < 0).all()                          # reversed: cartilage vs meniscus
    assert "log2_fold_change_shrunk" in t.columns                      # meniscus became the base level
    assert not list(out.glob("*_gsea.csv"))
    prov = json.loads((out / "provenance.json").read_text())
    assert DOI["gsea"] not in {x["doi"] for x in prov["citations"]}


def test_sparc_pools(fx, tmp_path, monkeypatch):
    """Per-library htseq files in pool- folders, metadata from SPARC samples.xlsx + subjects.xlsx."""
    # the experimental group already holds the sex: the design is confounded -> bad parameters
    assert run(fx / "sparc", tmp_path, monkeypatch, fx / "genesets",
               DESIGN="~ sex + subject_experimental_group") == 2
    assert run(fx / "sparc", tmp_path / "o", monkeypatch, fx / "genesets", DESIGN="~ sex") == 0
    out = tmp_path / "o" / "sex"
    de = json.loads((out / "de.json").read_text())
    assert de["layout"] == "two-column" and de["n_samples"] == 8 and de["sample_sheet"] == "samples.xlsx"
    assert de["sample_sheet_column"] == "pool_id" and de["subjects_sheet"] == "subjects.xlsx"
    assert de["organism"] == "mouse" and de["organism_from"] == "sample sheet species"
    res = pd.read_csv(out / "sex-male-vs-female.csv", keep_default_na=False, na_values=[""])
    sig = set(res.loc[res["padj"] < 0.05, "gene_id"])
    assert {mouse_names()[i] for i in SEX} <= sig
    assert not res["gene_id"].str.startswith("__").any()
    dois = {x["doi"] for x in json.loads((out / "provenance.json").read_text())["citations"]}
    assert DOI["msigdb_mouse"] in dois


def test_sparc_condition_from_group(fx, tmp_path, monkeypatch):
    """Default DESIGN ~ condition falls back to subject_experimental_group; every level vs the first."""
    assert run(fx / "sparc", tmp_path, monkeypatch, fx / "genesets", GENE_SETS="none") == 0
    out = tmp_path / "condition"
    design = json.loads((out / "design.json").read_text())
    assert design["contrasts"] == ["condition:Female-Tumor:Female-Control", "condition:Male-Control:Female-Control",
                                   "condition:Male-Tumor:Female-Control"]


def test_no_gene_set_matches(fx, tmp_path, monkeypatch):
    """No gene set overlaps the ranked genes: GSEA writes an empty table and the run succeeds."""
    sets = tmp_path / "sets"
    shutil.copytree(fx / "genesets", sets)
    (sets / "hallmark.human.gmt.gz").write_bytes(gzip.compress(b"SET_X\tx\tNOPE1\tNOPE2\n"))
    assert run(fx / "rejoin", tmp_path / "out", monkeypatch, sets, DESIGN="~ tissue", GENE_SETS="hallmark") == 0
    gs = pd.read_csv(tmp_path / "out" / "tissue" / "tissue-meniscus-vs-cartilage_gsea.csv")
    assert gs.empty and "nes" in gs.columns
    summary = json.loads((tmp_path / "out" / "tissue" / "summary.json").read_text())
    assert summary["contrasts"][0]["gsea"]["hallmark"]["n_sets"] == 0


def test_sample_sheet_csv(fx, tmp_path, monkeypatch):
    """A user sample sheet next to the counts; one sample without a value is left out."""
    inp = tmp_path / "in"
    shutil.copytree(fx / "sparc" / "primary", inp / "primary")
    rows = [{"library": p.name.removeprefix("pool-"), "sex": p.name[5], "condition": p.name[6]}
            for p in sorted((inp / "primary").iterdir())]
    rows[0]["condition"] = "NA"
    pd.DataFrame(rows).to_csv(inp / "design.csv", index=False)
    assert run(inp, tmp_path / "out", monkeypatch, fx / "genesets", DESIGN="~ sex + condition",
               CONTRASTS="condition:T:N,sex:M:F", SAMPLE_SHEET="design.csv", ORGANISM="mouse") == 0
    out = tmp_path / "out" / "condition"
    design = json.loads((out / "design.json").read_text())
    assert len(design["samples"]) == 7 and design["samples_left_out"][0]["sample_id"] == "pool-FN1"
    t = pd.read_csv(out / "condition-T-vs-N.csv", keep_default_na=False, na_values=[""])
    assert {mouse_names()[i] for i in UP} <= set(t.loc[t["padj"] < 0.05, "gene_id"])
    s = pd.read_csv(out / "samples.csv")
    assert s["used"].sum() == 7 and (s["library_size"] > 0).all()


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def test_other_layouts(tmp_path):
    genes = [f"ENSMUSG{i:011d}" for i in range(120)]
    rng = np.random.default_rng(0)
    # featureCounts
    lines = ["# Program:featureCounts v2.0.6; Command:...", "Geneid\tChr\tStart\tEnd\tStrand\tLength\t"
             "/x/a.Aligned.sortedByCoord.out.bam\t/x/b.bam"]
    lines += [f"{g}\tchr1\t1\t100\t+\t100\t{rng.integers(0, 50)}\t{rng.integers(0, 50)}" for g in genes]
    f = _write(tmp_path / "fc" / "counts.txt", "\n".join(lines) + "\n")
    c = C.load([f], tmp_path / "fc")
    assert c.layout == "featureCounts" and list(c.matrix.columns) == ["a", "b"]
    assert c.matrix.shape == (120, 2)
    # STAR, reverse-stranded
    files = []
    for s in ("s1", "s2"):
        body = "N_unmapped\t5\t5\t5\nN_multimapping\t9\t9\t9\nN_noFeature\t1\t50\t2\nN_ambiguous\t1\t1\t1\n"
        body += "".join(f"{g}\t{10 + i}\t{i % 3}\t{10 + i}\n" for i, g in enumerate(genes))
        files.append(_write(tmp_path / "star" / s / f"{s}ReadsPerGene.out.tab", body))
    c = C.load(files, tmp_path / "star")
    assert c.layout == "STAR" and list(c.matrix.columns) == ["s1", "s2"] and c.matrix.iloc[5, 0] == 15
    assert any("reverse-stranded" in n for n in c.notes)
    # RSEM, rounded
    files = [_write(tmp_path / "rsem" / f"{s}.genes.results",
                    "gene_id\ttranscript_id(s)\tlength\teffective_length\texpected_count\tTPM\tFPKM\n"
                    + "".join(f"{g}\tT\t100\t80\t{i + 0.4}\t1\t1\n" for i, g in enumerate(genes)))
             for s in ("r1", "r2", "r3")]
    c = C.load(files, tmp_path / "rsem")
    assert c.layout == "RSEM" and c.matrix.shape == (120, 3) and c.matrix.iloc[3, 0] == 3
    # rnaseq-quant gene_counts.csv wins over anything else, gene_tpm.csv is not counts
    q = pd.DataFrame({"gene_id": genes, "gene_name": [f"G{i}" for i in range(120)], "a": 5, "b": 6})
    q.to_csv(tmp_path / "rsem" / "gene_counts.csv", index=False)
    q.assign(a=0.5).to_csv(tmp_path / "rsem" / "gene_tpm.csv", index=False)
    c = C.load(sorted((tmp_path / "rsem").iterdir()), tmp_path / "rsem")
    assert c.layout == "rnaseq-quant" and c.gene_names.iloc[1] == "G1"


def test_negative_values_rejected(tmp_path):
    rows = "\n".join(f"g{i},{i},{-1 if i == 5 else i}" for i in range(150))
    f = _write(tmp_path / "m.csv", "gene,a,b\n" + rows + "\n")
    with pytest.raises(C.NoCounts):
        C.load([f], tmp_path)


def test_symlink_escape_ignored(fx, tmp_path, monkeypatch):
    inp = tmp_path / "in"
    shutil.copytree(fx / "rejoin", inp)
    outside = tmp_path / "secret"
    outside.mkdir()
    (outside / "other.csv").write_text("x,y\n1,2\n")
    (inp / "counts" / "link").symlink_to(outside)
    from processor.tables import walk

    names = [p.name for p in walk(inp)]
    assert "other.csv" not in names and "knee_spine_counts.csv" in names


def test_bad_parameters(fx, tmp_path, monkeypatch):
    g = fx / "genesets"
    r = fx / "rejoin"
    assert run(r, tmp_path, monkeypatch, g, DESIGN="tissue") == 2
    assert run(r, tmp_path, monkeypatch, g, DESIGN="~ tissue * donor") == 2
    assert run(r, tmp_path, monkeypatch, g, DESIGN="~ nope") == 2
    assert run(r, tmp_path, monkeypatch, g, DESIGN="~ tissue", CONTRASTS="tissue:bone:cartilage") == 2
    assert run(r, tmp_path, monkeypatch, g, DESIGN="~ tissue", CONTRASTS="donor:donor1:donor2") == 2
    assert run(r, tmp_path, monkeypatch, g, DESIGN="~ tissue", REFERENCE_LEVELS="tissue:bone") == 2
    assert run(r, tmp_path, monkeypatch, g, DESIGN="~ tissue", GENE_SETS="kegg") == 2
    assert run(r, tmp_path, monkeypatch, g, DESIGN="~ tissue", FDR="1.5") == 2
    assert run(r, tmp_path, monkeypatch, g, DESIGN="~ tissue", LABEL="Bad_Label") == 2
    assert run(r, tmp_path, monkeypatch, g, DESIGN="~ tissue", SAMPLE_SHEET="../x.csv") == 2
    assert run(r, tmp_path, monkeypatch, g, DESIGN="~ tissue", LFC_SHRINK="maybe") == 2


def test_missing_inputs(fx, tmp_path, monkeypatch):
    (tmp_path / "in").mkdir()
    assert run(tmp_path / "in", tmp_path / "out", monkeypatch, fx / "genesets") == 1


def test_handler(fx, tmp_path, monkeypatch):
    for k in PARAMS:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("BULK_DE_THREADS", "1")
    monkeypatch.setenv("GENESETS_DIR", str(fx / "genesets"))
    monkeypatch.setenv("FDR", "0.5")          # must be cleared by the handler
    from processor.handler import handler

    r = handler({"inputDir": str(fx / "rejoin"), "outputDir": str(tmp_path), "DESIGN": "~ tissue", "LABEL": "lam",
                 "LFC_SHRINK": False, "GENE_SETS": "none"}, None)
    assert r["status"] == "ok"
    summary = json.loads((tmp_path / "lam" / "summary.json").read_text())
    assert summary["fdr"] == 0.05 and not summary["contrasts"][0]["lfc_shrunk"]
    with pytest.raises(RuntimeError):
        handler({"inputDir": str(fx / "rejoin"), "outputDir": str(tmp_path), "DESIGN": "~ nope"}, None)

