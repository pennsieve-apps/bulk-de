"""
Synthetic fixtures with planted effects, plus a tiny gene-set folder (the image's
MSigDB files are not in the repository).

  rejoin/   the shape of RE-JOIN's deposited knee/spine tables: one genes x samples
            CSV (versioned Ensembl ids + gene names, a title row above the header) and
            a separate sample sheet sample_info.xlsx (Sample ID, Tissue, OA grade, Donor,
            Age). 12 samples: 6 donors x 2 tissues (cartilage, meniscus); OA grade high in
            donors 1-3. Planted: UP genes 4x higher in meniscus, DOWN genes 4x lower;
            GRADE genes 3x higher with high OA grade.
  sparc/    a SPARC-layout dataset like Discover 427: per-library htseq-style files
            primary/pool-<id>/<id>_Genes_ReadCount.txt (gene symbols, __ summary rows),
            samples.xlsx (two animals per pool: sample id, subject id, pool id) and
            subjects.xlsx (sex, subject experimental group, species). 8 pools: sex F/M x
            Control/Tumor x 2. Planted: UP genes 4x higher in Tumor, SEX genes 8x in male.
  genesets/ hallmark.{human,mouse}.gmt.gz with SET_UP (the planted UP genes + 10 others), SET_DOWN and
            two random sets; map.{human,mouse,rat}.tsv.gz; genesets.json.

    python -m tests.make_fixture data/fixture
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

N_GENES = 400
UP = list(range(0, 20))
DOWN = list(range(20, 32))
GRADE = list(range(40, 50))
SEX = list(range(60, 66))
BASE = np.exp(np.random.default_rng(0).normal(np.log(200.0), 1.2, N_GENES))


def human_names() -> list[str]:
    return [f"GENE{i:03d}" for i in range(N_GENES)]


def mouse_names() -> list[str]:
    return [f"Gene{i:03d}" for i in range(N_GENES)]


def ensembl_ids() -> list[str]:
    return [f"ENSG{i:011d}.{1 + i % 3}" for i in range(N_GENES)]


def _counts(mu: np.ndarray, rng: np.random.Generator, disp: float = 0.05) -> np.ndarray:
    n = 1.0 / disp
    return rng.negative_binomial(n, n / (n + mu))


def make_rejoin(root: Path, rng: np.random.Generator) -> None:
    rows, cols = [], {}
    for d in range(1, 7):
        for tissue in ("cartilage", "meniscus"):
            sid = f"D{d}_{tissue[:3].upper()}"
            grade = "high" if d <= 3 else "low"
            mu = BASE * np.exp(rng.normal(0, 0.15, N_GENES)) * np.exp(rng.normal(0, 0.2))
            if tissue == "meniscus":
                mu[UP] *= 4.0
                mu[DOWN] /= 4.0
            if grade == "high":
                mu[GRADE] *= 3.0
            cols[sid] = _counts(mu, rng)
            rows.append({"Sample ID": sid, "Tissue": tissue, "OA grade": grade, "Donor": f"donor{d}",
                         "Age": 50 + 3 * d})
    (root / "counts").mkdir(parents=True)
    (root / "metadata").mkdir(parents=True)
    df = pd.DataFrame(cols)
    df.insert(0, "Gene name", human_names())
    df.insert(0, "Ensembl ID", ensembl_ids())
    with open(root / "counts" / "knee_spine_counts.csv", "w") as fh:
        fh.write("RE-JOIN knee and spine tissues,,raw counts\n")
        df.to_csv(fh, index=False)
    pd.DataFrame(rows).to_excel(root / "metadata" / "sample_info.xlsx", index=False)


def make_sparc(root: Path, rng: np.random.Generator) -> None:
    samples, subjects = [], []
    n = 0
    for sex in ("F", "M"):
        for group, code in (("Control", "N"), ("Tumor", "T")):
            for rep in (1, 2):
                pool = f"{sex}{code}{rep}"
                mu = BASE * np.exp(rng.normal(0, 0.15, N_GENES)) * np.exp(rng.normal(0, 0.2))
                if group == "Tumor":
                    mu[UP] *= 4.0
                if sex == "M":
                    mu[SEX] *= 8.0
                c = _counts(mu, rng)
                d = root / "primary" / f"pool-{pool}"
                d.mkdir(parents=True)
                lines = [f"{g}\t{v}" for g, v in zip(mouse_names(), c)]
                lines += ["__no_feature\t1000", "__ambiguous\t20", "__alignment_not_unique\t5000"]
                (d / f"{pool}_Genes_ReadCount.txt").write_text("\n".join(lines) + "\n")
                for _ in range(2):
                    n += 1
                    samples.append({"sample id": f"sam-M{n}-{pool}", "subject id": f"sub-M{n}",
                                    "pool id": f"pool-{pool}", "sample type": "tissue"})
                    subjects.append({"subject id": f"sub-M{n}", "subject experimental group":
                                     f"{'Male' if sex == 'M' else 'Female'}-{group}",
                                     "sex": "male" if sex == "M" else "female", "species": "Mus musculus",
                                     "age": f"{6 + n % 3} weeks"})
    pd.DataFrame(samples).to_excel(root / "samples.xlsx", index=False)
    pd.DataFrame(subjects).to_excel(root / "subjects.xlsx", index=False)


def make_genesets(root: Path, rng: np.random.Generator) -> None:
    root.mkdir(parents=True)
    for species, names in (("human", human_names()), ("mouse", mouse_names())):
        others = [i for i in range(100, N_GENES)]
        sets = {
            "SET_UP": [names[i] for i in UP] + [names[i] for i in others[:10]],
            "SET_DOWN": [names[i] for i in DOWN] + [names[i] for i in others[10:15]],
            "SET_RANDOM_A": [names[i] for i in rng.choice(others, 30, replace=False)],
            "SET_RANDOM_B": [names[i] for i in rng.choice(others, 40, replace=False)],
        }
        text = "".join(f"{k}\tx\t" + "\t".join(v) + "\n" for k, v in sets.items())
        for coll in ("hallmark", "go-bp"):
            (root / f"{coll}.{species}.gmt.gz").write_bytes(gzip.compress(text.encode()))
    maps = {"human": [(e.split(".")[0].lower(), n) for e, n in zip(ensembl_ids(), human_names())]
            + [(n.lower(), n) for n in human_names()],
            "mouse": [(n.lower(), n) for n in mouse_names()],
            "rat": [(n.lower(), n) for n in mouse_names()]}
    for org, pairs in maps.items():
        (root / f"map.{org}.tsv.gz").write_bytes(gzip.compress("".join(f"{k}\t{v}\n" for k, v in pairs).encode()))
    (root / "genesets.json").write_text(json.dumps({"source": "synthetic test sets", "licence": "CC0"}) + "\n")


def make_fixture(root: Path) -> Path:
    rng = np.random.default_rng(1)
    make_rejoin(root / "rejoin", rng)
    make_sparc(root / "sparc", rng)
    make_genesets(root / "genesets", rng)
    return root


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    a = ap.parse_args()
    make_fixture(Path(a.out))
    print(f"fixture written to {a.out}")
