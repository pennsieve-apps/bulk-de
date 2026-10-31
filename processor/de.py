"""
Design, contrasts and the DESeq2 model (PyDESeq2).

DESIGN is an additive formula of sample-sheet columns, '~ a + b + c'; every
term is a categorical factor (numbers are levels too). Each factor's base level
is its REFERENCE_LEVELS entry, else the reference of its first contrast, else
the first level alphabetically. Columns and levels are encoded as tokens (F0,
L00, ...) before they reach formulaic, so any level name works, and the base
level is L00.

PyDESeq2 (Muzellec et al. 2023) follows DESeq2 (Love et al. 2014): median-of-
ratios size factors, gene-wise dispersions shrunk towards a parametric trend,
a negative binomial GLM, Cook's distance outlier refitting (>= 7 replicates),
Wald tests, independent filtering and Benjamini-Hochberg. Shrunk log2 fold
changes use PyDESeq2's lfc_shrink, an apeglm-type Cauchy prior (Zhu et al.
2019); it can shrink a coefficient of the model, so a contrast against the
factor's base level gets one and any other contrast does not. The PCA is on the
variance-stabilising transformation (Anders & Huber 2010), blind to the design,
of the 500 most variable genes, as DESeq2's plotPCA.
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

TERM = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
PCA_TOP_GENES = 500


class DesignError(ValueError):
    pass


def parse_design(formula: str) -> list[str]:
    f = formula.strip()
    if not f.startswith("~"):
        raise DesignError(f"DESIGN must start with '~', got {formula!r}")
    terms = [t.strip() for t in f[1:].split("+")]
    for t in terms:
        if not TERM.match(t):
            raise DesignError(f"DESIGN term {t!r}: only column names joined by '+' (no interactions, "
                              f"transformations or intercept terms)")
    if len(set(terms)) != len(terms):
        raise DesignError("DESIGN names a column twice")
    return terms


def parse_pairs(text: str, what: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in filter(None, (s.strip() for s in text.split(","))):
        col, sep, level = item.partition(":")
        if not sep or not col.strip() or not level.strip():
            raise DesignError(f"{what}: {item!r} is not column:level")
        out[col.strip()] = level.strip()
    return out


@dataclass
class Contrast:
    column: str
    level: str
    reference: str

    @property
    def name(self) -> str:
        return re.sub(r"[^A-Za-z0-9._-]+", "_", f"{self.column}-{self.level}-vs-{self.reference}").strip("_")

    @property
    def text(self) -> str:
        return f"{self.column}:{self.level}:{self.reference}"


def parse_contrasts(text: str) -> list[Contrast]:
    out = []
    for item in filter(None, (s.strip() for s in text.split(","))):
        parts = [p.strip() for p in item.split(":")]
        if len(parts) != 3 or not all(parts):
            raise DesignError(f"CONTRASTS: {item!r} is not column:level:reference")
        if parts[1] == parts[2]:
            raise DesignError(f"CONTRASTS: {item!r} compares a level with itself")
        out.append(Contrast(*parts))
    return out


@dataclass
class Design:
    terms: list[str]
    levels: dict[str, list[str]]                 # base level first
    contrasts: list[Contrast]
    tokens: dict[str, str] = field(default_factory=dict)              # column -> F<i>
    level_tokens: dict[str, dict[str, str]] = field(default_factory=dict)

    @property
    def base(self) -> dict[str, str]:
        return {c: lv[0] for c, lv in self.levels.items()}

    def frame(self, meta: pd.DataFrame) -> pd.DataFrame:
        return pd.DataFrame({self.tokens[c]: meta[c].map(self.level_tokens[c]) for c in self.terms},
                            index=meta.index)

    @property
    def formula(self) -> str:
        return "~ " + " + ".join(self.tokens[c] for c in self.terms)


def build(terms: list[str], meta: pd.DataFrame, contrasts: list[Contrast], references: dict[str, str]) -> Design:
    for col in references:
        if col not in terms:
            raise DesignError(f"REFERENCE_LEVELS names {col!r}, which is not a DESIGN term")
    for c in contrasts:
        if c.column not in terms:
            raise DesignError(f"contrast {c.text}: {c.column!r} is not a DESIGN term")
    levels: dict[str, list[str]] = {}
    for col in terms:
        present = sorted(set(meta[col]))
        if len(present) < 2:
            raise DesignError(f"DESIGN term {col!r} has one level ({present[0] if present else 'none'}) "
                              f"among the samples used")
        base = references.get(col) or next((c.reference for c in contrasts if c.column == col), None) or present[0]
        if base not in present:
            raise DesignError(f"reference level {base!r} is not a level of {col} (levels: {', '.join(present)})")
        levels[col] = [base] + [lv for lv in present if lv != base]
    last = terms[-1]
    if not contrasts:
        contrasts = [Contrast(last, lv, levels[last][0]) for lv in levels[last][1:]]
    for c in contrasts:
        for lv in (c.level, c.reference):
            if lv not in levels[c.column]:
                raise DesignError(f"contrast {c.text}: {lv!r} is not a level of {c.column} "
                                  f"(levels: {', '.join(sorted(levels[c.column]))})")
    names = [c.name for c in contrasts]
    clash = sorted({n for n in names if names.count(n) > 1})
    if clash:
        raise DesignError(f"contrasts with the same file name after replacing special characters: {', '.join(clash)}")
    d = Design(terms, levels, contrasts)
    d.tokens = {c: f"F{i}" for i, c in enumerate(terms)}
    d.level_tokens = {c: {lv: f"L{j:02d}" for j, lv in enumerate(lvs)} for c, lvs in levels.items()}
    X = pd.concat([pd.Series(1.0, index=meta.index)] +
                  [pd.get_dummies(meta[c], dtype=float).drop(columns=levels[c][0]) for c in terms], axis=1)
    rank = np.linalg.matrix_rank(X.to_numpy())
    if rank < X.shape[1]:
        raise DesignError("the design is not of full rank: two terms are confounded (e.g. every tumour sample "
                          "is also female); drop one term")
    if X.shape[0] - rank < 1:
        raise DesignError(f"{X.shape[0]} samples for {rank} coefficients leave no residual degrees of freedom")
    return d


def filter_genes(counts: pd.DataFrame, min_count: int, min_samples: int) -> pd.DataFrame:
    """counts: samples x genes. Keep genes with >= min_count counts in >= min_samples samples."""
    keep = ((counts >= min_count).sum(axis=0) >= max(1, min_samples)) & (counts.sum(axis=0) > 0)
    return counts.loc[:, keep]


@dataclass
class Fit:
    results: dict[str, pd.DataFrame]             # contrast name -> gene table
    shrunk: dict[str, bool]
    size_factors: pd.Series
    vst: pd.DataFrame                            # samples x genes


def fit(counts: pd.DataFrame, meta: pd.DataFrame, design: Design, fdr: float, shrink: bool, n_cpus: int) -> Fit:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    md = design.frame(meta)
    results, shrunk = {}, {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dds = DeseqDataSet(counts=counts, metadata=md, design=design.formula, refit_cooks=True, n_cpus=n_cpus,
                           quiet=True)
        dds.deseq2()
        for c in design.contrasts:
            f, lt = design.tokens[c.column], design.level_tokens[c.column]
            ds = DeseqStats(dds, contrast=[f, lt[c.level], lt[c.reference]], alpha=fdr, n_cpus=n_cpus, quiet=True)
            ds.summary()
            r = ds.results_df.copy()
            can_shrink = shrink and c.reference == design.base[c.column]
            if can_shrink:
                ds.lfc_shrink(coeff=f"{f}[T.{lt[c.level]}]")
                r["log2_fold_change_shrunk"] = ds.results_df["log2FoldChange"]
            results[c.name] = r.rename(columns={"baseMean": "base_mean", "log2FoldChange": "log2_fold_change",
                                                "lfcSE": "lfc_se"})
            shrunk[c.name] = can_shrink
        size_factors = pd.Series(np.asarray(dds.obs["size_factors"]), index=counts.index)
        dds.vst(use_design=False)
        vst = pd.DataFrame(np.asarray(dds.layers["vst_counts"]), index=counts.index, columns=counts.columns)
    return Fit(results, shrunk, size_factors, vst)


def pca(vst: pd.DataFrame, n_top: int = PCA_TOP_GENES) -> tuple[pd.DataFrame, list[float]]:
    top = vst.var(axis=0).sort_values(ascending=False).index[:n_top]
    x = vst[top].to_numpy()
    x = x - x.mean(axis=0)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    k = min(len(s), max(x.shape[0] - 1, 1), 10)
    scores = u[:, :k] * s[:k]
    var = (s ** 2) / max((s ** 2).sum(), 1e-12)
    return (pd.DataFrame(scores, index=vst.index, columns=[f"PC{i + 1}" for i in range(k)]),
            [float(v) for v in var[:k]])
