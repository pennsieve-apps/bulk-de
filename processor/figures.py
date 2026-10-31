"""Figures: PCA, volcano and MA plots (PNG, 150 dpi)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

DPI = 150
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]
UP, DOWN, NS = "#e34948", "#2a78d6", "#b9b8b2"
TEXT, TEXT2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": TEXT2, "axes.labelcolor": TEXT,
                         "xtick.color": TEXT2, "ytick.color": TEXT2, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.facecolor": "white"})
    return plt


def pca(scores: pd.DataFrame, var: list[float], meta: pd.DataFrame, colour: str, shape: str | None, title: str,
        path: Path) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    colours = sorted(meta[colour].unique())
    shapes = sorted(meta[shape].unique()) if shape else [None]
    y = scores["PC2"] if "PC2" in scores else pd.Series(0.0, index=scores.index)
    for i, c in enumerate(colours):
        for j, s in enumerate(shapes):
            sel = (meta[colour] == c) & ((meta[shape] == s) if shape else True)
            ids = meta.index[sel]
            if len(ids) == 0:
                continue
            ax.scatter(scores.loc[ids, "PC1"], y.loc[ids], s=46, color=SERIES[i % len(SERIES)],
                       marker=MARKERS[j % len(MARKERS)], edgecolor="white", linewidth=0.8,
                       label=f"{c}" + (f" · {s}" if shape else ""), zorder=3)
    if len(scores) <= 30:
        for sid in scores.index:
            ax.annotate(str(sid), (scores.at[sid, "PC1"], y.at[sid]), fontsize=6, color=TEXT2, xytext=(3, 3),
                        textcoords="offset points")
    ax.set_xlabel(f"PC1 ({var[0] * 100:.0f} % of variance)")
    ax.set_ylabel(f"PC2 ({var[1] * 100:.0f} % of variance)" if len(var) > 1 else "PC2")
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_title(title, fontsize=9, color=TEXT)
    ax.legend(frameon=False, fontsize=7, title=colour + (f" · {shape}" if shape else ""), title_fontsize=7,
              loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path


def _label(r: pd.DataFrame) -> pd.Series:
    return r["gene_name"].where(r["gene_name"].astype(str) != "", r["gene_id"]).astype(str)


def volcano(res: pd.DataFrame, fdr: float, title: str, path: Path, n_labels: int = 12) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    r = res.dropna(subset=["pvalue", "log2_fold_change"]).copy()
    r["nlp"] = -np.log10(np.clip(r["pvalue"].to_numpy(dtype=float), 1e-300, 1))
    sig = r["padj"].fillna(1) < fdr
    up, down = sig & (r["log2_fold_change"] > 0), sig & (r["log2_fold_change"] < 0)
    ax.scatter(r.loc[~sig, "log2_fold_change"], r.loc[~sig, "nlp"], s=5, color=NS, linewidth=0, label="not significant",
               rasterized=True)
    ax.scatter(r.loc[up, "log2_fold_change"], r.loc[up, "nlp"], s=10, color=UP, linewidth=0, label=f"up ({int(up.sum())})")
    ax.scatter(r.loc[down, "log2_fold_change"], r.loc[down, "nlp"], s=10, color=DOWN, linewidth=0,
               label=f"down ({int(down.sum())})")
    top = r[sig].nsmallest(n_labels, "pvalue")
    for (_, g), name in zip(top.iterrows(), _label(top)):
        ax.annotate(name, (g["log2_fold_change"], g["nlp"]), fontsize=7, color=TEXT, xytext=(3, 2),
                    textcoords="offset points")
    ax.axvline(0, color=GRID, linewidth=0.8, zorder=0)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 p")
    ax.set_title(title, fontsize=9, color=TEXT)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path


def ma(res: pd.DataFrame, fdr: float, title: str, path: Path) -> Path:
    plt = _plt()
    col = "log2_fold_change_shrunk" if "log2_fold_change_shrunk" in res else "log2_fold_change"
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    r = res[(res["base_mean"] > 0)].dropna(subset=[col])
    sig = r["padj"].fillna(1) < fdr
    x = np.log10(r["base_mean"].astype(float))
    ax.scatter(x[~sig], r.loc[~sig, col], s=5, color=NS, linewidth=0, rasterized=True, label="not significant")
    ax.scatter(x[sig & (r[col] > 0)], r.loc[sig & (r[col] > 0), col], s=9, color=UP, linewidth=0, label="up")
    ax.scatter(x[sig & (r[col] < 0)], r.loc[sig & (r[col] < 0), col], s=9, color=DOWN, linewidth=0, label="down")
    ax.axhline(0, color=TEXT2, linewidth=0.8, zorder=0)
    ax.set_xlabel("log10 mean of normalized counts")
    ax.set_ylabel("log2 fold change" + (" (shrunk)" if col.endswith("shrunk") else ""))
    ax.set_title(title, fontsize=9, color=TEXT)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path
