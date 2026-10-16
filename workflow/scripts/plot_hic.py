#!/usr/bin/env python3
"""
plot_hic.py — publication-style figures from the pipeline's own text outputs
(juicer_tools dump / pearsons / eigenvector, Arrowhead and HiCCUPS bedpe).
Nothing here recomputes Hi-C statistics; it only draws what juicer_tools wrote.

Figures (each saved as PNG at 600 dpi and as PDF, fonts embedded as TrueType):
    fig1_contact_maps   a) whole-chromosome contact map  b) zoom with TADs/loops
    fig2_compartments   eigenvector track above the Pearson correlation matrix
    fig3_tads_loops     TAD size / count and loop distance / count summaries

Usage:
    python3 plot_hic.py --sample hic1 --chrom chr1 --chrom-sizes hg38.chrom.sizes \
        --matrix chr1_contacts_25000.txt --matrix-res 25000 \
        --pearsons chr1_pearsons_100000.txt --pearsons-res 100000 \
        --eigenvector chr1_eigenvector_250000.txt --eigenvector-res 250000 \
        --tads 25000_blocks.bedpe --loops postprocessed_pixels_5000.bedpe \
        --zoom-start 60000000 --zoom-end 90000000 --outdir results/plots
"""
import argparse
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

# Okabe-Ito colour-blind-safe palette
BLUE, VERMILLION, GREY, BLACK = "#0072B2", "#D55E00", "#7F7F7F", "#000000"
DPI = 600
WIDTH = 7.2  # inches, a two-column journal figure

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 8,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def save(fig, outdir, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def panel_label(ax, letter):
    ax.text(-0.02, 1.04, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="bottom", ha="right")


def chrom_length(path, chrom):
    with open(path) as fh:
        for line in fh:
            name, size = line.split()[:2]
            if name == chrom:
                return int(size)
    raise SystemExit(f"{chrom} not found in {path}")


def load_dense_contacts(path, res, length):
    """juicer_tools dump gives 'start_bp end_bp value' for the upper triangle;
    return a symmetric float32 matrix with zeros and NaNs as NaN."""
    df = pd.read_csv(path, sep=r"\s+", header=None, names=["i", "j", "v"])
    n = int(np.ceil(length / res))
    mat = np.full((n, n), np.nan, dtype=np.float32)
    i = (df["i"].to_numpy() // res).astype(int)
    j = (df["j"].to_numpy() // res).astype(int)
    v = df["v"].to_numpy(dtype=np.float32)
    keep = (i < n) & (j < n)
    mat[i[keep], j[keep]] = v[keep]
    mat[j[keep], i[keep]] = v[keep]
    mat[mat <= 0] = np.nan
    return mat


def read_bedpe(path):
    """First six bedpe columns as numbers; header and comment lines drop out."""
    cols = ["chr1", "x1", "x2", "chr2", "y1", "y2"]
    try:
        df = pd.read_csv(path, sep="\t", header=None, usecols=range(6),
                         names=cols, dtype=str)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=cols)
    for c in ("x1", "x2", "y1", "y2"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna(subset=["x1", "x2", "y1", "y2"]).reset_index(drop=True)


def chrom_key(name):
    m = re.match(r"chr(\d+)$", name)
    if m:
        return (0, int(m.group(1)))
    return (1, name)


def contact_scale(mat):
    logv = np.log10(mat[np.isfinite(mat)])
    return np.percentile(logv, 2), np.percentile(logv, 99.5)


def fig_contact_maps(mat, res, chrom, tads, loops, zoom, sample, outdir):
    fig = plt.figure(figsize=(WIDTH, 3.6))
    gs = GridSpec(1, 2, figure=fig, wspace=0.45)
    vmin, vmax = contact_scale(mat)
    logm = np.log10(mat)
    cmap = plt.get_cmap("YlOrRd").copy()
    cmap.set_bad("white")

    # a) whole chromosome
    ax = fig.add_subplot(gs[0])
    mb = mat.shape[0] * res / 1e6
    im = ax.imshow(logm, cmap=cmap, vmin=vmin, vmax=vmax,
                   extent=[0, mb, mb, 0], interpolation="nearest")
    ax.set_xlabel(f"{chrom} position (Mb)")
    ax.set_ylabel(f"{chrom} position (Mb)")
    ax.set_title(f"{sample}: {chrom}, {res // 1000} kb bins")
    panel_label(ax, "a")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("log$_{10}$ normalized contacts")
    cb.outline.set_linewidth(0.4)

    # b) zoom with TADs (squares on the diagonal) and loops (upper triangle)
    zs, ze = zoom
    s_bin, e_bin = zs // res, min(ze // res, mat.shape[0])
    ax = fig.add_subplot(gs[1])
    im = ax.imshow(logm[s_bin:e_bin, s_bin:e_bin], cmap=cmap, vmin=vmin, vmax=vmax,
                   extent=[zs / 1e6, ze / 1e6, ze / 1e6, zs / 1e6],
                   interpolation="nearest")
    n_tad = n_loop = 0
    if tads is not None and len(tads):
        sel = tads[(tads["chr1"] == chrom) & (tads["x1"] >= zs) & (tads["x2"] <= ze)]
        for _, r in sel.iterrows():
            ax.add_patch(Rectangle((r.x1 / 1e6, r.x1 / 1e6), (r.x2 - r.x1) / 1e6,
                                   (r.x2 - r.x1) / 1e6, fill=False,
                                   edgecolor=BLACK, linewidth=0.5))
        n_tad = len(sel)
    if loops is not None and len(loops):
        sel = loops[(loops["chr1"] == chrom) & (loops["x1"] >= zs) & (loops["y2"] <= ze)]
        # anchor 1 is the upstream anchor: column = anchor 2, row = anchor 1
        ax.scatter((sel.y1 + sel.y2) / 2e6, (sel.x1 + sel.x2) / 2e6, s=18,
                   facecolors="none", edgecolors=BLUE, linewidths=0.7)
        n_loop = len(sel)
    ax.set_xlim(zs / 1e6, ze / 1e6)
    ax.set_ylim(ze / 1e6, zs / 1e6)
    ax.set_xlabel(f"{chrom} position (Mb)")
    ax.set_ylabel(f"{chrom} position (Mb)")
    ax.set_title(f"{zs / 1e6:g}-{ze / 1e6:g} Mb: {n_tad} TADs (black), {n_loop} loops (blue)")
    panel_label(ax, "b")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("log$_{10}$ normalized contacts")
    cb.outline.set_linewidth(0.4)
    save(fig, outdir, "fig1_contact_maps")


def fig_compartments(pearsons, p_res, eig, e_res, chrom, sample, outdir):
    fig = plt.figure(figsize=(WIDTH * 0.62, 5.6))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 6], width_ratios=[40, 1],
                  hspace=0.05, wspace=0.05)
    mb = pearsons.shape[0] * p_res / 1e6

    ax_t = fig.add_subplot(gs[0, 0])
    x = np.arange(len(eig)) * e_res / 1e6
    ev = np.nan_to_num(eig, nan=0.0)
    ax_t.fill_between(x, 0, ev, where=ev > 0, step="post", color=VERMILLION, linewidth=0)
    ax_t.fill_between(x, 0, ev, where=ev < 0, step="post", color=BLUE, linewidth=0)
    ax_t.axhline(0, color=BLACK, linewidth=0.4)
    ax_t.set_xlim(0, mb)
    ax_t.set_ylabel(f"PC1\n({e_res // 1000} kb)")
    ax_t.tick_params(labelbottom=False)
    ax_t.set_title(f"{sample}: {chrom} compartments")
    for side in ("top", "right"):
        ax_t.spines[side].set_visible(False)
    panel_label(ax_t, "a")

    ax = fig.add_subplot(gs[1, 0])
    im = ax.imshow(pearsons, cmap="RdBu_r", vmin=-1, vmax=1,
                   extent=[0, mb, mb, 0], interpolation="nearest")
    ax.set_xlabel(f"{chrom} position (Mb)")
    ax.set_ylabel(f"{chrom} position (Mb)")
    panel_label(ax, "b")
    cax = fig.add_subplot(gs[1, 1])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(f"Pearson correlation ({p_res // 1000} kb bins)")
    cb.outline.set_linewidth(0.4)
    save(fig, outdir, "fig2_compartments")


def fig_tads_loops(tads, loops, sample, outdir):
    fig, axes = plt.subplots(2, 2, figsize=(WIDTH, 4.6))
    fig.subplots_adjust(hspace=0.5, wspace=0.35)

    def counts_bar(ax, df, colour, title):
        c = df["chr1"].value_counts()
        c = c.reindex(sorted(c.index, key=chrom_key))
        ax.bar([n.replace("chr", "") for n in c.index], c.values, color=colour, width=0.8)
        ax.set_xlabel("Chromosome")
        ax.set_ylabel("Count")
        ax.set_title(title)
        ax.tick_params(axis="x", labelsize=5.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    if tads is not None and len(tads):
        size_kb = (tads["x2"] - tads["x1"]) / 1e3
        ax = axes[0, 0]
        ax.hist(size_kb, bins=40, color=GREY, edgecolor="white", linewidth=0.3)
        ax.axvline(size_kb.median(), color=BLACK, linestyle="--", linewidth=0.7)
        ax.set_xlabel("TAD size (kb)")
        ax.set_ylabel("Number of TADs")
        ax.set_title(f"TAD sizes (n = {len(tads)}, median {size_kb.median():.0f} kb)")
        counts_bar(axes[0, 1], tads, GREY, "TADs per chromosome")
    else:
        axes[0, 0].axis("off")
        axes[0, 1].axis("off")
    if loops is not None and len(loops):
        dist_kb = ((loops["y1"] + loops["y2"]) / 2 - (loops["x1"] + loops["x2"]) / 2).abs() / 1e3
        ax = axes[1, 0]
        ax.hist(dist_kb, bins=40, color=BLUE, edgecolor="white", linewidth=0.3)
        ax.axvline(dist_kb.median(), color=BLACK, linestyle="--", linewidth=0.7)
        ax.set_xlabel("Loop anchor distance (kb)")
        ax.set_ylabel("Number of loops")
        ax.set_title(f"Loop spans (n = {len(loops)}, median {dist_kb.median():.0f} kb)")
        counts_bar(axes[1, 1], loops, BLUE, "Loops per chromosome")
    else:
        axes[1, 0].axis("off")
        axes[1, 1].axis("off")
    for ax, letter in zip(axes.flat, "abcd"):
        if ax.axison:
            panel_label(ax, letter)
    fig.suptitle(f"{sample}: Arrowhead TADs and HiCCUPS loops", fontsize=8, y=0.98)
    save(fig, outdir, "fig3_tads_loops")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--sample", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--matrix-res", type=int, required=True)
    ap.add_argument("--pearsons", required=True)
    ap.add_argument("--pearsons-res", type=int, required=True)
    ap.add_argument("--eigenvector", required=True)
    ap.add_argument("--eigenvector-res", type=int, required=True)
    ap.add_argument("--tads", required=True)
    ap.add_argument("--loops")
    ap.add_argument("--zoom-start", type=int, required=True)
    ap.add_argument("--zoom-end", type=int, required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    length = chrom_length(a.chrom_sizes, a.chrom)
    zoom = (a.zoom_start, min(a.zoom_end, length))
    if zoom[0] >= zoom[1]:
        raise SystemExit(f"zoom window {a.zoom_start}-{a.zoom_end} is outside {a.chrom} ({length} bp)")
    tads = read_bedpe(a.tads)
    loops = read_bedpe(a.loops) if a.loops else None

    mat = load_dense_contacts(a.matrix, a.matrix_res, length)
    fig_contact_maps(mat, a.matrix_res, a.chrom, tads, loops, zoom, a.sample, a.outdir)
    del mat

    pearsons = pd.read_csv(a.pearsons, sep=r"\s+", header=None).to_numpy(dtype=np.float32)
    pearsons = pearsons[:, :pearsons.shape[0]]  # drop a trailing-tab empty column
    eig = pd.read_csv(a.eigenvector, header=None)[0].to_numpy(dtype=float)
    fig_compartments(pearsons, a.pearsons_res, eig, a.eigenvector_res, a.chrom, a.sample, a.outdir)

    fig_tads_loops(tads, loops, a.sample, a.outdir)


if __name__ == "__main__":
    main()
