#!/usr/bin/env python3
"""
plot_hic.py — publication-style figures from the pipeline's own text outputs
(juicer_tools dump / pearsons / eigenvector, Arrowhead and HiCCUPS bedpe).
Nothing here recomputes Hi-C statistics; it only draws what juicer_tools wrote.

Figures (each saved as PNG at 300 dpi and as PDF, fonts embedded as TrueType):
    fig1_contact_maps   a) whole-chromosome map (500 kb bins) with TAD and PC1 tracks
                        b) fine-resolution zoom with a TAD track and loop squares,
                        in the style of the tutorial's Juicebox screenshot
    fig2_compartments   eigenvector track above the Pearson correlation matrix
    fig3_tads_loops     TAD size / count and loop distance / count summaries

Usage:
    python3 plot_hic.py --sample hic1 --chrom chr1 --chrom-sizes hg38.chrom.sizes \
        --matrix chr1_contacts_25000.txt --matrix-res 25000 \
        --zoom-matrix chr1_zoom_92000000_95000000_5000.txt --zoom-res 5000 \
        --pearsons chr1_pearsons_100000.txt --pearsons-res 100000 \
        --eigenvector chr1_eigenvector_250000.txt --eigenvector-res 250000 \
        --tads 25000_blocks.bedpe --loops postprocessed_pixels_5000.bedpe \
        --zoom-start 92000000 --zoom-end 95000000 --outdir results/plots
"""
import argparse
import os
import re
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

# Okabe-Ito colour-blind-safe palette
BLUE, VERMILLION, GREY, BLACK = "#0072B2", "#D55E00", "#7F7F7F", "#000000"
DPI = 300  # also sets the raster resolution of the heatmaps inside the PDFs
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
    ax.text(-0.14, 1.04, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="bottom", ha="right")


def chrom_length(path, chrom):
    with open(path) as fh:
        for line in fh:
            name, size = line.split()[:2]
            if name == chrom:
                return int(size)
    raise SystemExit(f"{chrom} not found in {path}")


def load_dense_contacts(path, res, length, origin=0):
    """juicer_tools dump gives 'start_bp end_bp value' for the upper triangle;
    return a symmetric float32 matrix with zeros and NaNs as NaN. `origin` is
    the bp coordinate of the first bin (non-zero for a dumped region)."""
    df = pd.read_csv(path, sep=r"\s+", header=None, names=["i", "j", "v"])
    n = int(np.ceil(length / res))
    mat = np.full((n, n), np.nan, dtype=np.float32)
    i = ((df["i"].to_numpy() - origin) // res).astype(int)
    j = ((df["j"].to_numpy() - origin) // res).astype(int)
    v = df["v"].to_numpy(dtype=np.float32)
    keep = (i >= 0) & (j >= 0) & (i < n) & (j < n)
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


WHITE_RED = LinearSegmentedColormap.from_list("white_red", ["#FFFFFF", "#D7191C"])
WHITE_RED.set_bad("white")
GREEN = "#009E73"


def coarsen(mat, factor):
    """Mean of measured (non-NaN) bins over factor x factor blocks."""
    n = mat.shape[0] // factor * factor
    blocks = mat[:n, :n].reshape(n // factor, factor, n // factor, factor)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(blocks, axis=(1, 3))


def cap(values, pct):
    return float(np.percentile(values[np.isfinite(values)], pct))


def domain_track(ax, tads, chrom, start, end):
    """Green blocks for TADs, like the tutorial's Juicebox annotation track."""
    sel = tads[(tads["chr1"] == chrom) & (tads["x2"] >= start) & (tads["x1"] <= end)]
    for _, r in sel.iterrows():
        ax.add_patch(Rectangle((max(r.x1, start) / 1e6, 0), (min(r.x2, end) - max(r.x1, start)) / 1e6,
                               1, facecolor=GREEN, edgecolor="white", linewidth=0.3))
    ax.set_xlim(start / 1e6, end / 1e6)
    ax.set_ylim(0, 1)
    ax.axis("off")
    return len(sel)


def eigen_track(ax, eig, e_res, xmax_mb):
    x = np.arange(len(eig)) * e_res / 1e6
    ev = np.nan_to_num(eig, nan=0.0)
    ax.fill_between(x, 0, ev, where=ev > 0, step="post", color=VERMILLION, linewidth=0)
    ax.fill_between(x, 0, ev, where=ev < 0, step="post", color=BLUE, linewidth=0)
    ax.set_xlim(0, xmax_mb)
    ax.axis("off")


def fig_contact_maps(mat, res, zmat, zres, eig, e_res, chrom, tads, loops, zoom, sample, outdir):
    """a) whole chromosome, binned coarsely, with TAD and eigenvector tracks;
    b) the zoom window at fine resolution with a TAD track and loop squares."""
    fig = plt.figure(figsize=(WIDTH, 3.9))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[0.35, 0.35, 10], hspace=0.04,
                  wspace=0.22)

    # a) whole chromosome at ~500 kb, linear colour scale capped at the 95th percentile
    factor = max(1, 500000 // res)
    coarse = coarsen(mat, factor)
    mb = coarse.shape[0] * factor * res / 1e6
    ax_t = fig.add_subplot(gs[0, 0])
    domain_track(ax_t, tads, chrom, 0, int(mb * 1e6))
    ax_e = fig.add_subplot(gs[1, 0])
    eigen_track(ax_e, eig, e_res, mb)
    ax = fig.add_subplot(gs[2, 0])
    im = ax.imshow(coarse, cmap=WHITE_RED, vmin=0, vmax=cap(coarse, 95),
                   extent=[0, mb, mb, 0], interpolation="nearest")
    ax.set_xlabel(f"{chrom} position (Mb)")
    ax.set_ylabel(f"{chrom} position (Mb)")
    ax_t.set_title(f"{sample}: {chrom}, {factor * res // 1000} kb bins "
                   "(TADs green; PC1 red/blue)", fontsize=7)
    panel_label(ax_t, "a")
    cax = ax.inset_axes([1.03, 0, 0.035, 1])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("normalized contacts (mean per bin)")
    cb.outline.set_linewidth(0.4)

    # b) zoom at fine resolution
    zs, ze = zoom
    ax_t = fig.add_subplot(gs[0, 1])
    n_tad = domain_track(ax_t, tads, chrom, zs, ze)
    fig.add_subplot(gs[1, 1]).axis("off")
    ax = fig.add_subplot(gs[2, 1])
    im = ax.imshow(zmat, cmap=WHITE_RED, vmin=0, vmax=cap(zmat, 98),
                   extent=[zs / 1e6, ze / 1e6, ze / 1e6, zs / 1e6],
                   interpolation="nearest")
    n_loop = 0
    if loops is not None and len(loops):
        sel = loops[(loops["chr1"] == chrom) & (loops["x1"] >= zs) & (loops["y2"] <= ze)]
        # a loop is drawn on both sides of the diagonal, as Juicebox does
        for a, b in ((sel.y1 + sel.y2, sel.x1 + sel.x2), (sel.x1 + sel.x2, sel.y1 + sel.y2)):
            ax.scatter(a / 2e6, b / 2e6, s=30, marker="s", facecolors="none",
                       edgecolors="#00B7C7", linewidths=0.8)
        n_loop = len(sel)
    ax.set_xlim(zs / 1e6, ze / 1e6)
    ax.set_ylim(ze / 1e6, zs / 1e6)
    ax.set_xlabel(f"{chrom} position (Mb)")
    ax.set_ylabel(f"{chrom} position (Mb)")
    ax_t.set_title(f"{zs / 1e6:g}-{ze / 1e6:g} Mb, {zres // 1000} kb bins: "
                   f"{n_tad} TADs, {n_loop} loops (squares)", fontsize=7)
    panel_label(ax_t, "b")
    cax = ax.inset_axes([1.03, 0, 0.035, 1])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("normalized contacts")
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
    # Off-diagonal correlations are small; scale to their 98th percentile
    # (symmetric about 0) so the plaid pattern is visible. The diagonal clips.
    lim = float(np.nanpercentile(np.abs(pearsons), 98))
    im = ax.imshow(pearsons, cmap="RdBu_r", vmin=-lim, vmax=lim,
                   extent=[0, mb, mb, 0], interpolation="nearest")
    ax.set_xlabel(f"{chrom} position (Mb)")
    ax.set_ylabel(f"{chrom} position (Mb)")
    panel_label(ax, "b")
    cax = fig.add_subplot(gs[1, 1])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(f"Pearson correlation ({p_res // 1000} kb bins; clipped at ±{lim:.2f})")
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
        ax.tick_params(axis="x", labelsize=5.5, labelrotation=90)
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
    ap.add_argument("--zoom-matrix", required=True,
                    help="juicer_tools dump of the zoom window only")
    ap.add_argument("--zoom-res", type=int, required=True)
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

    eig = pd.read_csv(a.eigenvector, header=None)[0].to_numpy(dtype=float)
    mat = load_dense_contacts(a.matrix, a.matrix_res, length)
    zmat = load_dense_contacts(a.zoom_matrix, a.zoom_res, zoom[1] - zoom[0], origin=zoom[0])
    fig_contact_maps(mat, a.matrix_res, zmat, a.zoom_res, eig, a.eigenvector_res,
                     a.chrom, tads, loops, zoom, a.sample, a.outdir)
    del mat

    pearsons = pd.read_csv(a.pearsons, sep=r"\s+", header=None).to_numpy(dtype=np.float32)
    pearsons = pearsons[:, :pearsons.shape[0]]  # drop a trailing-tab empty column
    fig_compartments(pearsons, a.pearsons_res, eig, a.eigenvector_res, a.chrom, a.sample, a.outdir)

    fig_tads_loops(tads, loops, a.sample, a.outdir)


if __name__ == "__main__":
    main()
