#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.cm import viridis
from scipy.interpolate import interp1d
import csv
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute distances for a query strain against a distance TSV"
    )

    parser.add_argument(
        "--query_rsp",
        help="Query RSP identifier (e.g., RSP000123)"
    )

    parser.add_argument(
        "--query_strain",
        help="Query strain name or ID"
    )

    parser.add_argument(
        "--distance_tsv_input",
        help="Path to input distance TSV file"
    )

    return parser.parse_args()

if __name__ == "__main__":

    args = parse_args()

    print("Query RSP:", args.query_rsp)
    print("Query strain:", args.query_strain)
    print("Distance TSV:", args.distance_tsv_input)

    nn20 = []
    with open(args.distance_tsv_input, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            nn20.append({
                "name": row["name"],
                "distance": float(row["distance"]),
                "id": row["id"],
                "pval": float(row["pval"])
            })

    # -----------------------------
    # Prepare arrays
    # -----------------------------
    distances = [d["distance"] for d in nn20]
    pvals = [d["pval"] for d in nn20]
    labels = [f"{d['name']} ({d['id']})" for d in nn20]

    norm = Normalize(vmin=min(pvals), vmax=max(pvals))

    # -----------------------------
    # Compute distance-aware bins
    # -----------------------------
    NUM_BINS = 5
    dist_array = np.array(distances)
    pval_array = np.array(pvals)

    bin_edges = np.quantile(dist_array, np.linspace(0, 1, NUM_BINS + 1))
    bins = []
    for i in range(NUM_BINS):
        lo = bin_edges[i]
        hi = bin_edges[i + 1]
        mask = (dist_array >= lo) & (dist_array < hi)
        if mask.any():
            bin_p = np.median(pval_array[mask])
        else:
            bin_p = np.nan
        bins.append((lo, hi, bin_p))

    # Interpolation for smooth strip
    bin_centers = np.array([(lo + hi) / 2 for lo, hi, _ in bins])
    bin_pvals = np.array([p for _, _, p in bins])
    interp_func = interp1d(bin_centers, bin_pvals, kind='linear', fill_value="extrapolate")

    # -----------------------------
    # Figure + margins
    # -----------------------------
    # dynamically set figure height
    POINT_HEIGHT = 0.35  # inches per point
    fig_height = max(6, len(labels) * POINT_HEIGHT)
    fig, ax = plt.subplots(figsize=(10, fig_height))

    #fig, ax = plt.subplots(figsize=(10, 6))
    fig.subplots_adjust(left=0.35, right=0.92, top=0.90, bottom=0.38)

    # -----------------------------
    # Plot points
    # -----------------------------
    for y, (x, p) in enumerate(zip(distances, pvals)):
        ax.scatter(
            x, y,
            s=120,
            color=viridis(norm(p)),
            edgecolor="black",
            zorder=3
        )
        OFFSET = (max(distances) - min(distances)) * 0.025
        ax.text(
            x + OFFSET,
            y,
            f"p={p:.2f}",
            va="center",
            ha="left",
            fontsize=8
        )

    # -----------------------------
    # Axes formatting
    # -----------------------------
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    #BUFFER = 0.010
    x_range = max(distances) - min(distances)
    BUFFER = x_range * 0.18   # try 5–10%
    ax.set_xlim(min(distances) - BUFFER, max(distances) + BUFFER)
    ###

    ax.set_xlim(min(distances) - BUFFER, max(distances) + BUFFER)
    ax.set_xticks(np.round(np.linspace(min(distances), max(distances), 5), 3))
    ax.set_xlabel("Genetic distance (1 - IBS)")
    ax.set_title("Nearest genetic distances to %s (%s)"%(args.query_rsp, args.query_strain))
    ax.grid(axis="x", linestyle="--", alpha=0.5)
    ax.set_axisbelow(True)

    # -----------------------------
    # Smooth blended p-value strip
    # -----------------------------
    x_vals = np.linspace(min(distances), max(distances), 500)
    ymin, ymax = -0.28, -0.18

    for i in range(len(x_vals)-1):
        seg_start = x_vals[i]
        seg_end = x_vals[i+1]
        seg_p = float(interp_func((seg_start + seg_end)/2))
        ax.axvspan(
            seg_start, seg_end,
            ymin=ymin, ymax=ymax,
            transform=ax.get_xaxis_transform(),
            color=viridis(norm(seg_p)),
            clip_on=False
        )

    # Strip label
    ax.text(
        min(distances) + 0.001,
        -0.33,
        "Empirical p-value",
        transform=ax.get_xaxis_transform(),
        ha="left",
        va="top",
        fontsize=10,
        fontweight="bold"
    )

    # Min / Max labels at edges
    label_y = -0.16
    ax.text(
        min(distances),
        label_y,
        f"min: {min(pvals):.2f}",
        transform=ax.get_xaxis_transform(),
        ha="left",
        va="bottom",
        fontsize=8,
        color="black"
    )
    ax.text(
        max(distances),
        label_y,
        f"max: {max(pvals):.2f}",
        transform=ax.get_xaxis_transform(),
        ha="right",
        va="bottom",
        fontsize=8,
        color="black"
    )

    # -----------------------------
    # Output
    # -----------------------------
    plt.savefig("%s-genetic_distance_with_pvals_blended.png"%(args.query_rsp), dpi=300)
    plt.close()
