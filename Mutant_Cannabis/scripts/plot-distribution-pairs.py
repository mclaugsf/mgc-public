#!/usr/bin/env python3

import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys
import hashlib
from scipy.stats import percentileofscore
from matplotlib.patches import Rectangle

# -----------------------------
# Color handling
# -----------------------------
COLOR_MAP = {
    "RSP11211": "#1f77b4",
    "RSP11212": "#ff7f0e",
    "RSP11213": "#2ca02c",
    "RSP12876": "#e377c2",
    "RSP13487": "#9467bd",
    "RSP13488": "#8c564b",
    "RSP13489": "#d62728",
}

def auto_color(key, cmap_name="tab20"):
    cmap = plt.get_cmap(cmap_name)
    h = int(hashlib.md5(key.encode()).hexdigest(), 16)
    return cmap(h % cmap.N)

def get_color(sample):
    s = sample.upper()
    if s in COLOR_MAP:
        return COLOR_MAP[s]
    return auto_color(s)

# -----------------------------
# CLI
# -----------------------------
parser = argparse.ArgumentParser(
    description="Plot genetic distance distribution for one query vs multiple hit samples"
)
parser.add_argument("--file", required=True, help="Pairwise distance CSV from R")
parser.add_argument("--query-sample", required=True, help="Focal sample (SAMPLE1)")
parser.add_argument(
    "--hit-sample",
    action="append",
    required=True,
    help="Hit sample(s) (can be repeated)"
)
parser.add_argument("--output", default="distance_distribution.png")
parser.add_argument("--smooth", action="store_true")
parser.add_argument("--percentile_file", default=None)
args = parser.parse_args()

query = args.query_sample
hits = sorted(set(args.hit_sample))

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv(args.file)
required = {"Sample1", "Sample2", "Distance", "EmpiricalPvalue"}
missing = required - set(df.columns)
if missing:
    sys.exit(f"Missing required columns: {', '.join(missing)}")

# Remove self distances
df = df[df["Sample1"] != df["Sample2"]].copy()

# Canonicalize pairs (avoid duplicates)
df["pair_key"] = df.apply(
    lambda r: "_".join(sorted([r["Sample1"], r["Sample2"]])),
    axis=1
)
df = df.drop_duplicates("pair_key")

# -----------------------------
# Identify highlighted pairs
# -----------------------------
def is_hit(row):
    return (
        (row["Sample1"] == query and row["Sample2"] in hits) or
        (row["Sample2"] == query and row["Sample1"] == hits)
    )

df["highlight"] = df.apply(is_hit, axis=1)
highlighted = df[df["highlight"]].copy()
background = df[~df["highlight"]].copy()

if highlighted.empty:
    sys.exit("No query–hit pairs found in the data")

# Assign HitSample for coloring
highlighted["HitSample"] = highlighted.apply(
    lambda r: r["Sample2"] if r["Sample1"] == query else r["Sample1"],
    axis=1
)

# -----------------------------
# Compute stack levels for highlighted points
# -----------------------------
highlighted = highlighted.sort_values("Distance").copy()
tolerance = 1e-5  # distances within this difference are considered equal

stack_levels = []
bins = []

for d in highlighted["Distance"]:
    for i, b in enumerate(bins):
        if abs(d - b) <= tolerance:
            stack_levels.append(i)
            break
    else:
        bins.append(d)
        stack_levels.append(len(bins) - 1)

highlighted["stack_level"] = stack_levels

# -----------------------------
# Figure setup
# -----------------------------
base_width = 12
fig_width = base_width + len(hits) * 0.4
fig_height = 6
plt.figure(figsize=(fig_width, fig_height))
sns.set(style="whitegrid", context="talk")

# -----------------------------
# Background distribution
# -----------------------------
vals = background["Distance"]

if args.smooth:
    sns.kdeplot(vals, fill=True, color="lightblue", alpha=0.4)
    ylabel = "Density"
else:
    plt.hist(vals, bins=40, color="lightblue", alpha=0.4)
    ylabel = "Count"

ax = plt.gca()
ymin, ymax = ax.get_ylim()

# -----------------------------
# Colored stacked rugs (tight stacking)
# -----------------------------
#line_fraction = 0.35  # rug line height as fraction of stack height
rug_base = ymin + (ymax - ymin) * 0.01  # small margin above bottom

line_fraction = 0.5  # taller lines
rug_height = (ymax - ymin) * 0.08  # increase vertical space per stack level


#rug_height = (ymax - ymin) * 0.05       # vertical space per stack level

for _, row in highlighted.iterrows():
    color = get_color(row["HitSample"])
    x = row["Distance"]
    level = row["stack_level"]

    line_height = rug_height * line_fraction
    y0 = rug_base + level * line_height
    y1 = y0 + line_height

    plt.vlines(
        x=x,
        ymin=y0,
        ymax=y1,
        color=color,
        linewidth=6,
        alpha=0.95
    )

# -----------------------------
# Legend (hit samples only)
# -----------------------------
legend_elements = [
    Rectangle((0, 0), 0.3, 1, color=get_color(h), label=h)
    for h in hits
]

plt.legend(handles=legend_elements, title="Hit samples",
           bbox_to_anchor=(1.02, 1), loc="upper left")

# -----------------------------
# Labels
# -----------------------------
plt.xlabel("Genetic distance (1 − IBS)")
plt.ylabel(ylabel)
plt.title(f"Genetic distance distribution\nQuery sample: {query}")
plt.tight_layout()
plt.savefig(args.output, dpi=300)
plt.close()
print(f"Plot saved to {args.output}", file=sys.stderr)

# -----------------------------
# Percentiles (optional)
# -----------------------------
if args.percentile_file:
    rows = []
    all_distances = df["Distance"]

    for _, row in highlighted.iterrows():
        rows.append({
            "Query": query,
            "HitSample": row["HitSample"],
            "Distance": row["Distance"],
            "Percentile": percentileofscore(all_distances, row["Distance"], kind="rank"),
            "EmpiricalPvalue": row["EmpiricalPvalue"]
        })

    out = pd.DataFrame(rows).sort_values("Distance")
    out.to_csv(args.percentile_file, sep="\t", index=False)
    print(f"Percentiles written to {args.percentile_file}", file=sys.stderr)
