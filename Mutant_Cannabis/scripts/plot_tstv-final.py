#!/usr/bin/env python3
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from matplotlib.lines import Line2D
from scipy.stats import percentileofscore

def auto_color(key, cmap_name="tab20"):
    cmap = plt.get_cmap(cmap_name)
    h = int(hashlib.md5(key.encode()).hexdigest(), 16)
    return cmap(h % cmap.N)

def get_color(wgs_rsp):
    wgs_rsp = wgs_rsp.upper()
    if wgs_rsp in COLOR_MAP:
        return COLOR_MAP[wgs_rsp]
    return auto_color(wgs_rsp)

##make it so these are always the colors IF these RSPs appear (ones we are interested in presenting in stable colors)
COLOR_MAP = {
    "RSP11211": "#1f77b4",
    "RSP11212": "#ff7f0e",
    "RSP11213": "#2ca02c",
    "RSP12876": "#e377c2",
    "RSP13487": "#9467bd",
    "RSP13488": "#8c564b",
    "RSP13489": "#d62728",
}

# -----------------------------
# Command-line arguments
# -----------------------------
parser = argparse.ArgumentParser(description="Plot Ts/Tv distribution with highlighted IDs.")
parser.add_argument("--file", type=str, required=True, help="Path to TSV file")
parser.add_argument("--highlight", type=str, default="", help="Comma-separated WGS_RSP IDs to highlight")
parser.add_argument("--privacy", type=str, nargs="+", default=["all"], help="Privacy values to include: private, public, all (default all)")
parser.add_argument("--output", type=str, default="tstv_plot.png", help="Output figure file")
parser.add_argument("--smooth", action="store_true", help="Enable KDE smoothing (default off)")
parser.add_argument("--percentile_file", type=str, default=None,
                    help="Optional TSV file to save percentile ranks of highlighted samples")
args = parser.parse_args()

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv(args.file, sep="\t", comment=None, encoding="utf-8")

# Strip whitespace from column names
df.columns = df.columns.str.strip()

# Strip all string columns
for col in df.select_dtypes(include='object'):
    df[col] = df[col].astype(str).str.strip()

# Force WGS_RSP column to uppercase for matching
df['WGS_RSP_upper'] = df['WGS_RSP'].str.upper()

# -----------------------------
# Parse highlight list
highlight_ids = [x.strip().upper() for x in args.highlight.split(",") if x.strip()]

# -----------------------------
# Separate highlighted rows (always keep)
highlight_df = df[df['WGS_RSP_upper'].isin(highlight_ids)].copy()
highlight_df['highlight'] = True

# Filter non-highlighted rows by privacy
if "all" not in [p.lower() for p in args.privacy]:
    df_filtered = df[~df['WGS_RSP_upper'].isin(highlight_ids) & df['Privacy'].isin(args.privacy)].copy()
else:
    df_filtered = df[~df['WGS_RSP_upper'].isin(highlight_ids)].copy()

# Combine highlighted rows back
plot_df = pd.concat([highlight_df, df_filtered], ignore_index=True)

# Ensure highlight column is boolean (fix ~ error)
plot_df['highlight'] = plot_df['highlight'].fillna(False).astype(bool)

# -----------------------------
# Warn if any requested highlights not found
found_ids = plot_df.loc[plot_df['highlight'], 'WGS_RSP_upper'].tolist()
missing = [x for x in highlight_ids if x not in found_ids]
if missing:
    print(f"Warning: The following IDs were not found in the data: {', '.join(missing)}", file=sys.stderr)

# -----------------------------
# Print sample counts
total_samples = len(plot_df)
highlighted_samples = plot_df['highlight'].sum()
non_highlighted_samples = total_samples - highlighted_samples

print(f"Total samples plotted: {total_samples}", file=sys.stderr)
print(f"Highlighted samples: {highlighted_samples}", file=sys.stderr)
print(f"Non-highlighted samples: {non_highlighted_samples}", file=sys.stderr)

# -----------------------------
# Use RSP_Plus_Strain for legend labels
plot_df['legend_label'] = plot_df['RSP_Plus_Strain']

# -----------------------------
# Determine figure size dynamically
base_width = 12
width_per_highlight = 0.3
fig_width = base_width + highlighted_samples * width_per_highlight
fig_height = 6
plt.figure(figsize=(fig_width, fig_height))

# -----------------------------
# Plotting
sns.set(style="whitegrid", context="talk")

# Color palette
palette = sns.color_palette("tab10", n_colors=max(10, highlighted_samples))

# Plot distribution: KDE if smooth, else histogram
background_vals = plot_df.loc[~plot_df['highlight'], 'Ts/Tv']
if args.smooth:
    sns.kdeplot(background_vals, fill=True, color='lightblue', alpha=0.4, label='Other samples')
    ylabel = "Density"
else:
    plt.hist(background_vals, bins=30, color='lightblue', alpha=0.4, label='Other samples')
    ylabel = "Count"

# -----------------------------
# Plot highlighted points as stacked thick rug lines
highlighted = plot_df[plot_df['highlight']].copy()
highlighted = highlighted.sort_values('Ts/Tv')
highlighted['stack_level'] = highlighted.groupby('Ts/Tv').cumcount()

#rug_height = 0.06
#rug_height = 1
#rug_base = -0.05

if args.smooth:
    # Keep the nice fixed ticks you love
    rug_height = 1
    rug_base = -0.05
else:
    # Scale rug ticks relative to histogram
    counts, bins, patches = plt.hist(plot_df.loc[~plot_df['highlight'], 'Ts/Tv'],
                                     bins=30, color='lightblue', alpha=0.4, label='Other samples')
    ymax = counts.max()
    rug_height = ymax * 0.1       # 10% of axis height per stack level
    rug_base = 0                   # start at bottom of histogram

# Make sure stacked points get a vertical offset
highlighted['stack_level'] = highlighted.groupby('Ts/Tv').cumcount()

for i, (_, row) in enumerate(highlighted.iterrows()):
    color = get_color(row["WGS_RSP"])
    x = row['Ts/Tv']
    level = row['stack_level']
    #color = palette[i % len(palette)]
    y0 = rug_base + level * rug_height
    y1 = y0 + rug_height * 0.85
    plt.vlines(x=x, ymin=y0, ymax=y1, color=color, linewidth=6, alpha=0.9)

# -----------------------------
# Legend
#legend_elements = [Line2D([0], [0], color=palette[i % len(palette)], lw=6,
#                          label=row['legend_label'])
#                   for i, (_, row) in enumerate(highlighted.iterrows())]

from matplotlib.patches import Rectangle

highlighted_sorted = highlighted.sort_values("WGS_RSP")

legend_elements = [
    Rectangle((0,0), width=0.2, height=1, color=get_color(row["WGS_RSP"]),
              label=row['legend_label'])
    for i, (_, row) in enumerate(highlighted_sorted.iterrows())
]

plt.legend(handles=legend_elements, title="Highlighted Samples", bbox_to_anchor=(1.05,1), loc='upper left')

plt.xlabel("Ts/Tv Ratio")
plt.ylabel(ylabel)
plt.title("Distribution of Ts/Tv Ratios")
plt.ylim(bottom=rug_base - 0.03)
plt.tight_layout()
plt.savefig(args.output, dpi=300)
print(f"Plot saved to {args.output}", file=sys.stderr)

# -----------------------------
# Compute percentile ranks and counts if requested
if args.percentile_file:
    ts_tv_all = plot_df['Ts/Tv']
    percentiles = []
    for _, row in highlighted.iterrows():
        value = row['Ts/Tv']
        perc = percentileofscore(ts_tv_all, value, kind='rank')
        num_lower = (ts_tv_all < value).sum()
        num_higher = (ts_tv_all > value).sum()
        percentiles.append({
            'WGS_RSP': row['WGS_RSP'],
            'Ts/Tv': value,
            'Percentile': perc,
            'Num_Lower': num_lower,
            'Num_Higher': num_higher,
            'RSP_Plus_Strain': row['legend_label']
        })
    percentile_df = pd.DataFrame(percentiles)
    percentile_df.to_csv(args.percentile_file, sep="\t", index=False)
    print(f"Percentile ranks and counts saved to {args.percentile_file}", file=sys.stderr)
