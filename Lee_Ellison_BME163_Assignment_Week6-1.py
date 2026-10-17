# run command used during testing:
# * python3 Lee_Ellison_BME163_Assignment_Week6.py -o Lee_Ellison_BME163_Assignment_Week6.png -p BME163_Input_Data_Assignment6.phase -e BME163_Input_Data_Assignment6.exp -g Arntl,Clock,Dbp,Per1,Per2,Per3,Nr1d1,Insig2,Cry1

import matplotlib.pyplot as plt
import numpy as np
import argparse

# ─── Colormap (viridis) ──────────────────────────────────────────────────────
# Builds a list of 101 (R,G,B) tuples spanning the viridis color gradient.
# Index 0 = darkest purple (Min), index 100 = bright yellow (Max).
# colormaps.py is assumed to be in the same directory at runtime;
# the function is reproduced here so the script is self-contained.
def viridis():
    viridis5 = (253/255, 231/255, 37/255)
    viridis4 = (94/255, 201/255, 98/255)
    viridis3 = (33/255, 145/255, 140/255)
    viridis2 = (59/255, 82/255, 139/255)
    viridis1 = (68/255,   1/255, 84/255)

    R1 = np.linspace(viridis1[0], viridis2[0], 26)
    G1 = np.linspace(viridis1[1], viridis2[1], 26)
    B1 = np.linspace(viridis1[2], viridis2[2], 26)

    R2 = np.linspace(viridis2[0], viridis3[0], 26)
    G2 = np.linspace(viridis2[1], viridis3[1], 26)
    B2 = np.linspace(viridis2[2], viridis3[2], 26)

    R3 = np.linspace(viridis3[0], viridis4[0], 26)
    G3 = np.linspace(viridis3[1], viridis4[1], 26)
    B3 = np.linspace(viridis3[2], viridis4[2], 26)

    R4 = np.linspace(viridis4[0], viridis5[0], 26)
    G4 = np.linspace(viridis4[1], viridis5[1], 26)
    B4 = np.linspace(viridis4[2], viridis5[2], 26)

    R = np.concatenate((R1[:-1], R2[:-1], R3[:-1], R4), axis=None)
    G = np.concatenate((G1[:-1], G2[:-1], G3[:-1], G4), axis=None)
    B = np.concatenate((B1[:-1], B2[:-1], B3[:-1], B4), axis=None)

    viridis_cmap = []
    for index in range(0, 101, 1):
        viridis_cmap.append((R[index], G[index], B[index]))
    return viridis_cmap

# ─── Argument Parsing ────────────────────────────────────────────────────────
# -p : path to .phase file (Ensembl_ID → Peak_phase CT)
# -e : path to .exp file   (gene expression FPKM values)
# -g : comma-separated gene symbols to mark with y-axis tick labels
# -o : output PNG filename
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--phase',  default='BME163_Input_Data_Assignment6.phase',
                    help='Path to .phase file')
parser.add_argument('-e', '--exp',    default='BME163_Input_Data_Assignment6.exp',
                    help='Path to .exp file')
parser.add_argument('-g', '--genes',  default='Arntl,Clock,Dbp,Per1,Per2,Per3,Nr1d1,Insig2,Cry1',
                    help='Comma-separated list of genes for tick labels')
parser.add_argument('-o', '--output', default='Lee_Ellison_BME163_Assignment_Week6.png',
                    help='Output PNG filename')
args = parser.parse_args()

# Build list of gene symbols to label; empty list if -g not provided
label_genes = [g.strip() for g in args.genes.split(',') if g.strip()] if args.genes else []

# ─── Load Colormap ───────────────────────────────────────────────────────────
viridis_cmap = viridis()   # list of 101 (R,G,B) tuples, indexed 0-100

# ─── Read Phase File ─────────────────────────────────────────────────────────
# Builds a dict: Ensembl_ID -> Peak_phase(CT) as float.
# Used to (a) filter genes to those with a known phase, and (b) sort them.
phase_dict = {}
with open(args.phase, 'r') as f:
    f.readline()  # skip header line
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            phase_dict[parts[0]] = float(parts[1])

# ─── Read Expression File ────────────────────────────────────────────────────
# Collects gene symbol, Ensembl ID, and FPKM values (CT0..CT21) for each gene.
# FPKM columns are identified dynamically by the FPKM_CT* header prefix.
gene_symbols = []   # gene symbol strings
ensembl_ids  = []   # corresponding Ensembl IDs
fpkm_matrix  = []   # list of 8-value FPKM lists per gene

with open(args.exp, 'r') as f:
    header    = f.readline().strip().split('\t')
    fpkm_cols = [i for i, h in enumerate(header) if h.startswith('FPKM_CT')]

    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < max(fpkm_cols) + 1:
            continue
        gene_symbols.append(parts[0])
        ensembl_ids.append(parts[1])
        fpkm_matrix.append([float(parts[i]) for i in fpkm_cols])

gene_symbols = np.array(gene_symbols)
ensembl_ids  = np.array(ensembl_ids)
fpkm_matrix  = np.array(fpkm_matrix)   # shape: (N_genes, 8)

# ─── Per-Gene Normalization ──────────────────────────────────────────────────
# Each gene's 8 FPKM values are independently scaled to [0, 100].
# Formula: normalized = ((value - row_min) / (row_max - row_min)) * 100
# Maps the gene's minimum expression to 0 and maximum to 100,
# used as an integer index into the 101-entry viridis colormap.
normalized = np.zeros_like(fpkm_matrix)
for i in range(fpkm_matrix.shape[0]):
    row = fpkm_matrix[i]
    rmin, rmax = row.min(), row.max()
    if rmax - rmin == 0:
        normalized[i] = 0.0
    else:
        normalized[i] = ((row - rmin) / (rmax - rmin)) * 100.0

# ─── Filter and Sort by Peak Phase ──────────────────────────────────────────
# Only genes present in the .phase file are kept.
# Sorted ascending by Peak_phase(CT): lowest phase -> index 0 -> bottom of
# the linear heatmap and clockwise start of the circular heatmap.
# For the example -g list this yields (bottom to top):
#   Arntl, Cry1, Per2, Per1, Per3, Dbp, Nr1d1, Insig2, Clock
phases = np.array([phase_dict.get(eid, np.nan) for eid in ensembl_ids])

has_phase    = ~np.isnan(phases)
gene_symbols = gene_symbols[has_phase]
ensembl_ids  = ensembl_ids[has_phase]
normalized   = normalized[has_phase]
phases       = phases[has_phase]

sort_order   = np.argsort(phases)           # ascending phase order
sort_order = sort_order[::-1]
gene_symbols = gene_symbols[sort_order]
normalized   = normalized[sort_order]       # shape: (N_sorted, 8)
phases       = phases[sort_order]

N      = len(gene_symbols)     # total genes after filtering
n_time = normalized.shape[1]   # number of time points (8)

# ─── Apply BME163 Style Sheet ────────────────────────────────────────────────
plt.style.use('BME163')

# ─── Figure and Panel Layout ─────────────────────────────────────────────────
# Exact geometry as specified in the assignment (figure-fraction units).
fig_1  = plt.figure(figsize=(5, 3))
panel1 = plt.axes([0.7/5,  0.3/3,  0.75/5, 2.5/3], frameon=True)   # linear heatmap
panel2 = plt.axes([1.5/5,  1.45/3, 0.1/5,  0.2/3], frameon=True)   # colorbar legend
# panel3 = plt.axes([2.2/5,  0.3/3,  2.5/5,  2.5/3], frameon=False)  # circular heatmap

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 1 — Linear Heatmap
# ═══════════════════════════════════════════════════════════════════════════════
# Layout:
#   x-axis -> time (8 columns, each occupying width 1, from x=0 to x=8)
#   y-axis -> genes sorted by phase (row gene_idx spans y=[gene_idx, gene_idx+1])
#   color  -> viridis colormap indexed by normalized FPKM value (0-100)
#
# Each cell is drawn as a horizontal .plot() line at the row midpoint (y+0.5),
# spanning the full column width [t, t+1].

panel1.set_xlim(0, n_time)
panel1.set_ylim(0, N)

for gene_idx in range(N):
    for t in range(n_time):
        color = viridis_cmap[int(round(normalized[gene_idx, t]))]
        # Horizontal segment spanning the full width of this time-point's column,
        # centered vertically in this gene's row.
        panel1.plot([t, t + 1], [gene_idx + 0.5, gene_idx + 0.5],
                    color=color, linewidth=0.5, solid_capstyle='butt')

# ── Panel 1 X-axis ───────────────────────────────────────────────────────────
# 8 time-points (CT0, CT3, CT6, ..., CT21) each occupy one x-unit.
# Major labeled ticks are centered on CT0, CT6, CT12, CT18 columns -> x = 0.5, 2.5, 4.5, 6.5.
# Minor unlabeled ticks (same size) are centered on the intervening columns
# (CT3, CT9, CT15, CT21) -> x = 1.5, 3.5, 5.5, 7.5.
# Visual pattern along the axis: 0, tick, 6, tick, 12, tick, 18, tick.

major_x = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]

panel1.set_xticks(major_x)
panel1.set_xticklabels(['0', '', '6', '', '12', '', '18', ''])
panel1.set_xlabel('CT')

# ── Panel 1 Y-axis ───────────────────────────────────────────────────────────
# Tick marks and labels only for genes in the -g list.
# Each gene's tick is at its row midpoint (sorted_index + 0.5).
# Ascending phase sort: lowest phase -> bottom, highest -> top.
# Example order bottom to top: Arntl, Cry1, Per2, Per1, Per3, Dbp, Nr1d1, Insig2, Clock.

ytick_positions = []
ytick_labels    = []
for gname in label_genes:
    idxs = np.where(gene_symbols == gname)[0]
    if len(idxs) > 0:
        ytick_positions.append(idxs[0] + 0.5)
        ytick_labels.append(gname)

panel1.set_yticks(ytick_positions)
panel1.set_yticklabels(ytick_labels)
panel1.yaxis.set_tick_params(length=2)

# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 2 — Colorbar Legend
# ═══════════════════════════════════════════════════════════════════════════════
# A small vertical gradient bar drawn with stacked horizontal .plot() lines.
# ylim is [0, 100]; each step i draws the color for viridis index i.
# Tick marks are shown at y=0 (Min) and y=100 (Max) with text labels to the right.

n_legend_steps = 100
panel2.set_xlim(0, 1)
panel2.set_ylim(0, n_legend_steps)
panel2.set_xticks([])

# Draw the gradient: one colored line per viridis step
for i in range(n_legend_steps):
    panel2.plot([0, 1], [i + 0.5, i + 0.5],
                color=viridis_cmap[i], linewidth=0.5, solid_capstyle='butt')

# Tick marks at bottom (Min) and top (Max) only; labels suppressed in favor of text
panel2.set_yticks([0, n_legend_steps])
panel2.set_yticklabels(['Min', 'Max'])
panel2.tick_params(axis='y', which='both', 
                   left=False,      # Hide left ticks
                   labelleft=False, # Hide left labels
                   right=True,      # Show right ticks
                   labelright=True, # Show right labels
                   length=2)
"""
# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 3 — Circular Heatmap
# ═══════════════════════════════════════════════════════════════════════════════
# Same gene order and normalized values as Panel 1, mapped onto a circle:
#   Angular position -> gene (ascending phase, stepping clockwise from top)
#   Radial band      -> time-point (innermost = CT0, outermost = CT21)
#   Color            -> viridis-indexed normalized FPKM (same as Panel 1)
#
# Geometry:
#   inner_radius = 24  (inner edge of CT0 band)
#   band_width   = 1   (each time-point adds 1 radial unit)
#   outer_radius = 32  (24 + 8 time-points)
#
# Each arc segment is drawn as a .plot() line at the radial midpoint of its band.

panel3.set_xlim(-36.1, 36.1)
panel3.set_ylim(-36.1, 36.1)
panel3.set_aspect('equal')
panel3.axis('off')

inner_radius = 24
band_width   = 1

# Clockwise angular step per gene, starting from the top (pi/2)
dtheta = 2 * np.pi / N

for gene_idx in range(N):
    theta_start = np.pi / 2 - gene_idx * dtheta
    theta_end   = np.pi / 2 - (gene_idx + 1) * dtheta

    for t in range(n_time):
        color = viridis_cmap[int(round(normalized[gene_idx, t]))]
        r_mid = inner_radius + t * band_width + band_width / 2.0
        thetas = np.linspace(theta_start, theta_end, 3)
        panel3.plot(r_mid * np.cos(thetas), r_mid * np.sin(thetas),
                    color=color, linewidth=3.6, solid_capstyle='butt')

# ── Circular Panel Gene Labels ───────────────────────────────────────────────
# Text annotations just outside the outermost band for each gene in -g.
# Rotated tangentially so text reads along the arc.
label_radius = inner_radius + n_time * band_width + 1.5

for gname in label_genes:
    idxs = np.where(gene_symbols == gname)[0]
    if len(idxs) == 0:
        continue
    gene_idx  = idxs[0]
    theta_mid = np.pi / 2 - (gene_idx + 0.5) * dtheta

    panel3.text(label_radius * np.cos(theta_mid),
                label_radius * np.sin(theta_mid),
                gname,
                ha='center', va='center',
                fontsize=5,
                rotation=np.degrees(theta_mid) - 90,
                rotation_mode='anchor')
"""
# ─── Save Figure ─────────────────────────────────────────────────────────────
# No bbox_inches argument: figsize=(5,3) at dpi=600 gives exactly 3000x1800 px,
# matching the answer6.png dimensions.
plt.savefig(args.output, dpi=600)
print(f"Saved figure to {args.output}")