# run command:
# python3 Lee_Ellison_BME163_Assignment_Week4.py -i BME163_Input_Data_Week4.ident -c BME163_Input_Data_Week4.cov -o Lee_Ellison_BME163_Assignment_Week4.png

# import libraries
import matplotlib.pyplot as plt
import numpy as np
import argparse

# create commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument('--identity_file', '-i', type=str, action='store', help='identity file')
parser.add_argument('--coverage_file', '-c', type=str, action='store', help='coverage file')
parser.add_argument('--output_file', '-o', type=str, action='store', help='output file')
args = parser.parse_args()
ident_file = args.identity_file
cov_file = args.coverage_file
out_file = args.output_file

# create figure
plt.style.use('BME163')
figure_width = 6
figure_height = 2.5
plt.figure(figsize=(figure_width, figure_height))

# create panel
panel_width = 4.5
panel_height = 1.5
relative_panel_width = panel_width / figure_width
relative_panel_height = panel_height / figure_height
panel = plt.axes([0.5 / figure_width, 0.15, relative_panel_width, relative_panel_height])
panel.set_xlabel('Subread Coverage')
panel.set_ylabel('Identity (%)')
panel.set_xticks([1, 2, 3, 4])
panel.set_xticklabels(['1-3','4-6','7-9','>=10'])
panel.set_xlim(0.5, 4.5)
panel.set_yticks([75, 80, 85, 90, 95, 100])
panel.set_ylim(75, 100)

# parse coverage file, cov_dict = {'sequence id':int(coverage)}
with open(cov_file, 'r') as f:
    cov_contents = f.read()
cov_contents = cov_contents.split('\n')
cov_contents.pop()
cov_dict = dict()
for i in range(len(cov_contents)):
    line = cov_contents[i]
    line = line.split('\t')
    if line[0] not in cov_dict.keys():
        cov_dict[line[0]] = int(line[1])

# parse identity file, id_dict = {'sequence id':float(percent identity)}
with open(ident_file, 'r') as f:
    id_contents = f.read()
id_contents = id_contents.split('\n')
id_contents.pop()
id_dict = dict()
for i in range(len(id_contents)):
    line = id_contents[i]
    line = line.split('\t')
    if line[0] not in id_dict.keys():
        id_dict[line[0]] = float(line[1])

# group data points into bins (4 swarms for 4 coverage ranges)
bins = {"1-3": [], "4-6": [], "7-9": [], ">=10": []}
for subread in id_dict:
    if subread in cov_dict:
        cov = cov_dict[subread]
        ident = id_dict[subread]
        if 1 <= cov <= 3:
            bins["1-3"].append(ident)
        elif 4 <= cov <= 6:
            bins["4-6"].append(ident)
        elif 7 <= cov <= 9:
            bins["7-9"].append(ident)
        elif cov >= 10:
            bins[">=10"].append(ident)

# Colors for each bin (matched to image)
colors = {
    "1-3":  "#4472C4",
    "4-6":  "#2B7A7A",
    "7-9":  "#ECA020",
    ">=10": "#D03820"
}
bin_positions = {"1-3": 1, "4-6": 2, "7-9": 3, ">=10": 4}

# plot swarms
markersize     = 2
min_dist       = markersize / 72
panel_w_in     = 4.5
panel_h_in     = 1.5
x_range        = 4.0
y_range        = 25.0
x_scale        = panel_w_in / x_range
y_scale        = panel_h_in / y_range
max_half_width = 0.45
def swarm(panel, y_vals, x_center, color):
    y_vals = list(y_vals)
    if len(y_vals) > 500:
        y_vals = list(np.random.choice(y_vals, 500, replace=False))
    sorted_y = sorted(y_vals)
    placed = []
    not_plotted = 0
    for y in sorted_y:
        forbidden = []
        for px, py in placed:
            dy_in = abs(y - py) * y_scale
            if dy_in < min_dist:
                half_gap = np.sqrt(min_dist**2 - dy_in**2) / x_scale
                forbidden.append((px - half_gap, px + half_gap))
        candidates = [x_center]
        for lo, hi in forbidden:
            candidates.append(lo - 1e-9)
            candidates.append(hi + 1e-9)
        candidates = [x for x in candidates
                      if abs(x - x_center) <= max_half_width]
        candidates.sort(key=lambda x: abs(x - x_center))
        placed_flag = False
        for x_try in candidates:
            if all(not (lo <= x_try <= hi) for lo, hi in forbidden):
                placed.append((x_try, y))
                placed_flag = True
                break
        if not placed_flag:
            not_plotted += 1
    if not_plotted > 0:
        print(f"{not_plotted} points could not be plotted at position {x_center}")
    if placed:
        xs = [p[0] for p in placed]
        ys = [p[1] for p in placed]
        panel.plot(xs, ys, 'o',
                   color=color,
                   markersize=markersize,
                   markeredgewidth=0,
                   linestyle='none')
def stratified_subsample(y_vals, n=500, n_strata=50):
    """
    Divide y_vals into n_strata equal-width bins.
    Sample from each stratum proportional to its count,
    so the density shape is preserved across the full range.
    """
    y_arr   = np.array(y_vals)
    y_min, y_max = y_arr.min(), y_arr.max()
    edges   = np.linspace(y_min, y_max, n_strata + 1)
    # assign each point to a stratum
    strata  = []
    for i in range(n_strata):
        lo, hi = edges[i], edges[i + 1]
        # include right edge only on last bin
        if i < n_strata - 1:
            mask = (y_arr >= lo) & (y_arr < hi)
        else:
            mask = (y_arr >= lo) & (y_arr <= hi)
        strata.append(y_arr[mask])
    # compute proportional quota for each stratum
    counts  = np.array([len(s) for s in strata])
    total   = counts.sum()
    quotas  = np.floor(counts / total * n).astype(int)
    # distribute any rounding remainder to the largest strata
    remainder = n - quotas.sum()
    if remainder > 0:
        largest = np.argsort(counts)[::-1]
        for i in range(remainder):
            quotas[largest[i % len(largest)]] += 1
    # draw from each stratum
    sampled = []
    for s, q in zip(strata, quotas):
        if len(s) == 0 or q == 0:
            continue
        q = min(q, len(s))   # can't sample more than available
        sampled.extend(np.random.choice(s, q, replace=False).tolist())
    return sampled
for bin_name, position in bin_positions.items():
    swarm(panel, bins[bin_name], position, colors[bin_name])

# save figure
plt.savefig(out_file, dpi=600)