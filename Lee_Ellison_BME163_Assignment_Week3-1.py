# run command:
#     python3 Lee_Ellison_BME163_Assignment_Week3.py -p BME163_Input_Data_Week3.position.tsv -c BME163_Input_Data_Week3.celltype.tsv -o Lee_Ellison_BME163_Assignment_Week3.png

# import libraries
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.patheffects as patheffects
import numpy as np
import argparse

# create commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument('--position_file', '-p', type=str, action='store', help='position file')
parser.add_argument('--cell_type_file', '-c', type=str, action='store', help='cell type file')
parser.add_argument('--out_file', '-o', type=str, action='store', help='output file')
args = parser.parse_args()
pos_file = args.position_file
ct_file = args.cell_type_file
out_file = args.out_file

# create figure: 3in x 3in, BME163 style sheet
plt.style.use('BME163')
figure_width = 5
figure_height = 3
plt.figure(figsize=(figure_width, figure_height))

# create panels
panel_width = 1.5
panel_height = 1.5
# left panel is cell type cluster plot
left_panel = plt.axes([0.5 / figure_width, 0.5 / figure_height, panel_width / figure_width, panel_height / figure_height]) # i will have to likely figure out the max and min in my data and plot the axes ticks at some small distance above these maxima
left_panel.set_xlabel('tSNE 2')
left_panel.set_ylabel('tSNE 1')
left_panel.set_xticks([-20, 0, 20])
left_panel.set_yticks([-40, -20, 0, 20])
left_panel.set_xlim(-30, 30)
left_panel.set_ylim(-40, 30)
# middle panel is cell position cluster plot
middle_panel = plt.axes([2.5 / figure_width, 0.5 / figure_height, panel_width / figure_width, panel_height / figure_height]) # i will have to likely figure out the max and min in my data and plot the axes ticks at some small distance above these maxima
middle_panel.set_xlabel('tSNE 2')
middle_panel.set_ylabel('tSNE 1')
middle_panel.set_xticks([-20, 0, 20])
middle_panel.set_yticks([-40, -20, 0, 20])
middle_panel.set_xlim(-30, 30)
middle_panel.set_ylim(-40, 30)
# right panel is heatmap legend
right_panel = plt.axes([4.0 / figure_width, 1.1 / figure_height, 0.1 / figure_width, 0.3 / figure_height])
right_panel.tick_params(right=True, labelright=True, left=False, labelleft=False, bottom=False, labelbottom=False)
y_min, y_max = right_panel.get_ylim()
right_panel.set_yticks([y_min, y_max], labels=['Min', 'Max'])

# parse cell type file: make a dictionary named cells_dict where keys are barcodes (all seem to be unique in file) and values are list, first value in each list is cell type as string
# ex: {'CGATCG':['monocyte'], 'CGATATC':['tCell']}
with open(ct_file, 'r') as f:
    ct_contents = f.read()
ct_contents = ct_contents.split('\n')
ct_contents.pop(0)
ct_contents.pop()
for line in range(len(ct_contents)):
    ct_contents[line] = ct_contents[line].split('\t')
cells_dict = dict()
for cell in ct_contents:
    this_ct_cell_barcode = cell[2]
    this_ct_cell_type = cell[1]
    cells_dict[this_ct_cell_barcode] = [this_ct_cell_type]

# parse position file: update cells_dict, now all values contain that cell's x position as the second element and y position as the third element
# ex: {'CGATCG':['monocyte', '-6', '7'], 'CGATATC':['tCell', '10', '3']}
# note: all list elements are string type for now
with open(pos_file, "r") as f:
    pos_contents = f.read()
    pos_contents = pos_contents.split('\n')
pos_contents.pop()
for line in range(len(pos_contents)):
    pos_contents[line] = pos_contents[line].split('\t')
for cell in pos_contents:
    this_pos_cell_barcode = cell[0]
    this_pos_cell_x_crd = cell[1]
    this_pos_cell_y_crd = cell[2]
    cells_dict[this_pos_cell_barcode].append(this_pos_cell_x_crd)
    cells_dict[this_pos_cell_barcode].append(this_pos_cell_y_crd)

# plot points in left panel, build data structures to calculate coordinates of text legends in left panel
# these dictionaries are {'barcode':[4.0, 3,0]} notice now coords are float
monocytes_x_crds = []
monocytes_y_crds = []
tCells_x_crds = []
tCells_y_crds = []
bCells_x_crds = []
bCells_y_crds = []
for cell in cells_dict.keys():
    this_cell_type = cells_dict[cell][0]
    this_cell_x_pos = float(cells_dict[cell][1])
    this_cell_y_pos = float(cells_dict[cell][2])
    if this_cell_type == 'monocyte':
        dot_clr = (180/255, 180/255, 180/255)
        monocytes_x_crds.append(this_cell_x_pos)
        monocytes_y_crds.append(this_cell_y_pos)
    elif this_cell_type == 'tCell':
        dot_clr = (120/255, 172/255, 145/255)
        tCells_x_crds.append(this_cell_x_pos)
        tCells_y_crds.append(this_cell_y_pos)
    elif this_cell_type == 'bCell':
        dot_clr = (88/255, 85/255, 120/255)
        bCells_x_crds.append(this_cell_x_pos)
        bCells_y_crds.append(this_cell_y_pos)
    left_panel.plot(this_cell_x_pos, this_cell_y_pos, marker='o', markersize=4, markeredgewidth=0.1, markeredgecolor='black', markerfacecolor=dot_clr, linewidth=0)

# plot text in left panel
monocyte_label_x_pos = np.median(monocytes_x_crds)
monocyte_label_y_pos = np.median(monocytes_y_crds)
monocyte_label = left_panel.text(monocyte_label_x_pos, monocyte_label_y_pos, 'monocyte', ha='center', va='center')
monocyte_label.set_path_effects([patheffects.withStroke(linewidth=1, foreground='white')])
tCells_label_x_pos = np.median(tCells_x_crds)
tCells_label_y_pos = np.median(tCells_y_crds)
tCells_label = left_panel.text(tCells_label_x_pos, tCells_label_y_pos, 'tCell', ha='center', va='center')
tCells_label.set_path_effects([patheffects.withStroke(linewidth=1, foreground='white')])
bCells_label_x_pos = np.median(bCells_x_crds)
bCells_label_y_pos = np.median(bCells_y_crds)
bCells_label = left_panel.text(bCells_label_x_pos, bCells_label_y_pos, 'bCell', ha='center', va='center')
bCells_label.set_path_effects([patheffects.withStroke(linewidth=1, foreground='white')])

# make viridis color gradient
colors = np.array([
    (253/255, 231/255, 37/255),   # viridis5 (yellow)
    (94/255, 201/255, 98/255),    # viridis4
    (33/255, 145/255, 140/255),   # viridis3
    (59/255, 82/255, 139/255),    # viridis2
    (68/255, 1/255, 84/255)       # viridis1 (purple)
])
def interpolate_colormap(colors, n=100):
    # 4 segments between 5 control points
    segment_positions = np.linspace(0, 1, len(colors))
    # target positions for 100 colors
    t = np.linspace(0, 1, n)
    result = np.zeros((n, 3))
    for i in range(3):  # RGB channels
        result[:, i] = np.interp(t, segment_positions, colors[:, i])
    return [tuple(rgb) for rgb in result]
# generate list of 100 RGB values (first color is yellow, last is purple)
viridis_100 = interpolate_colormap(colors, 100)
viridis_100.reverse()

# plot middle panel points and text
cells_list = list(cells_dict.values())
x_scale = panel_width / (30 - (-30))    # 1.5 in / 60 data units = 0.025 in/unit
y_scale = panel_height / (30 - (-40))   # 1.5 in / 70 data units ≈ 0.02143 in/unit
for cell_1 in cells_list:
    x1 = float(cell_1[1])
    y1 = float(cell_1[2])
    nearby_cells = 0
    for cell_2 in cells_list:
        if cell_1 is cell_2:    # skip self-comparison
            continue
        x2 = float(cell_2[1])
        y2 = float(cell_2[2])
        distance = np.sqrt(((x2 - x1) * x_scale) ** 2 + ((y2 - y1) * y_scale) ** 2) # distance converted to inches
        if distance < 8 / 72: # 'nearby' threshold that makes colors match answer is ratio 8 / 72
            nearby_cells += 1
    if nearby_cells > 99:
        nearby_cells = 99    # cap at 99, not 100
    middle_panel.plot(x1, y1, marker='o', markersize=4, markeredgewidth=0, markerfacecolor=viridis_100[nearby_cells], linewidth=0)
middle_panel.text(-20.5, -35, 'Density', ha='center', va='center')

# fill right panel
x_min, x_max = right_panel.get_xlim()
rect_height = (y_max - y_min) / len(viridis_100)
rect_width = x_max - x_min
for i, color in enumerate(viridis_100):
    rectangle = mplpatches.Rectangle(
        (x_min, y_min + i * rect_height), width=rect_width, height=rect_height, facecolor=color, linewidth=0)
    right_panel.add_patch(rectangle)
    
# save figure
plt.savefig(out_file, dpi=600)