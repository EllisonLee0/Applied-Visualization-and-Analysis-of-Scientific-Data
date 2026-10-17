# run command
# python3 Lee_Ellison_BME163_Assignment_Week2.py -i BME163_Input_Data_1.txt -o Lee_Ellison_BME163_Assignment_Week2.png

# import libraries
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse

# create commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument('--inputFile', '-i', type=str, action='store', help='input file')
parser.add_argument('--outFile', '-o', type=str, action='store', help='output file')
args = parser.parse_args()
in_file = args.inputFile
out_file = args.outFile

# create figure: 3in x 3in, BME163 style sheet
plt.style.use('BME163')
figure_width = 3
figure_height = 3
plt.figure(figsize=(figure_width, figure_height))

# create panels
# main: 1.5in x 1.5in, tick marks bottom 0 - 15
main_left = 0.7
main_bottom = 0.3
main_width = 1.5
main_height = 1.5
main_panel = plt.axes([main_left / figure_width,
                       main_bottom / figure_height,
                       main_width / figure_width,
                       main_height / figure_height])
main_panel.tick_params(bottom=True,labelbottom=True,
                      left=False,labelleft=False,
                      right=False,labelright=False,
                      top=False,labeltop=False)
main_panel.set_xlim(0,15)
main_panel.set_ylim(0,15)
# left: 0.25in wide, 1.5in high, ticks left 0-15, ticks bottom 20-0
left_left = 0.38
left_bottom = 0.3
left_width = 0.25
left_height = 1.5
left_panel = plt.axes([left_left / figure_width,
                       left_bottom / figure_height,
                       left_width / figure_width,
                       left_height / figure_height],
                       zorder=2)
left_panel.tick_params(bottom=True,labelbottom=True,
                      left=True,labelleft=True,
                      right=False,labelright=False,
                      top=False,labeltop=False)
left_panel.set_xlim(20,0)
left_panel.set_ylim(0,15)
left_panel.set_yticks(np.arange(0, 16, 5)) # left ticks every 5, not every 2.5 which is default for some reason
# top: 1.5in wide, 0.25in high, 1.5in high, ticks left 0-20
top_left = 0.7
top_bottom = 1.87
top_width = 1.5
top_height = 0.25
top_panel = plt.axes([top_left / figure_width,
                      top_bottom / figure_height,
                      top_width / figure_width,
                      top_height / figure_height],
                      zorder=2)
top_panel.tick_params(bottom=False,labelbottom=False,
                      left=True,labelleft=True,
                      right=False,labelright=False,
                      top=False,labeltop=False)
top_panel.set_xlim(0,15)
top_panel.set_ylim(0,20)
# note: numbers on axes tick marks are slightly not as tall as answer key right now, fix later

# read input file data
with open(in_file, "r") as f:
    file_contents = f.read()
file_lines = file_contents.split("\n") # each list item is 1 line of file
file_line_items = []
for line in file_lines:
    line_items_list = line.split("\t")
    file_line_items.append(line_items_list) # each list item is a nested list with 3 items, each nested list item is 1 tab delimited item in the file
file_line_items.pop() # get rid of last nested list to prevent indexing errors
x_vals = []
y_vals = []
for line in range(len(file_line_items)):
    x_vals.append(file_line_items[line][1]) # each list item is the x value of that file line
    y_vals.append(file_line_items[line][2]) # each list item is the y value of that file line

# do +1 and log2 conversion of x and y coords
for num in range(len(x_vals)):
    x_vals[num] = float(x_vals[num])
    x_vals[num] += 1
    x_vals[num] = np.log2(x_vals[num])
    y_vals[num] = float(y_vals[num])
    y_vals[num] += 1
    y_vals[num] = np.log2(y_vals[num])

# plot points on main panel
main_panel.plot(x_vals,
                y_vals,
                marker='o',
                markersize=3,
                markerfacecolor=(88/255,85/255,120/255),
                markeredgewidth=0,
                linewidth=0,
                alpha=0.1
                )

# make top histogram
xmin = min(x_vals)
xmax = max(x_vals)
num_bins = 31
bin_ticks = np.linspace(0, 15, num_bins)
x_counts, bin_ticks = np.histogram(x_vals, bin_ticks)
for i in range(len(x_counts)-1):
    rectLeft = bin_ticks[i]
    rectBottom = 0
    rectWidth = bin_ticks[i+1] - bin_ticks[i]
    rectHeight = x_counts[i] + 1
    rectHeight = np.log2(rectHeight)
    rectangle=mplpatches.Rectangle((rectLeft, rectBottom),
                                   rectWidth,
                                   rectHeight,
                                   facecolor=(120/255,172/255,145/255),
                                   edgecolor='black',
                                   linewidth=0.3,
                                   zorder=1,
                                   )
    top_panel.add_patch(rectangle)

# make left histogram
ymin = min(y_vals)
ymax = max(y_vals)
num_bins = 31
bin_ticks = np.linspace(0, 15, num_bins)
y_counts, bin_ticks = np.histogram(y_vals, bin_ticks)
for i in range(len(y_counts)-1):
    rectHeight = bin_ticks[i+1] - bin_ticks[i]
    rectWidth = np.log2(y_counts[i] + 1)
    rectLeft = 0
    rectBottom = bin_ticks[i]
    rectangle = mplpatches.Rectangle((rectLeft, rectBottom),
                                     rectWidth,
                                     rectHeight,
                                     facecolor='grey',
                                     edgecolor='black',
                                     linewidth=0.3,
                                     zorder=1)
    left_panel.add_patch(rectangle)

# save figure
plt.savefig(out_file, dpi=600)