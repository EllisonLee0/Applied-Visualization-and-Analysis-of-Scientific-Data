# while developing this code, I always ran it using a conda environment where I had matplotlib installed, to achieve this:
    # conda create -n bme163env python=3.11
    # conda activate bme163env
    # conda install matplotlib

# terminal command to run code: python3 Lee_Ellison_BME163_Assignment_Week1.py -o Lee_Ellison_BME163_Assignment_Week1.png

import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import argparse

# set commandline argument: name of output file
parser = argparse.ArgumentParser()
parser.add_argument('--outFile','-o' ,type=str,action='store',help='output file')
args = parser.parse_args()
outFile = args.outFile

# create figure, 5in width, 2in height
width = 5
height = 2
plt.figure(figsize=(width,height))

# to install bme163 style sheet using terminal, I did:
    # print(mpl.get_configdir()) outputs: /Users/ellisonlee/.matplotlib
    # mkdir -p /Users/ellisonlee/.matplotlib/stylelib
    # mv /Users/ellisonlee/Downloads/BME163.mplstyle /Users/ellisonlee/.matplotlib/stylelib/

# use bme 163 style sheet without using absolute path
plt.style.use('BME163')

# create panel
left1 = 0.2 / width
bottom1 = 0.2 / height
width1 = 2 / width
height1 = 1 / height
panel1=plt.axes([left1,bottom1,width1,height1],zorder=3)

# set panel axes scales and ranges
panel1.set_xlim(0,101)
panel1.set_ylim(0,200)

# remove axes tick marks
panel1.set_xticks([])
panel1.set_yticks([])

# add horizontal line down middle of panel
panel1.plot(
    [0, 101],
    [100, 100],
    color='black',
    linewidth=1,
    zorder=3
)

# create gradient of colors seen in top rectangle
viridis5 = (253/255, 231/255, 37/255)
viridis4 = (94/255, 201/255, 98/255)
viridis3 = (33/255, 145/255, 140/255)
viridis2 = (59/255, 82/255, 139/255)
viridis1 = (68/255, 1/255, 84/255)
colors = [viridis1, viridis2, viridis3, viridis4, viridis5]
gradient = []
n_total = 101 # how many thin, vertical rectangles should be in the gradient
n_segments = len(colors) - 1  # the total gradient is not uniform, it is made of 4 distinct 'segments' aka sub-gradients defined between every 2 successive viridis rgb values given
for i in range(n_total):
    # figure out which segment of the total gradient we are in 
    t_global = i / (n_total - 1)
    segment = min(int(t_global * n_segments), n_segments - 1)
    t_local = (t_global * n_segments) - segment # figure out which rectangle of this current segment we are at
    # compute the rgb values for this rectangle
    c0 = colors[segment]
    c1 = colors[segment + 1]
    r = c0[0] + (c1[0] - c0[0]) * t_local
    g = c0[1] + (c1[1] - c0[1]) * t_local
    b = c0[2] + (c1[2] - c0[2]) * t_local
    gradient.append((r, g, b))

# draw top gradient using thin rectangles
n_top = len(gradient)
for i in range(n_top):
    width = 1 / n_top
    x_center = (i + 0.5) / n_top
    rectangle1 = mplpatches.Rectangle(
    (x_center * 101 - (101/n_top)/2, 100),
    101 / n_top,
    100,
    color=gradient[i],
    zorder=2
    )
    panel1.add_patch(rectangle1)

# create gradient of colors seen in bottom diagonal
plasma5 = (237/255, 252/255, 27/255)
plasma4 = (245/255, 135/255, 48/255)
plasma3 = (190/255, 48/255, 101/255)
plasma2 = (87/255, 0/255, 151/255)
plasma1 = (15/255, 0/255, 118/255)
gradient = []
colors = [plasma1, plasma2, plasma3, plasma4, plasma5]
n_segments = len(colors) - 1
n = 101
for i in range(n):
    t_global = i / (n - 1)
    segment = min(int(t_global * n_segments), n_segments - 1)
    t_local = (t_global * n_segments) - segment
    c0 = colors[segment]
    c1 = colors[segment + 1]
    r = c0[0] + (c1[0] - c0[0]) * t_local
    g = c0[1] + (c1[1] - c0[1]) * t_local
    b = c0[2] + (c1[2] - c0[2]) * t_local
    gradient.append((r, g, b))

# draw bottom gradient using left-edge aligned dots
n = 101
for i in range(n):
    x = i * (101 / n)
    y = x
    panel1.plot(
        x, y,
        marker='.',
        color=gradient[i],
        linestyle='None',
        zorder=1
    )

# save figure, 600 dpi
plt.savefig(outFile,dpi=600)