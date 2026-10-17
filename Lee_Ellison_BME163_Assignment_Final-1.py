# run command used
'''
Lee_Ellison_BME163_Assignment_Final.py -f1 BME163_Input_Data_Final_1.psl -f2 BME163_Input_Data_Final_2.psl -g gencode.vM12.annotation.gtf -c chr7:45232000-45241000 -o Lee_Ellison_BME163_Assignment_Final.png
'''

# import libraries
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mplpatches
import argparse

# function to parse psl files into a data structure suitable for panel plotting
def parse_psl_file(file_path, chromosome, left, right, sort_order, contained_only):
    """
    Takes .PSL file and region of interest as input then outputs a list of tuples containg all info needed to
    fill the corresponding psl file panel.
    """

    # get file contents into file_contents as list of nested lists, each nested list is 1 line of the file, each element in a nested list is the item in 1 column of that file line
    with open(file_path, 'r') as f:
        file_contents = f.read().strip()
    file_contents = file_contents.split('\n')
    for i in range(len(file_contents)):
        file_contents[i] = file_contents[i].strip().split('\t')
    
    # filter to keep only reads within the chromosome and roi specified in command arguments
    reads_of_interest = []
    for read in file_contents:
        if read[13] != chromosome:
            continue
        this_read_start = int(read[15])
        this_read_end   = int(read[16])
        if this_read_start >= right or this_read_end <= left:
            continue
        if contained_only and (this_read_start <= left or this_read_end >= right):
            continue
        this_read_blockstarts = [int(x) for x in read[20].rstrip(',').split(',')]
        this_read_blocksizes  = [int(x) for x in read[18].rstrip(',').split(',')]

        if not any(b_start + b_size > left and b_start < right for b_start, b_size in zip(this_read_blockstarts, this_read_blocksizes)):
            continue

        reads_of_interest.append((this_read_start,
                                  this_read_end,
                                  this_read_blockstarts,
                                  this_read_blocksizes))
    
    # sort by start/end as specified
    if sort_order == 'by_start':
        reads_of_interest.sort(key=lambda r: r[0])
    elif sort_order == 'by_end':
        reads_of_interest.sort(key=lambda r: r[1])

    # sorted list of tuples used to plot psl panels
    return reads_of_interest

def draw_psl_panel(panel, reads, color, line_width, bottom_pad, top_pad_inches):
    """
    Draws sorted reads onto the panel using greedy interval packing.
    Each read: a thin connector (height 0.05) spanning start->end, plus one rectangle (height 0.5) per block.
    """
    track_ends = []  # track_ends[i] = rightmost x currently used on track i

    for start, end, blockstarts, blocksizes in reads:

        # find first track this read fits on
        y = None
        for i, track_end in enumerate(track_ends):
            if track_end < start:
                track_ends[i] = end
                y = i
                break
        if y is None:
            y = len(track_ends)
            track_ends.append(end)

        # thin connector across the full read span
        panel.add_patch(mplpatches.Rectangle(
            (start, y - 0.05/2),
            end - start,
            0.05,
            facecolor=color,
            edgecolor='none',
            linewidth=line_width,
        ))

        # one rectangle per block
        for b_start, b_size in zip(blockstarts, blocksizes):
            panel.add_patch(mplpatches.Rectangle(
                (b_start, y - 0.5/2),
                b_size,
                0.5,
                facecolor=color,
                edgecolor='none',
                linewidth=line_width,
            ))

    # ylim: pin top track to top_pad_inches below the panel top
    N = len(track_ends)
    top_track_edge = N - 0.75
    ymin = -0.25 - bottom_pad

    panel_h_inches = panel.get_position().height * panel.figure.get_figheight()
    r = top_pad_inches / panel_h_inches
    ymax = (top_track_edge - r * ymin) / (1 - r)

    panel.set_ylim(ymin, ymax)

def parse_gtf_file(file_path, chromosome, left, right, sort_order):
    """
    Takes .GTF file and region of interest as input, returns a sorted list of
    transcripts whose features overlap the ROI. Each transcript is a tuple:
        (start, end, exons, cds_list)
    where exons and cds_list are lists of (start, end) tuples.
    """

    # get file contents as list of lists
    with open(file_path, 'r') as f:
        file_contents = f.read().strip()
    file_contents = file_contents.split('\n')
    for i in range(len(file_contents)):
        file_contents[i] = file_contents[i].strip().split('\t')

    # remove comment lines from file_contents
    file_contents = file_contents[5:]

    # group exon/CDS rows by transcript_id, keeping only those in the ROI
    transcripts = {}  # tid -> {'exons': [(s,e), ...], 'cds': [(s,e), ...]}
    for row in file_contents:
        if row[0] != chromosome:
            continue
        feature_type = row[2]
        if feature_type != 'exon' and feature_type != 'CDS':
            continue
        feature_start = int(row[3])
        feature_end   = int(row[4])

        # extract transcript_id from the attributes column (col 8)
        attrs = row[8]
        tid_start = attrs.find('transcript_id "') + len('transcript_id "')
        tid_end   = attrs.find('"', tid_start)
        tid = attrs[tid_start:tid_end]

        if tid not in transcripts:
            transcripts[tid] = {'exons': [], 'cds': []}
        if feature_type == 'exon':
            transcripts[tid]['exons'].append((feature_start, feature_end))
        else:  # CDS
            transcripts[tid]['cds'].append((feature_start, feature_end))

    # build the list of transcripts to draw; compute span from exon bounds
    transcripts_of_interest = []
    for tid, data in transcripts.items():
        if not data['exons']:
            continue
        t_start = min(e[0] for e in data['exons'])
        t_end   = max(e[1] for e in data['exons'])
        if t_start >= right or t_end <= left:
            continue
        transcripts_of_interest.append((t_start, t_end, data['exons'], data['cds']))

    # sort
    if sort_order == 'by_start':
        transcripts_of_interest.sort(key=lambda t: t[0])
    elif sort_order == 'by_end':
        transcripts_of_interest.sort(key=lambda t: t[1])

    return transcripts_of_interest

def draw_gtf_panel(panel, transcripts, color, line_width, bottom_pad, top_pad_inches):
    """
    Draws sorted transcripts onto the panel using greedy interval packing.
    Each transcript: connector (height 0.05), exons (height 0.25),
    and CDS rectangles (height 0.5) layered on top.
    """
    track_ends = []  # track_ends[i] = rightmost x currently used on track i

    for start, end, exons, cds_list in transcripts:

        # find first track this transcript fits on
        y = None
        for i, track_end in enumerate(track_ends):
            if track_end < start:
                track_ends[i] = end
                y = i
                break
        if y is None:
            y = len(track_ends)
            track_ends.append(end)

        # thin connector across the full transcript span
        panel.add_patch(mplpatches.Rectangle(
            (start, y - 0.05/2),
            end - start,
            0.05,
            facecolor=color,
            edgecolor='black',
            linewidth=line_width,
        ))

        # exon rectangles (height 0.25)
        for e_start, e_end in exons:
            panel.add_patch(mplpatches.Rectangle(
                (e_start, y - 0.25/2),
                e_end - e_start,
                0.25,
                facecolor=color,
                edgecolor='black',
                linewidth=line_width,
            ))

        # CDS rectangles (height 0.5) drawn on top of exons
        for c_start, c_end in cds_list:
            panel.add_patch(mplpatches.Rectangle(
                (c_start, y - 0.5/2),
                c_end - c_start,
                0.5,
                facecolor=color,
                edgecolor='black',
                linewidth=line_width,
            ))

    # ylim: pin top track to top_pad_inches below the panel top
    N = len(track_ends)
    top_track_edge = N - 0.75
    ymin = -0.25 - bottom_pad

    panel_h_inches = panel.get_position().height * panel.figure.get_figheight()
    r = top_pad_inches / panel_h_inches
    ymax = (top_track_edge - r * ymin) / (1 - r)

    panel.set_ylim(ymin, ymax)

def draw_histogram_panel(panel, reads, color, left, right):
    """
    Draws a 1bp-bin read coverage histogram on the panel.
    Only block (aligned) regions contribute to coverage; connectors don't.
    Bars are drawn without an outline per spec.
    """
    # one bin per bp across the ROI
    coverage = np.zeros(right - left, dtype=int)

    # accumulate coverage from every block of every read
    for start, end, blockstarts, blocksizes in reads:
        for b_start, b_size in zip(blockstarts, blocksizes):
            b_end = b_start + b_size
            # clip block to the ROI before indexing
            clipped_start = max(b_start, left)
            clipped_end   = min(b_end, right)
            if clipped_start < clipped_end:
                coverage[clipped_start - left : clipped_end - left] += 1

    # plot as 1bp-wide bars, no outline
    x_positions = np.arange(left, right)
    panel.bar(x_positions, coverage, width=1, color=color,
              linewidth=0, edgecolor='none')

    # ylim fits the data tightly (panel 0 abuts panel 1, so peaks reach the boundary)
    panel.set_ylim(max(coverage.max(), 1), 0)

def make_figure():
    """
    """

    # command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_file'    ,'-o'   ,type=str   ,action='store'     ,help='output file'                       ,default='Lee_Ellison_BME163_Assignment_Final.png')
    parser.add_argument('--psl_file_1'  ,'-f1'  ,type=str   ,action='store'     ,help='input psl file 1'                  ,default='BME163_Input_Data_Final_1.psl')
    parser.add_argument('--psl_file_2'  ,'-f2'  ,type=str   ,action='store'     ,help='input psl file 2'                  ,default='BME163_Input_Data_Final_2.psl')
    parser.add_argument('--gtf_file'    ,'-g'   ,type=str   ,action='store'     ,help='input gencode annotations file'    ,default='gencode.vM12.annotation.gtf')
    parser.add_argument('--roi'         ,'-c'   ,type=str   ,action='store'     ,help='input region of interest'          ,default='chr7:45232000-45241000')
    args = parser.parse_args()
    out_file   = args.out_file
    psl_file_1 = args.psl_file_1
    psl_file_2 = args.psl_file_2
    gtf_file   = args.gtf_file
    roi        = args.roi

    # process the -c argument into the 3 pieces of info we want from it
    chromosome = roi[0                  :   roi.find(':')]      # str
    left       = roi[roi.find(':')+1    :   roi.find('-')]      # int
    right      = roi[roi.find('-')+1    :                ]      # int
    left = int(left)
    right = int(right)

    # generate figure
    plt.style.use('BME163')
    figure_height = 6
    figure_width  = 5
    fig_1 = plt.figure(figsize=(figure_width, figure_height))

    # panel_0: bottom-most panel, histogram fused to bottom of panel 1
    panel_0 = plt.axes([0.1/figure_width, 0.1/figure_height, 4/figure_width, 0.4/figure_height], frameon=True)
    panel_0.set_xlim(left, right)
    panel_0.set_xticks([])
    panel_0.set_yticks([])

    # panel_1: panel above panel_0, blue
    panel_1 = plt.axes([0.1/figure_width, 0.5/figure_height, 4/figure_width, 1.5/figure_height], frameon=True)
    panel_1.set_xlim(left, right)
    panel_1.set_xticks([])
    panel_1.set_yticks([])

    # panel_2: panel above panel_1, orange
    panel_2 = plt.axes([0.1/figure_width, 2.2/figure_height, 4/figure_width, 1.5/figure_height], frameon=True)
    panel_2.set_xlim(left, right)
    panel_2.set_xticks([])
    panel_2.set_yticks([])

    # panel_3: top-most panel, made from gtf file, gray
    panel_3 = plt.axes([0.1/figure_width, 3.9/figure_height, 4/figure_width, 1.5/figure_height], frameon=True)
    panel_3.set_xlim(left, right)
    panel_3.set_xticks([])
    panel_3.set_yticks([])

    # fill panels: 0 and 1 blue, 2 orange
    panel_2_data        = parse_psl_file(psl_file_1, chromosome, left, right, sort_order='by_end', contained_only=True)
    panel_0_and_1_data  = parse_psl_file(psl_file_2, chromosome, left, right, sort_order='by_start', contained_only=False)
    iOrange = (230/255, 87/255, 43/255)
    iBlue   = (88/255, 85/255, 120/255)
    draw_psl_panel(panel_2, panel_2_data,       iOrange,    0,    bottom_pad=0.75,      top_pad_inches=0.1389)
    draw_psl_panel(panel_1, panel_0_and_1_data, iBlue,      0.05, bottom_pad=0.76,    top_pad_inches=0.1514)
    panel_3_data = parse_gtf_file(gtf_file, chromosome, left, right, sort_order='by_start')
    draw_gtf_panel(panel_3, panel_3_data, 'Grey', 0.25, bottom_pad=0.75, top_pad_inches=0.237)
    draw_histogram_panel(panel_0, panel_0_and_1_data, iBlue, left, right)

    # save figure
    plt.savefig(out_file, dpi=2400)
    print('done')
    
make_figure()