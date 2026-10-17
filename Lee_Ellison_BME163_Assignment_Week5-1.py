# 163 run command:
# python3 Lee_Ellison_BME163_Assignment_Week5.py -s BME163_Input_Data_Week5.fasta -o Lee_Ellison_BME163_Assignment_Week5.png -A A.png -T T.png -G G.png -C C.png

# ── Imports ───────────────────────────────────────────────────────────────────
import matplotlib.pyplot as plt
import numpy as np
import argparse
import matplotlib.patches as mplpatches

# ── Command-line arguments ────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument('--splice_file', '-s', type=str, action='store',
                    help='FASTA file of splice-site sequences',
                    default='BME163_Input_Data_Week5.fasta')
parser.add_argument('--output_file', '-o', type=str, action='store',
                    help='Output PNG figure path',
                    default='Lee_Ellison_BME163_Assignment_Week5.png')
parser.add_argument('--A_file', '-A', type=str, action='store',
                    help='PNG image of the letter A', default='A.png')
parser.add_argument('--C_file', '-C', type=str, action='store',
                    help='PNG image of the letter C', default='C.png')
parser.add_argument('--G_file', '-G', type=str, action='store',
                    help='PNG image of the letter G', default='G.png')
parser.add_argument('--T_file', '-T', type=str, action='store',
                    help='PNG image of the letter T', default='T.png')
args = parser.parse_args()

splice_file = args.splice_file
out_file    = args.output_file

# ── Load letter images ────────────────────────────────────────────────────────
# plt.imread returns an HxWx4 RGBA numpy array for PNG files.
# Each image is a coloured letter on a transparent background.
letter_images = {
    'A': plt.imread(args.A_file),
    'T': plt.imread(args.T_file),
    'G': plt.imread(args.G_file),
    'C': plt.imread(args.C_file),
}

# ── Parse FASTA and split by splice-site type ─────────────────────────────────
# Headers start with either ">5'" or ">3'".
# Every sequence is exactly 20 bp, representing the 10 bp on each side of
# the splice junction (positions -10 to -1 = upstream, 0 to +9 = downstream).
five_prime_seqs  = []   # sequences belonging to 5' splice sites
three_prime_seqs = []   # sequences belonging to 3' splice sites

current_header = None
current_seq    = ''

with open(splice_file) as fh:
    for line in fh:
        line = line.strip()
        if line.startswith('>'):
            # When we hit a new header, save the previous record first
            if current_header is not None and current_seq:
                if current_header.startswith(">5'"):
                    five_prime_seqs.append(current_seq)
                elif current_header.startswith(">3'"):
                    three_prime_seqs.append(current_seq)
            current_header = line
            current_seq    = ''
        else:
            current_seq += line.upper()   # accumulate bases; force uppercase

    # Save the very last record in the file (no trailing header to trigger it)
    if current_header is not None and current_seq:
        if current_header.startswith(">5'"):
            five_prime_seqs.append(current_seq)
        elif current_header.startswith(">3'"):
            three_prime_seqs.append(current_seq)

# ── Information-content calculation ──────────────────────────────────────────
def compute_logo_stacks(sequences):
    """
    For every column in the alignment, compute per-base drawing heights
    (in bits) for the sequence logo.

    Math summary
    ------------
    1.  f_{b,i}  = count of base b at position i  /  total sequences
                   (relative frequency)

    2.  H_i      = -sum( f_{b,i} * log2(f_{b,i}) )
                   (Shannon entropy; 0 = perfectly conserved, 2 = random)

    3.  e_n      = (s - 1) / (2 * n * ln2)
                   (small-sample correction; s=4 for nucleotides, n=# seqs)

    4.  R_i      = log2(4) - (H_i + e_n)   [clamped to >= 0]
                   (information content = total stack height in bits)

    5.  height   = f_{b,i} * R_i
                   (fraction of the stack occupied by base b)

    Returns
    -------
    stacks : list of dicts, one per alignment column
        'x'            : x-axis coordinate mapped as  i - 10  (so -10 to +9)
        'sorted_bases' : bases ordered LEAST to MOST frequent
                         (drawn bottom to top so the most-frequent letter
                          ends up on top, matching the expected figure)
        'heights'      : dict {base: drawing_height_in_bits}
    """
    bases = ['A', 'C', 'G', 'T']
    n     = len(sequences)      # number of sequences in the alignment
    s     = 4                   # nucleotide alphabet size

    # Small-sample correction: compensates for over-estimating conservation
    # in small datasets.  Formula: (s-1) / (2*n*ln2)
    e_n = (s - 1) / (2 * n * np.log(2))

    seq_len = len(sequences[0])   # 20 bp per sequence
    stacks  = []

    for i in range(seq_len):

        # Step 1 — count each base at this column
        counts = {b: 0 for b in bases}
        for seq in sequences:
            base = seq[i]
            if base in counts:
                counts[base] += 1

        total = sum(counts.values())
        freqs = {b: counts[b] / total for b in bases}   # relative frequencies

        # Step 2 — Shannon entropy
        H_i = 0.0
        for b in bases:
            if freqs[b] > 0:
                H_i -= freqs[b] * np.log2(freqs[b])

        # Step 3+4 — information content (total stack height, 0 to 2 bits)
        R_i = max(0.0, np.log2(s) - (H_i + e_n))

        # Step 5 — per-base drawing heights
        heights = {b: freqs[b] * R_i for b in bases}

        # Sort ascending by frequency:
        #   least-frequent base drawn first (at the bottom of the stack)
        #   most-frequent base drawn last   (at the top  of the stack)
        sorted_bases = sorted(bases, key=lambda b: freqs[b])

        # Map sequence index i (0-19) to x-axis value (-10 to +9):
        #   index 0  -> x = -10  (10 bp upstream of junction)
        #   index 10 -> x =   0  (splice junction)
        #   index 19 -> x =  +9
        x = i - 10

        stacks.append({
            'x':            x,
            'sorted_bases': sorted_bases,
            'heights':      heights,
        })

    return stacks


# ── Plotting helper ───────────────────────────────────────────────────────────
def draw_logo(panel, sequences, letter_images):
    """
    Draw a sequence logo on *panel* using pre-loaded RGBA letter images.

    Letter alignment — LEFT-ALIGNED
    --------------------------------
    Each letter image occupies the horizontal window  [x, x+1]  in data
    coordinates, where x is the integer position (-10, -9, ..., +9).
    The LEFT edge of the glyph sits exactly on the tick mark x; the glyph
    extends one full unit to the RIGHT.

    Why this gives left-alignment:
      In the previous (centred) version the extent was [x-0.5, x+0.5],
      meaning the glyph was symmetric around the tick.  Shifting the window
      right by 0.5 units — to [x, x+1] — anchors the left edge at the tick,
      which is the standard definition of left-alignment on a discrete axis.

    Vertical stacking
    -----------------
    A running `y_bottom` cursor starts at 0 and advances by each letter's
    height (in bits) as we draw from least-frequent (bottom) to
    most-frequent (top).

    Rendering
    ---------
    `imshow` with `extent` maps the image pixels into data coordinates.
    `aspect='auto'` lets matplotlib freely stretch the image to fill the
    requested rectangle.  The alpha channel in the RGBA PNG handles
    transparency automatically — no white boxes behind the letters.

    Vertical line at x = 0
    -----------------------
    `axvline(x=0)` draws a thin vertical line across the full height of the
    panel at the splice-junction position.  It is drawn at zorder=3, which
    is higher than the letter images (zorder=2), so the line sits visually
    on top of any letter that might overlap x=0.
    """
    stacks = compute_logo_stacks(sequences)

    for stack in stacks:
        x            = stack['x']
        sorted_bases = stack['sorted_bases']   # least to most frequent
        heights      = stack['heights']

        y_bottom = 0.0   # running baseline; moves up as we add letters

        for base in sorted_bases:
            h = heights[base]
            if h < 1e-9:       # skip bases whose height is negligibly small
                continue

            img = letter_images[base]   # RGBA numpy array

            # LEFT-ALIGNED placement:
            #   x_left  = x       (left edge sits on the integer tick)
            #   x_right = x + 1   (extends one full unit to the right)
            #   y range = [y_bottom, y_bottom + h]  in bits
            panel.imshow(
                img,
                extent=[x, x + 1, y_bottom, y_bottom + h],
                aspect='auto',
                zorder=2,                 # above axes background/grid
                interpolation='bilinear',
            )

            y_bottom += h   # raise baseline for the next (more-frequent) letter

    # ── Vertical line marking the splice junction at x = 0 ───────────────
    # axvline draws a line that spans the full y-range of the panel.
    # zorder=3 puts it on top of the letter images (zorder=2).
    panel.axvline(
        x=0,
        color='black',
        linewidth=0.5,
        zorder=3,
    )


# ── Figure and panels ─────────────────────────────────────────────────────────
plt.style.use('BME163')

figure_width  = 5
figure_height = 2
fig_1 = plt.figure(figsize=(figure_width, figure_height))

panel_width  = 1.5
panel_height = 0.5

# Left panel — 5' splice sites
panel1 = plt.axes(
    [0.5 / figure_width, 0.3, panel_width / figure_width, panel_height / figure_height],
    frameon=True,
)
panel1.set_ylim(0, 2)
panel1.set_xlim(-10, 10)
panel1.set_title("5'SS")
panel1.set_ylabel('Bits')
panel1.set_xlabel('Distance to\nSplice Site')

# Right panel — 3' splice sites
panel2 = plt.axes(
    [2.2 / figure_width, 0.3, panel_width / figure_width, panel_height / figure_height],
    frameon=True,
)
panel2.set_ylim(0, 2)
panel2.set_yticks([])       # hide y ticks — scale already shown on left panel
panel2.set_xlim(-10, 10)
panel2.set_title("3'SS")
panel2.set_xlabel('Distance to\nSplice Site')

# ── Draw logos ────────────────────────────────────────────────────────────────
draw_logo(panel1, five_prime_seqs,  letter_images)
draw_logo(panel2, three_prime_seqs, letter_images)

# ── Save ──────────────────────────────────────────────────────────────────────
plt.savefig(out_file, dpi=600)
print(f"Saved figure to {out_file}")