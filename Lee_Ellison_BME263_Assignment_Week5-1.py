# 263 run command:
# python3 Lee_Ellison_BME263_Assignment_Week5.py -b BME163_Input_Data_Week5.bed -g Mus_musculus.GRCm38.dna.primary_assembly.fa -o Lee_Ellison_BME263_Assignment_Week5.png -A A.png -T T.png -G G.png -C C.png

# ── Imports (only the same libraries used previously) ─────────────────────────
import matplotlib.pyplot as plt
import numpy as np
import argparse
import matplotlib.patches as mplpatches

# ══════════════════════════════════════════════════════════════════════════════
# COMMAND-LINE ARGUMENTS
# ══════════════════════════════════════════════════════════════════════════════
parser = argparse.ArgumentParser()
parser.add_argument('--bed_file',    '-b', type=str, action='store',
                    help='BED file of splice-junction coordinates')
parser.add_argument('--genome_file', '-g', type=str, action='store',
                    help='Genome FASTA (Mus_musculus GRCm38 primary assembly)')
parser.add_argument('--output_file', '-o', type=str, action='store',
                    help='Output PNG figure path',
                    default='Lee_Ellison_BME263_Assignment_Week5.png')
parser.add_argument('--A_file', '-A', type=str, action='store',
                    help='PNG image of the letter A', default='A.png')
parser.add_argument('--C_file', '-C', type=str, action='store',
                    help='PNG image of the letter C', default='C.png')
parser.add_argument('--G_file', '-G', type=str, action='store',
                    help='PNG image of the letter G', default='G.png')
parser.add_argument('--T_file', '-T', type=str, action='store',
                    help='PNG image of the letter T', default='T.png')
args = parser.parse_args()

bed_file    = args.bed_file
genome_file = args.genome_file
out_file    = args.output_file

# ══════════════════════════════════════════════════════════════════════════════
# LOAD LETTER IMAGES
# ══════════════════════════════════════════════════════════════════════════════
# plt.imread returns an H×W×4 RGBA numpy array.
# Each PNG is a coloured letter on a transparent background.
letter_images = {
    'A': plt.imread(args.A_file),
    'T': plt.imread(args.T_file),
    'G': plt.imread(args.G_file),
    'C': plt.imread(args.C_file),
}

# ══════════════════════════════════════════════════════════════════════════════
# LOAD GENOME FASTA → dictionary  { 'chr1': 'ACGT...', 'chr2': ..., ... }
# ══════════════════════════════════════════════════════════════════════════════
#
# The FASTA headers look like:
#   >1 dna:chromosome chromosome:GRCm38:1:1:195471971:1 REF
#   >MT dna:chromosome chromosome:GRCm38:MT:1:16299:1 REF
#
# The BED file already uses the "chr" prefix (chr1, chrM, chrX …).
# We therefore convert the raw FASTA chromosome names at parse time:
#   "1"  → "chr1"
#   "MT" → "chrM"
#   "X"  → "chrX"   (same pattern: prepend "chr", but MT→M is the exception)
#
# Strategy: read line by line, accumulating sequence lines in a list.
# When a new ">" header is encountered, join the accumulated list into a
# single string and store it, then start a fresh accumulator for the new chr.
# Joining once at the end of each chromosome is far faster than repeated
# string concatenation (strings are immutable in Python; += would be O(n²)).

print("Loading genome… (this may take a minute for the full mouse assembly)")

genome = {}          # { chr_name: sequence_string }
current_chr  = None  # name of the chromosome being assembled
current_lines = []   # lines of sequence accumulated for current_chr

def fasta_chr_to_bed_chr(raw_name):
    """
    Convert a raw FASTA chromosome name to the BED-file "chr" convention.

    Rules (from the FASTA headers seen in the assembly):
      "MT"          → "chrM"    (mitochondrial chromosome special case)
      anything else → "chr" + raw_name  (e.g. "1" → "chr1", "X" → "chrX")
    """
    if raw_name == 'MT':
        return 'chrM'
    return 'chr' + raw_name

with open(genome_file) as fh:
    for line in fh:
        line = line.rstrip('\n')   # remove newline only; keep any trailing spaces
        if line.startswith('>'):
            # ── Save the chromosome we just finished ──────────────────────
            if current_chr is not None:
                genome[current_chr] = ''.join(current_lines).upper()
                current_lines = []

            # ── Parse the new chromosome name from the header ─────────────
            # Header format: ">1 dna:chromosome ..."
            # We only need the first word after ">".
            raw_name    = line[1:].split()[0]          # e.g. "1", "MT", "X"
            current_chr = fasta_chr_to_bed_chr(raw_name)  # e.g. "chr1", "chrM"

        else:
            # Sequence line — append to current accumulator
            current_lines.append(line)

# Save the very last chromosome (no following ">" to trigger the save above)
if current_chr is not None and current_lines:
    genome[current_chr] = ''.join(current_lines).upper()

print(f"Loaded {len(genome)} chromosomes/scaffolds.")

# ══════════════════════════════════════════════════════════════════════════════
# REVERSE-COMPLEMENT HELPER
# ══════════════════════════════════════════════════════════════════════════════
# A complement table: every base maps to its Watson-Crick partner.
# We then reverse the string to get the reverse complement.
#
# This is needed for 3' splice sites (see extraction section below).
_COMPLEMENT = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence string."""
    return seq.translate(_COMPLEMENT)[::-1]

# ══════════════════════════════════════════════════════════════════════════════
# PARSE BED FILE AND EXTRACT SPLICE-JUNCTION SEQUENCES
# ══════════════════════════════════════════════════════════════════════════════
#
# BED columns (tab-separated, 0-based coordinates):
#   col 0 : chromosome  (already "chr1", "chrM" etc.)
#   col 1 : junction position  (start; same as col 2 for point features)
#   col 2 : junction position  (end;   same as col 1)
#   col 3 : name, e.g. "5'_172505320_..." or "3'_125405983_..."
#   col 5 : score (not used here)
#
# The junction position (col 1) is the 0-based genomic coordinate of the
# first base AFTER the exon (for 5' SS) or the first base OF the exon
# (for 3' SS). In both cases we want the 10 bp immediately upstream and
# the 10 bp immediately downstream, giving a 20 bp window:
#
#   window = genome[chrom][ pos-10 : pos+10 ]
#            ↑ 10 bases before           ↑ 10 bases after
#
# Mapped to our x-axis (x = i − 10):
#   i = 0  → x = −10  (10 bp upstream of junction)
#   i = 9  → x =  −1
#   i = 10 → x =   0  (the junction itself)
#   i = 19 → x =  +9
#
# WHY REVERSE-COMPLEMENT ALL 3' SEQUENCES?
# ─────────────────────────────────────────
# A 3' splice site (acceptor) lies at the end of an intron on whichever
# strand the gene is transcribed.  For genes on the minus (−) strand the
# intron runs in the direction of decreasing genomic coordinates, so the
# "upstream" (intronic) bases are actually to the RIGHT of `pos` in the
# genome string.  Extracting genome[pos-10:pos+10] and then taking the
# reverse complement flips the orientation so that:
#   • positions −10 … −1  always show the intronic side  (ending in AG)
#   • positions   0 …  +9 always show the exonic side
# This normalisation makes all 3' junctions directly comparable regardless
# of which strand they are on.

WINDOW = 10   # bases on each side of the junction

five_prime_seqs  = []
three_prime_seqs = []

with open(bed_file) as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue

        cols  = line.split('\t')
        chrom = cols[0]                # e.g. "chr1"
        pos   = int(cols[1])           # 0-based junction position
        name  = cols[3]                # e.g. "5'_172505320_..."

        # Skip chromosomes not present in our loaded genome
        # (e.g. scaffolds that might appear in BED but not in our dict)
        if chrom not in genome:
            continue

        chrom_seq = genome[chrom]
        chrom_len = len(chrom_seq)

        # Make sure the window doesn't fall outside the chromosome
        if pos - WINDOW < 0 or pos + WINDOW > chrom_len:
            continue

        # Extract the raw 20 bp window (0-based Python slicing)
        raw_seq = chrom_seq[pos - WINDOW : pos + WINDOW]

        # Classify and store
        if name.startswith("5'"):
            five_prime_seqs.append(raw_seq)

        elif name.startswith("3'"):
            # Reverse-complement to normalise orientation (see explanation above)
            three_prime_seqs.append(reverse_complement(raw_seq))

print(f"5' splice sites: {len(five_prime_seqs)}")
print(f"3' splice sites: {len(three_prime_seqs)}")

# ══════════════════════════════════════════════════════════════════════════════
# INFORMATION-CONTENT CALCULATION
# ══════════════════════════════════════════════════════════════════════════════
def compute_logo_stacks(sequences):
    """
    For every column in the alignment compute per-base drawing heights (bits).

    Math (from the Schneider & Stephens sequence-logo paper):
    ──────────────────────────────────────────────────────────
    1.  f_{b,i} = count(b at column i) / n          relative frequency
    2.  H_i     = −∑ f_{b,i} · log₂(f_{b,i})       Shannon entropy
    3.  e_n     = (s−1) / (2·n·ln2)                 small-sample correction
                  s = 4 (nucleotides), n = # sequences
    4.  R_i     = log₂(4) − (H_i + e_n)  [≥ 0]     info content (total stack height)
    5.  height  = f_{b,i} · R_i                     per-base drawing height

    Returns
    ───────
    stacks : list of dicts  (one per alignment column)
        'x'            : data-coord x value (column index − 10)
        'sorted_bases' : bases LEAST → MOST frequent
                         (draw bottom→top so most-frequent ends on top)
        'heights'      : {base: height_in_bits}
    """
    bases = ['A', 'C', 'G', 'T']
    n     = len(sequences)
    s     = 4

    # Small-sample correction: reduces apparent conservation for small n
    e_n = (s - 1) / (2 * n * np.log(2))

    seq_len = len(sequences[0])   # 20 bp
    stacks  = []

    for i in range(seq_len):

        # Count each base at column i
        counts = {b: 0 for b in bases}
        for seq in sequences:
            b = seq[i]
            if b in counts:
                counts[b] += 1

        total = sum(counts.values())
        freqs = {b: counts[b] / total for b in bases}

        # Shannon entropy
        H_i = 0.0
        for b in bases:
            if freqs[b] > 0:
                H_i -= freqs[b] * np.log2(freqs[b])

        # Information content (stack height, 0–2 bits)
        R_i = max(0.0, np.log2(s) - (H_i + e_n))

        # Per-base drawing heights
        heights = {b: freqs[b] * R_i for b in bases}

        # Sort ascending by frequency so least-frequent is drawn first (bottom)
        sorted_bases = sorted(bases, key=lambda b: freqs[b])

        # i=0 → x=−10, i=10 → x=0 (splice junction), i=19 → x=+9
        x = i - WINDOW

        stacks.append({'x': x, 'sorted_bases': sorted_bases, 'heights': heights})

    return stacks


# ══════════════════════════════════════════════════════════════════════════════
# LOGO DRAWING
# ══════════════════════════════════════════════════════════════════════════════
def draw_logo(panel, sequences, letter_images):
    """
    Render a sequence logo on *panel* using RGBA letter images.

    Placement:
      LEFT-ALIGNED: each letter image occupies [x, x+1] horizontally,
      so its left edge is anchored at the integer tick x.

    Stacking:
      Letters are drawn bottom→top (least→most frequent).
      `y_bottom` is a running cursor that advances by each letter's height.

    Vertical junction line:
      axvline(x=0) marks the splice boundary; drawn at zorder=3 (above
      letter images at zorder=2) so it is never hidden.
    """
    stacks = compute_logo_stacks(sequences)

    for stack in stacks:
        x            = stack['x']
        sorted_bases = stack['sorted_bases']
        heights      = stack['heights']
        y_bottom     = 0.0

        for base in sorted_bases:
            h = heights[base]
            if h < 1e-9:
                continue

            # extent = [x_left, x_right, y_bottom, y_top] in data coordinates.
            # [x, x+1] gives left-alignment (left edge at tick x).
            panel.imshow(
                letter_images[base],
                extent=[x, x + 1, y_bottom, y_bottom + h],
                aspect='auto',
                zorder=2,
                interpolation='bilinear',
            )
            y_bottom += h

    # Vertical line at the splice junction (x = 0), on top of letters
    panel.axvline(x=0, color='black', linewidth=0.5, zorder=3)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE, PANELS, AND RENDERING
# ══════════════════════════════════════════════════════════════════════════════
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
panel2.set_yticks([])       # y-scale already shown on panel1
panel2.set_xlim(-10, 10)
panel2.set_title("3'SS")
panel2.set_xlabel('Distance to\nSplice Site')

# Draw logos on both panels
draw_logo(panel1, five_prime_seqs,  letter_images)
draw_logo(panel2, three_prime_seqs, letter_images)

# Save figure
plt.savefig(out_file, dpi=600)
print(f"Saved figure to {out_file}")