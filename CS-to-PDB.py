# ==============================================================================
# CS-to-PDB : Sequence Conservation to Protein Structure Mapper
# ==============================================================================
#
# OVERVIEW
# --------
# CS-to-PDB is a standalone application that calculates residue conservation
# from a Multiple Sequence Alignment (MSA) and maps evolutionary information
# onto three-dimensional protein structures by storing conservation scores in
# the PDB B-factor field.
#
# The generated structures can be visualized directly in PyMOL or ChimeraX
# using identical color palettes. The program also produces an interactive
# HTML alignment viewer to compare sequence conservation with structural data.
#
#
# FEATURES
# --------
# • Calculates residue conservation using multiple scoring methods:
#     - BLOSUM62 similarity
#     - Shannon Entropy
#     - Property Conservation
#     - Jensen-Shannon Divergence (JSD)
#     - Rate4Site-like weighted conservation
#
# • Optional gap penalty for poorly aligned regions.
#
# • Maps conservation scores onto PDB or mmCIF structures.
#
# • Supports all Matplotlib colormaps, including custom palettes
#   (e.g. MolNympheas), ensuring identical coloring across all outputs.
#
# • Generates:
#     - PDB structure with conservation values stored as B-factors
#     - Residue conservation table (.txt)
#     - PyMOL visualization script (.pml)
#     - ChimeraX visualization script (.cxc)
#     - Interactive HTML alignment viewer
#
#
# OUTPUT FILES
# ------------
# • [prefix].pdb              Protein structure with mapped conservation scores
# • [prefix].txt              Residue-by-residue conservation table
# • [prefix].pml              PyMOL visualization script
# • [prefix].cxc              ChimeraX visualization script
# • [prefix]_alignment.html   Interactive conservation-colored alignment viewer
#
#
# APPLICATIONS
# ------------
# • Identification of conserved functional regions
# • Structural interpretation of sequence conservation
# • Comparative analysis of protein families
# • Preparation of publication-quality figures
# • Interactive exploration of sequence–structure relationships
#
# ==============================================================================
import math
import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from collections import Counter
import warnings

# Direct imports (Removed the try/except block for better Linux compatibility)
from Bio import AlignIO, PDB
from Bio.Align import substitution_matrices
from Bio.SeqUtils import seq1
from Bio import BiopythonDeprecationWarning
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

# Suppress Biopython warnings
warnings.simplefilter('ignore', BiopythonDeprecationWarning)

# ============================================================
# Constants for the new scoring metrics
# ============================================================

# Physicochemical grouping used by the "property" metric. Six broad classes;
# feel free to refine (e.g. splitting "special" or separating His out of
# "positive") if your protein family needs finer granularity.
PROPERTY_GROUPS = {
    'A': 'aliphatic', 'V': 'aliphatic', 'L': 'aliphatic', 'I': 'aliphatic', 'M': 'aliphatic',
    'F': 'aromatic', 'W': 'aromatic', 'Y': 'aromatic',
    'S': 'polar', 'T': 'polar', 'N': 'polar', 'Q': 'polar', 'C': 'polar',
    'K': 'positive', 'R': 'positive', 'H': 'positive',
    'D': 'negative', 'E': 'negative',
    'G': 'special', 'P': 'special',
}

# Standard amino acid background frequencies (Robinson & Robinson, 1991),
# used as the reference distribution for the Jensen-Shannon Divergence metric.
BACKGROUND_FREQ = {
    'A': 0.078, 'R': 0.051, 'N': 0.041, 'D': 0.052, 'C': 0.024,
    'Q': 0.034, 'E': 0.059, 'G': 0.083, 'H': 0.025, 'I': 0.062,
    'L': 0.092, 'K': 0.056, 'M': 0.024, 'F': 0.044, 'P': 0.043,
    'S': 0.059, 'T': 0.055, 'W': 0.014, 'Y': 0.034, 'V': 0.072,
}
_bg_total = sum(BACKGROUND_FREQ.values())
BACKGROUND_FREQ = {aa: f / _bg_total for aa, f in BACKGROUND_FREQ.items()}

METRIC_LABELS = {
    "blosum": "BLOSUM62",
    "entropy": "Shannon Entropy",
    "property": "Property Conservation",
    "jsd": "Jensen-Shannon Divergence",
    "rate4site": "Rate4Site-like (proxy)",
}

# ============================================================
# Scoring helper functions
# ============================================================

def shannon_entropy_score(valid_residues):
    """Classic Shannon-entropy conservation score, 0-100 (100 = fully conserved)."""
    if not valid_residues:
        return 0
    counts = Counter(valid_residues)
    total = len(valid_residues)
    H = sum(-(c / total) * math.log2(c / total) for c in counts.values())
    return (1 - H / math.log2(20)) * 100


def blosum_column_score(valid_residues, matrix):
    """Average pairwise BLOSUM62 self/cross similarity for a column, 0-100."""
    valid = [aa for aa in valid_residues if aa in matrix.alphabet]
    if not valid:
        return 0
    raw_score = sum(matrix[a, a] for a in valid)
    for i, aa1 in enumerate(valid):
        for aa2 in valid[i + 1:]:
            raw_score += 2 * matrix[aa1, aa2]
    max_score = len(valid) ** 2 * max(matrix[a, a] for a in valid)
    return (raw_score / max_score) * 100 if max_score > 0 else 0


def property_conservation_score(valid_residues):
    """
    Groups residues into physicochemical classes (see PROPERTY_GROUPS) and
    scores a column by the fraction belonging to the majority class.
    A column that is 100% hydrophobic (even if the exact residue varies,
    e.g. L/I/V) scores near 100.
    """
    groups = [PROPERTY_GROUPS[aa] for aa in valid_residues if aa in PROPERTY_GROUPS]
    if not groups:
        return 0
    counts = Counter(groups)
    majority_fraction = counts.most_common(1)[0][1] / len(groups)
    return majority_fraction * 100


def jsd_conservation_score(valid_residues, background_freq=BACKGROUND_FREQ):
    """
    Jensen-Shannon Divergence between a column's amino-acid distribution and
    a fixed background distribution (Capra & Singh, 2007 style, without the
    sliding-window smoothing they apply). JSD is bounded in [0, 1] bit for
    two distributions, so the result is simply scaled by 100.
    """
    valid = [aa for aa in valid_residues if aa in background_freq]
    if not valid:
        return 0
    counts = Counter(valid)
    total = len(valid)
    col_freq = {aa: counts.get(aa, 0) / total for aa in background_freq}

    mixture = {aa: 0.5 * (col_freq[aa] + background_freq[aa]) for aa in background_freq}

    def kl_divergence(p, q):
        return sum(p[aa] * math.log2(p[aa] / q[aa]) for aa in p if p[aa] > 0 and q[aa] > 0)

    jsd = 0.5 * kl_divergence(col_freq, mixture) + 0.5 * kl_divergence(background_freq, mixture)
    return max(0.0, min(jsd, 1.0)) * 100


def compute_henikoff_weights(alignment):
    """
    Henikoff & Henikoff (1994) position-based sequence weighting.
    Down-weights sequences that are over-represented in the alignment
    (near-duplicates contribute less), which is the key ingredient the
    "rate4site-like" proxy borrows from real evolutionary-rate estimators.
    Returns one weight per sequence, normalized to sum to len(alignment).
    """
    n_seq = len(alignment)
    aln_len = alignment.get_alignment_length()
    weights = [0.0] * n_seq

    for col in range(aln_len):
        column = alignment[:, col]
        counts = Counter(column)
        distinct = len(counts)
        if distinct == 0:
            continue
        for i, aa in enumerate(column):
            weights[i] += 1.0 / (distinct * counts[aa])

    total = sum(weights)
    if total == 0:
        return [1.0] * n_seq
    return [w * n_seq / total for w in weights]


def rate4site_like_score(column, weights):
    """
    SIMPLIFIED PROXY for site-specific evolutionary rate, NOT the actual
    Rate4Site algorithm (which requires ML rate estimation on a phylogenetic
    tree). Here we compute a sequence-weighted Shannon entropy: residues from
    over-represented (near-identical) sequences count less, giving a rate
    estimate that is less biased by uneven taxon sampling than raw entropy.
    Returned as a conservation score, 0-100 (100 = most conserved / slowest).
    """
    weighted_counts = {}
    total_weight = 0.0
    for aa, w in zip(column, weights):
        if aa == "-":
            continue
        weighted_counts[aa] = weighted_counts.get(aa, 0.0) + w
        total_weight += w

    if total_weight == 0:
        return 0

    H = 0.0
    for c in weighted_counts.values():
        p = c / total_weight
        if p > 0:
            H += -p * math.log2(p)
    return (1 - H / math.log2(20)) * 100


def compute_scores(alignment, metric="blosum", gap_penalty=False):
    """
    Dispatches to the appropriate per-column scoring function.
    metric one of: "blosum", "entropy", "property", "jsd", "rate4site".
    All metrics are normalized to a common 0-100 scale so downstream
    B-factor mapping / color spectra work unchanged regardless of choice.
    """
    scores = []
    aln_len = alignment.get_alignment_length()

    matrix = substitution_matrices.load("BLOSUM62") if metric == "blosum" else None
    # Only needed for the rate4site-like proxy; cheap to skip otherwise.
    henikoff_weights = compute_henikoff_weights(alignment) if metric == "rate4site" else None

    for col in range(aln_len):
        column = alignment[:, col]
        valid = [aa for aa in column if aa != "-"]

        if metric == "entropy":
            score = shannon_entropy_score(valid)
        elif metric == "property":
            score = property_conservation_score(valid)
        elif metric == "jsd":
            score = jsd_conservation_score(valid)
        elif metric == "rate4site":
            score = rate4site_like_score(column, henikoff_weights)
        else:  # "blosum" (default/fallback)
            score = blosum_column_score(valid, matrix)

        if gap_penalty:
            score *= (1 - column.count("-") / len(column))
        scores.append(score)
    return scores

def get_cmap(name):
    if name.lower() == "molnympheas":
        colors = ['#1d3557', '#457b9d', '#a8d5ba', '#f1faee', '#fbb4b9', '#b56576']
        return mcolors.LinearSegmentedColormap.from_list("molnympheas", colors)
    try: return plt.get_cmap(name)
    except: return plt.get_cmap('viridis')

def generate_pml(out_pdb, out_pml, palette, step_value):
    cmap = get_cmap(palette)
    # Calculate the number of intervals based on the user's step value
    num_steps = int(100 / step_value) + 1
    intervals = np.linspace(0, 100, num_steps)
    
    color_defs = ""
    color_names = []
    for i, val in enumerate(intervals):
        rgba = cmap(val / 100)
        cname = f"c_{palette}_{i}"
        color_defs += f"set_color {cname}, [{rgba[0]:.3f}, {rgba[1]:.3f}, {rgba[2]:.3f}]\n"
        color_names.append(cname)
        
    # Use commas to separate values in the lists for PyMOL
    intervals_list = ", ".join(f"{v:.0f}" for v in intervals)
    colors_list = ", ".join(color_names)

    with open(out_pml, "w") as f:
        f.write(f"""load {os.path.basename(out_pdb)}, conservation_struct
bg_color white
hide all
show surface, conservation_struct

{color_defs}
# spectrum command uses spaces
spectrum b, {" ".join(color_names)}, conservation_struct, minimum=0, maximum=100

# ramp_new command MUST use commas
ramp_new conservation_ramp, conservation_struct, [{intervals_list}], [{colors_list}]
zoom conservation_struct
""")
def generate_alignment_html(alignment, scores, output_html, palette, metric,
                            block_size=60):
    """
    Generate an ESPript-like HTML alignment viewer.
    """

    cmap = get_cmap(palette)

    def score_to_hex(score):
        r, g, b, _ = cmap(score / 100.0)
        return "#{:02x}{:02x}{:02x}".format(
            int(r * 255),
            int(g * 255),
            int(b * 255)
        )

    html = []

    html.append(f"""
<!DOCTYPE html>
<html>

<head>

<meta charset="utf-8">

<title>Alignment viewer</title>

<style>

body{{
    font-family:Arial,sans-serif;
    margin:30px;
}}

h1{{
    margin-bottom:5px;
}}

.info{{
    margin-bottom:20px;
}}

.legend{{
    width:700px;
    height:20px;
    border:1px solid #999;
    background:linear-gradient(to right,
""")

    grad=[]

    for i in range(101):
        grad.append(f"{score_to_hex(i)} {i}%")

    html.append(",".join(grad))

    html.append("""
);
}

.scale{
width:700px;
display:flex;
justify-content:space-between;
margin-bottom:25px;
font-size:12px;
}

table{
border-collapse:collapse;
font-family:Consolas,monospace;
font-size:15px;
margin-bottom:18px;
}

td{
padding:0;
}

.name{
font-weight:bold;
padding-right:20px;
white-space:nowrap;
min-width:180px;
}

.num{
height:18px;
font-size:11px;
text-align:center;
color:#555;
}

.res{
width:16px;
height:20px;
text-align:center;
}

.ref{
font-weight:bold;
border-radius:2px;
}

</style>

</head>

<body>

<h1>Alignment Conservation</h1>

<div class="info">
<b>Metric:</b> """+metric+"""<br>
<b>Palette:</b> """+palette+"""
</div>

<div class="legend"></div>

<div class="scale">
<span>0</span>
<span>50</span>
<span>100</span>
</div>

""")
    # ============================================================
    # Alignment blocks
    # ============================================================

    aln_len = alignment.get_alignment_length()
    ref_number = 0

    for start in range(0, aln_len, block_size):

        end = min(start + block_size, aln_len)

        html.append("<table>")

        # --------------------------------------------------------
        # Numbering row
        # --------------------------------------------------------

        html.append("<tr>")
        html.append("<td class='name'></td>")

        tmp_number = ref_number

        for i in range(start, end):

            aa = alignment[0].seq[i]

            if aa != "-":
                tmp_number += 1

            if aa != "-" and tmp_number % 10 == 0:
                html.append(
                    f"<td class='num' colspan='1'>{tmp_number}</td>"
                )
            else:
                html.append("<td class='num'></td>")

        html.append("</tr>")

        # --------------------------------------------------------
        # Alignment
        # --------------------------------------------------------

        for seq_index, record in enumerate(alignment):

            html.append("<tr>")

            html.append(
                f"<td class='name'>{record.id}</td>"
            )

            current_ref = ref_number

            for i in range(start, end):

                aa = record.seq[i]

                score = scores[i]

                color = score_to_hex(score)

                if alignment[0].seq[i] != "-":
                    current_ref += 1

                if aa == "-":

                    html.append(
                        "<td class='res'>-</td>"
                    )

                else:

                    if seq_index == 0:

                        html.append(
                            f"<td class='res ref' "
                            f"style='background:{color};' "
                            f"title='Residue {current_ref}\n"
                            f"Score {score:.2f}'>"
                            f"{aa}</td>"
                        )

                    else:

                        html.append(
                            f"<td class='res' "
                            f"style='background:{color};' "
                            f"title='Reference position {current_ref}\n"
                            f"Score {score:.2f}'>"
                            f"{aa}</td>"
                        )

            html.append("</tr>")

        html.append("</table>")

        # Update reference numbering for next block

        for i in range(start, end):
            if alignment[0].seq[i] != "-":
                ref_number += 1

    # ============================================================
    # Footer
    # ============================================================

    html.append("""

</body>

</html>

""")

    with open(output_html, "w", encoding="utf-8") as f:
        f.write("".join(html))        



def generate_cxc(out_pdb, out_cxc, palette, step_value):
    """Generate ChimeraX script using exact Matplotlib colors."""
    cmap = get_cmap(palette)
    num_steps = int(100 / step_value) + 1
    intervals = np.linspace(0, 100, num_steps)

    # Convert matplotlib RGB to hex colors for ChimeraX
    hex_colors = []
    for val in intervals:
        rgba = cmap(val / 100.0)
        r = int(rgba[0] * 255)
        g = int(rgba[1] * 255)
        b = int(rgba[2] * 255)
        hex_colors.append(f"#{r:02x}{g:02x}{b:02x}")

    palette_string = ":".join(hex_colors)

    with open(out_cxc, "w") as f:
        f.write(f"""# ChimeraX visualization script
open {os.path.basename(out_pdb)}
hide all
show surface
color byattribute bfactor palette {palette_string} range 0,100
key
""")

def process_structure(alignment_file, structure_file, output_dir, metric, gap_penalty, palette, step_value):
    # Create directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    base_name = os.path.basename(output_dir)
    
    alignment = AlignIO.read(alignment_file, "fasta")
    ref_seq = str(alignment[0].seq)
    scores = compute_scores(alignment, metric, gap_penalty)

    parser = PDB.MMCIFParser(QUIET=True) if structure_file.lower().endswith(".cif") else PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("struct", structure_file)
    residues = [r for m in structure for c in m for r in c if PDB.is_aa(r, standard=True)]
    
    mapped_scores = [scores[i] for i, aa in enumerate(ref_seq) if aa != "-"]
    for residue, score in zip(residues, mapped_scores):
        for atom in residue: atom.set_bfactor(score)

    pdb_path = os.path.join(output_dir, f"{base_name}.pdb")
    txt_path = os.path.join(output_dir, f"{base_name}.txt")
    pml_path = os.path.join(output_dir, f"{base_name}.pml")
    cxc_path = os.path.join(output_dir, f"{base_name}.cxc")
    html_path = os.path.join(output_dir, f"{base_name}_alignment.html")

    io = PDB.PDBIO(); io.set_structure(structure); io.save(pdb_path)
    with open(txt_path, "w") as f:
        f.write(f"# Metric: {METRIC_LABELS.get(metric, metric)}\n")
        f.write(f"# Gap penalty applied: {gap_penalty}\n")
        f.write("Chain\tResID\tResName\tScore\n")
        for residue, score in zip(residues, mapped_scores):
            f.write(f"{residue.parent.id}\t{residue.id[1]}\t{residue.resname}\t{score:.3f}\n")

    generate_pml(pdb_path, pml_path, palette, step_value)
    generate_cxc(pdb_path, cxc_path, palette, step_value)
    generate_alignment_html(
    alignment,
    scores,
    html_path,
    palette,
    METRIC_LABELS.get(metric, metric)
)

    return {"dir": output_dir}

def validate_sequence(residues, ref_seq):
    structure_seq = "".join(seq1(r.get_resname(), custom_map={}) for r in residues)
    msa_seq = "".join(aa for aa in ref_seq if aa != "-")

    if structure_seq != msa_seq:
        raise ValueError("Sequence mismatch between structure and alignment.")

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("CS-to-PDB GUI v7")
        self.geometry("820x760") # Taller to accommodate the extended scoring info text
        self.columnconfigure(1, weight=1)
        self.create_widgets()

    def create_widgets(self):
        # File Inputs
        ttk.Label(self, text="Alignment (FASTA):").grid(row=0, column=0, padx=5, pady=5, sticky='e')
        self.aln = ttk.Entry(self); self.aln.grid(row=0, column=1, sticky="ew")
        ttk.Button(self, text="Browse", command=lambda: self.browse(self.aln)).grid(row=0, column=2, padx=5)

        ttk.Label(self, text="Structure (PDB/CIF):").grid(row=1, column=0, padx=5, pady=5, sticky='e')
        self.pdb = ttk.Entry(self); self.pdb.grid(row=1, column=1, sticky="ew")
        ttk.Button(self, text="Browse", command=lambda: self.browse(self.pdb)).grid(row=1, column=2, padx=5)

        ttk.Label(self, text="Output Prefix:").grid(row=2, column=0, padx=5, pady=5, sticky='e')
        self.prefix = ttk.Entry(self); self.prefix.insert(0, "CS-to-PDB-result"); self.prefix.grid(row=2, column=1, sticky="ew")

        # Config Frame
        frame = ttk.Frame(self); frame.grid(row=3, column=0, columnspan=3, pady=10)
        self.metric = ttk.Combobox(
            frame,
            values=["blosum", "entropy", "property", "jsd", "rate4site"],
            state="readonly", width=12
        )
        self.metric.current(0); self.metric.pack(side=tk.LEFT, padx=5)
        
        # Dynamically fetch all matplotlib colormaps and add custom ones to the front
        all_palettes = ["molnympheas"] + list(plt.colormaps())
        self.palette = ttk.Combobox(frame, values=all_palettes, state="readonly", width=15); self.palette.current(0); self.palette.pack(side=tk.LEFT, padx=5)
        self.palette.bind("<<ComboboxSelected>>", self.update_preview)
        
        # Gradient Step Controller
        ttk.Label(frame, text="Step %:").pack(side=tk.LEFT, padx=(5, 2))
        self.step_var = tk.IntVar(value=10)
        self.step_cb = ttk.Combobox(frame, textvariable=self.step_var, values=[5, 10, 20, 25, 50], state="readonly", width=4)
        self.step_cb.pack(side=tk.LEFT, padx=5)

        self.gap = tk.BooleanVar(); ttk.Checkbutton(frame, text="Penalize gaps", variable=self.gap).pack(side=tk.LEFT, padx=5)

        # Info Box
        info = ttk.LabelFrame(self, text="Scoring Information", padding=10)
        info.grid(row=4, column=0, columnspan=3, padx=10, pady=5, sticky="ew")
        info_text = (
            "• BLOSUM: Uses evolutionary substitution rates (BLOSUM62). High score = frequently conserved amino acids.\n\n"
            "• Entropy: Uses Shannon entropy to measure sequence diversity. High score = low diversity (highly conserved).\n\n"
            "• Property Conservation: Groups amino acids by physicochemical class (aliphatic, aromatic, polar, "
            "positive, negative, special). High score = column dominated by one chemical class, even if the exact residue varies.\n\n"
            "• JSD: Jensen-Shannon Divergence between a column's composition and a standard background distribution. "
            "High score = strongly biased/conserved column.\n\n"
            "• Rate4Site-like: SIMPLIFIED PROXY only \u2014 sequence-weighted entropy (Henikoff & Henikoff weights), "
            "not the true ML/phylogenetic Rate4Site algorithm.\n\n"
            "• Penalize Gaps: Reduces the final score of a column proportionally to the number of gaps ('-') it contains."
        )
        ttk.Label(info, text=info_text, wraplength=650, justify="left").pack(fill='x')

        self.preview = tk.Canvas(self, width=400, height=20, bd=1, relief="sunken")
        self.preview.grid(row=5, column=0, columnspan=3, pady=5)
        self.update_preview()

        ttk.Button(self, text="Run Mapping", command=self.run).grid(row=6, column=0, columnspan=3, pady=10)

    def browse(self, entry):
        file = filedialog.askopenfilename()
        if file: entry.delete(0, tk.END); entry.insert(0, file)

    def update_preview(self, event=None):
        self.preview.delete("all")
        cmap = get_cmap(self.palette.get())
        for x in range(400):
            rgba = cmap(x / 400)
            color = "#%02x%02x%02x" % (int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))
            self.preview.create_line(x, 0, x, 20, fill=color)

    def run(self):
        try:
            # Pass the step value into the processor
            stats = process_structure(
                self.aln.get(), 
                self.pdb.get(), 
                self.prefix.get(), 
                self.metric.get(), 
                self.gap.get(), 
                self.palette.get(),
                self.step_var.get()
            )
            messagebox.showinfo("Success", "Files generated successfully.")
        except Exception as e: 
            messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    app = App()
    app.mainloop()
