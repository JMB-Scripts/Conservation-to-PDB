# CS-to-PDB

**CS-to-PDB** is a standalone Python application that calculates residue conservation from a Multiple Sequence Alignment (MSA) and maps evolutionary information onto three-dimensional protein structures by storing conservation scores in the PDB B-factor field. The resulting structures can be visualized directly in **PyMOL** and **ChimeraX**, while an interactive **HTML alignment viewer** provides a convenient way to compare sequence conservation with structural data.

---

## Features

- Multiple conservation metrics:
  - BLOSUM62 similarity
  - Shannon Entropy
  - Property Conservation
  - Jensen-Shannon Divergence (JSD)
  - Rate4Site-like weighted conservation
- Optional gap penalty for poorly aligned regions.
- Supports **PDB** and **mmCIF** structures.
- Maps conservation scores directly into the **B-factor** field.
- Generates visualization scripts for **PyMOL** and **ChimeraX**.
- Generates an interactive **HTML alignment viewer**.
- Supports all **Matplotlib** colormaps, including custom palettes (e.g. MolNympheas).
- Uses identical color palettes for:
  - Protein structures
  - HTML alignment viewer
  - PyMOL
  - ChimeraX
- Simple graphical user interface.
- Cross-platform (Linux, macOS and Windows).

---

## Output files

```
Project/
│
├── protein.pdb
├── conservation.txt
├── visualization.pml
├── visualization.cxc
└── alignment.html
```

| File | Description |
|------|-------------|
| **.pdb** | Protein structure with conservation scores stored in the B-factor field |
| **.txt** | Residue-by-residue conservation table |
| **.pml** | PyMOL visualization script |
| **.cxc** | ChimeraX visualization script |
| **.html** | Interactive conservation-colored alignment viewer |

---

## Typical workflow

```
Multiple Sequence Alignment
             │
             ▼
   Conservation Calculation
             │
             ▼
      Score Normalization
             │
      ┌──────┴──────────┐
      ▼                 ▼
 Protein Structure   Alignment
      │                 │
      ▼                 ▼
  PDB (B-factor)    HTML Viewer
      │
      ├────────► PyMOL
      └────────► ChimeraX
```

---

## Requirements

- Python 3.9+
- Biopython
- NumPy
- Matplotlib
- Tkinter (included with most Python installations)

Install the required packages with:

```bash
pip install biopython numpy matplotlib
```

---

## Applications

- Visualizing evolutionary conservation on protein structures.
- Identifying conserved functional regions.
- Comparing sequence conservation with structural features.
- Preparing publication-quality figures.
- Teaching structural biology and sequence analysis.

---


## License
Do whatever you want with the scripts and keep in mind that it has been code by a stupid biochemist (for that matter me).
