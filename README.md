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
  
<img width="932" height="900" alt="image" src="https://github.com/user-attachments/assets/793883a3-88c7-414e-a560-fb1e184c8bba" />

- Uses identical color palettes for:
  - Protein structures
  - HTML alignment viewer

<img width="1264" height="710" alt="image" src="https://github.com/user-attachments/assets/91a47061-cbba-43ac-9bb8-6d3229864095" />


  - PyMOL
    
<img width="1232" height="982" alt="image" src="https://github.com/user-attachments/assets/b907a9b1-2b5b-40fe-8a1f-5faab9d3898e" />


  - ChimeraX

<img width="1139" height="982" alt="image" src="https://github.com/user-attachments/assets/16b1c41d-7325-420d-b475-69e3f7493c1a" />


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

## License
Do whatever you want with the scripts and keep in mind that it has been code by a stupid biochemist (for that matter me).
