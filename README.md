# barcodeEZ

**barcodeEZ** is a Python package for designing complex, combinatorial DNA barcode libraries for molecular cloning.

You create a site-aware `Barcodes` object, describe the structure of your library (sites, positions, restriction-enzyme boundaries, fixed sequences), and barcodeEZ draws orthogonal barcode sequence from a prebuilt corpus, assembles forward and reverse oligos, screens them for unwanted sequence content, and exports a ready-to-order oligo pool.

## Features

- Design a library from restriction-enzyme boundaries — built-in default panel or your own enzymes
- Add up to 8 combinatorial positions per site, with automatically assigned internal overhangs
- Generate barcodes of any length from a corpus of ~20,000 orthogonal 60-mers bundled with the package
- Attach fixed sequences (adapters, spacers, anchors) to either side of a site
- Validate against restriction sites and arbitrary motifs, automatically replacing contaminated barcodes
- Inspect the full design as a pandas DataFrame
- Export a single-stranded oligo order form as CSV

## Installation

barcodeEZ requires **Python ≥ 3.9**.

```bash
pip install barcodeEZ
```

Or install the latest development version from GitHub:

```bash
pip install git+https://github.com/goodez/barcodeEZ.git
```

The barcode corpus ships with the package; `biopython` and `pandas` are installed automatically.

## Quickstart

```python
from barcodeEZ import Barcodes

b = Barcodes(n_sites=2)              # two barcoded sites, default enzymes
b.generate_barcodes(bc_len=20, n_barcodes=96)
b.validate()                        # screen + replace contaminated barcodes
b.write_order_form('library.csv')   # export the oligo pool
```

For a combinatorial library with positions and adapters:

```python
b = Barcodes(n_sites=1)
b.add_positions(n_per_site=3)              # positions A, B, C per site
b.generate_barcodes(bc_len=18, n_barcodes=96)
b.add_fixed_sequence('AAGCTT', site=1, side='left')
b.validate()
print(b.view())
```

## API overview

| Method | Purpose |
|--------|---------|
| `Barcodes(n_sites, custom_enzymes=None)` | Create a library structure |
| `add_positions(*, n_per_site)` | Add combinatorial positions per site (call before generating) |
| `generate_barcodes(bc_len, n_barcodes)` | Draw barcodes and assemble oligos |
| `add_fixed_sequence(seq, site, side)` | Attach a fixed flanking sequence to a site |
| `validate(ignore_defaults=False, motifs=None)` | Screen and replace barcodes containing unwanted motifs |
| `view()` | Return the full library as a pandas DataFrame |
| `write_order_form(file)` | Export the oligo pool as CSV |
| `print_structure()` | Print the site/enzyme layout |

See the [full documentation](https://goodez.github.io/barcodeEZ/) for parameter details and examples.

## Requirements

- Python ≥ 3.9
- biopython ≥ 1.81, pandas (installed automatically with pip)
