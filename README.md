# barcodeEZ

<img src="docs/assets/barcodeEZ.png" alt="barcodeEZ" width="400"/>

**barcodeEZ** is a Python package for designing complex, combinatorial DNA barcode libraries for molecular cloning.

The package allows design and creation of a site-aware `Barcodes` object. It is built around sites, expandable positions, restriction-enzyme boundaries, optimized overhangs, and optional fixed sequences. barcodeEZ draws orthogonal barcode sequence from a prebuilt corpus, assembles forward and reverse oligos, screens them for unwanted sequence content, and exports a ready-to-order oligo pool.

## What it does

- **Design** a library structure from restriction-enzyme boundaries — use the built-in default panel or supply your own enzymes.
- **Add positions** within each site for combinatorial, multi-position barcoding (up to 8 positions per site), with automatically assigned optimized internal overhangs.
- **Generate barcodes** of any length from a corpus of ~20,000 orthogonal 60-mers shipped with the package.
- **Attach fixed sequences** to either end of a site.
- **Validate** against restriction sites and undesired motifs, automatically swapping out any contaminated barcode.
- **Inspect** the whole design as a pandas DataFrame.
- **Export** a single-stranded oligo order form as CSV.

## Installation

```bash
pip install barcodeEZ
```

The barcode corpus ships with the package; `biopython` and `pandas` are installed automatically.

## A minimal example

```python
from barcodeEZ import Barcodes

b = Barcodes(n_sites=2)             # two barcoded sites, default enzymes
b.generate_barcodes(bc_len=20, n_barcodes=96)
b.validate()                        # screen + replace contaminated barcodes
b.write_order_form('library.csv')   # export the oligo pool
```

## API

All functionality is accessed through the `Barcodes(n_sites, custom_enzymes=None)` object. Methods return `self` and can be chained.

| Method | Purpose |
|--------|---------|
| `add_positions(*, n_per_site)` | Add combinatorial positions per site (call before generating) |
| `generate_barcodes(bc_len, n_barcodes)` | Draw barcodes and assemble oligos |
| `add_fixed_sequence(seq, site, side)` | Attach a fixed flanking sequence to a site |
| `validate(ignore_defaults=False, motifs=None)` | Screen and replace barcodes containing unwanted motifs |
| `view()` | Return the full library as a pandas DataFrame |
| `write_order_form(file)` | Export the oligo pool as CSV |
| `print_structure()` | Print the site/enzyme layout |

See the [full documentation](https://goodez.github.io/barcodeEZ/) for parameter details and examples.

## Requirements

- Python ≥ 3.10
- biopython ≥ 1.81, pandas (installed automatically with pip)
