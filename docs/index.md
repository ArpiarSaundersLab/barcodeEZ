# barcodeEZ

**barcodeEZ** is a Python package for designing complex, combinatorial DNA barcode libraries for molecular cloning.

The package allows design and creation of a site-aware `Barcodes` object. It is built around sites, expandable positions, restriction-enzyme boundaries, optimized overhangs, and optional fixed sequences. BarcodeEZ draws orthogonal barcode sequence from a prebuilt corpus, assembles forward and reverse oligos, screens them for unwanted sequence content, and exports a ready-to-order oligo pool.

## What it does

- **Design** a library structure from restriction-enzyme boundaries — use the built-in default panel or supply your own enzymes.
- **Add positions** within each site for combinatorial, multi-position barcoding (up to 8 positions per site), with automatically assigned optimized internal overhangs.
- **Generate barcodes** of any length from a corpus of ~30,000 orthogonal 60-mers shipped with the package.
- **Attach fixed sequences** to either end of a site.
- **Validate** against restriction sites and undesired motifs, automatically swapping out any contaminated barcode.
- **Inspect** the whole design as a pandas DataFrame.
- **Export** a single-stranded oligo order form as CSV.

## A minimal example

```python
from barcodeEZ import Barcodes

b = Barcodes(n_sites=2)             # two barcoded sites, default enzymes
b.generate_barcodes(bc_len=20, n_barcodes=96)
b.validate()                        # screen + replace contaminated barcodes
b.write_order_form('library.csv')   # export the oligo pool
```

## Quick links

- [Installation](getting-started/installation.md)
- [Quickstart](getting-started/quickstart.md)
- [Key Concepts](usage/concepts.md)
- [API Reference](usage/api.md)
