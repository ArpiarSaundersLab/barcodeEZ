# Installation

barcodeEZ requires **Python ≥ 3.9**.

## Install with pip

```bash
pip install barcodeEZ
```

## Dependencies

These are installed automatically:

| Package | Purpose |
|---------|---------|
| [biopython](https://biopython.org/) (≥ 1.81) | Restriction-enzyme recognition sites |
| [pandas](https://pandas.pydata.org/) | Tabular views and CSV export |

The barcode corpus (~30,000 orthogonal 60-mers) ships with the package — no separate download is needed.

## Verify the installation

```python
from barcodeEZ import Barcodes

b = Barcodes(n_sites=1)
b.generate_barcodes(bc_len=20, n_barcodes=3)
print(b.view())
```
