# API Reference

The entire public API is the `Barcodes` class.

```python
from barcodeEZ import Barcodes
```

Most methods **return the `Barcodes` object itself**, so calls can be chained:

```python
b = (Barcodes(n_sites=2)
        .add_positions(n_per_site=4)
        .generate_barcodes(bc_len=20, n_barcodes=96)
        .validate())
```

::: barcodeEZ.Barcodes
