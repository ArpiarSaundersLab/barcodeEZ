# Quickstart

This walkthrough builds a small library end to end. For the full set of parameters, see the [API Reference](../usage/api.md).

## 1. Create a library

A `Barcodes` object is defined by the number of barcoded **sites**. Each site is bounded by a pair of restriction enzymes drawn from a default panel.

```python
from barcodeEZ import Barcodes

b = Barcodes(n_sites=2)
b.print_structure()
```

```text
----- EcoRI - SITE1 - BamHI - SITE2 - NheI -----
```

Each site sits between two enzymes; adjacent sites share a boundary enzyme.

## 2. Generate barcodes

Draw barcodes from the corpus. `bc_len` sets the length of each barcode; `n_barcodes` sets how many per site.

```python
b.generate_barcodes(bc_len=20, n_barcodes=2)
```

barcodeEZ assembles the forward and reverse oligos automatically, adding the correct enzyme overhang to each end.

## 3. Validate against unwanted motifs

`validate()` screens every barcode against a default panel of restriction sites and replaces any contaminated barcode with a fresh draw from the corpus.

```python
b.validate()
```

```text
Validation complete: no unwanted motifs found.
```

## 4. Inspect the design

`view()` returns a pandas DataFrame of the whole library.

```python
b.view()
```

```text
 site position              barcode              forward_oligo              reverse_oligo
    1        A CGTTCACGGTAACGCTACGT AATTCCGTTCACGGTAACGCTACGTG GATCCACGTAGCGTTACCGTGAACGG
    1        A ATTGTGCTCTCGCGCGGACC AATTCATTGTGCTCTCGCGCGGACCG GATCCGGTCCGCGCGAGAGCACAATG
    2        A AGCTTCAGGAGTCTCCATCG GATCCAGCTTCAGGAGTCTCCATCGG CTAGCCGATGGAGACTCCTGAAGCTG
    2        A GATCTGTTCGGAACTAATCC GATCCGATCTGTTCGGAACTAATCCG CTAGCGGATTAGTTCCGAACAGATCG
```

## 5. Export the order form

Write the single-stranded oligo pool to CSV. Each barcode yields two rows — a forward (`_f`) and reverse (`_r`) oligo — grouped by site.

```python
b.write_order_form('library.csv')
```

```text
opool_name,oligo_sequence
site1_f,AATTCCGTTCACGGTAACGCTACGTG
site1_r,GATCCACGTAGCGTTACCGTGAACGG
...
```

Pass `metadata=True` to also record which site, position, and barcode each oligo
came from — useful for tracking the pool back to the design after synthesis.

```python
b.write_order_form('library.csv', metadata=True)
```

```text
opool_name,oligo_sequence,site,position,barcode
site1_f,AATTCCGTTCACGGTAACGCTACGTG,1,A,CGTTCACGGTAACGCTACGT
site1_r,GATCCACGTAGCGTTACCGTGAACGG,1,A,CGTTCACGGTAACGCTACGT
...
```

## Going further: positions and fixed sequences

One can customize libraries further by adding internal positions within sites. Positions are connected via optimized 4 bp overhangs. This allows the total number of combined barcodes to expand beyond n = 6. Positions should be added before barcode generation, as the method will wipe clean any barcodes in the object. 

Fixed sequence can also be appended to either terminal end of a desired site. Users must provide the sequence, site number, and desired side to append to (left or right).

```python
b = Barcodes(n_sites=1)
b.add_positions(n_per_site=3)        # positions A, B, C per site
b.generate_barcodes(bc_len=18, n_barcodes=1)
b.add_fixed_sequence('AAGCTT', site=1, side='left')
b.validate()
print(b.view().to_string(index=False))
```

See [Key Concepts](../usage/concepts.md) for how positions, overhangs, and fixed sequences fit together.
