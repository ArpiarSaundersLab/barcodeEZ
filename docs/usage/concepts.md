# Key Concepts

This page explains the building blocks of a barcodeEZ library and the order in which they fit together.

## The design workflow

A library is built up in a fixed sequence of steps. Later steps depend on earlier ones.

```text
Barcodes(n_sites)          create the library structure
      │
add_positions()            (optional) combinatorial positions per site
      │
generate_barcodes()        draw barcodes + assemble oligos
      │
add_fixed_sequence()       (optional) attach fixed sequence per site
      │
validate()                 screen + replace contaminated barcodes
      │
view() / write_order_form()  inspect and export
```

!!! note "Order matters"
    `add_positions()` resets any barcodes, so call it **before** `generate_barcodes()`. `add_fixed_sequence()` and `validate()` come **after** generation. Re-running `generate_barcodes()` or `add_fixed_sequence()` marks the library as unvalidated again.

## Sites

A **site** is one barcoded locus in your construct. Every site is bounded on its left and right by a restriction enzyme. With *N* sites you have *N + 1* enzyme boundaries, and adjacent sites share the enzyme between them:

```text
EcoRI - SITE1 - BamHI - SITE2 - NheI
```

By default, sites use this enzyme panel in order:

```text
EcoRI, BamHI, NheI, XhoI, PluTI, AgeI, MluI
```

so the default configuration supports up to **6 sites**. Supply your own panel with `custom_enzymes` (you must provide exactly `n_sites + 1` enzyme names).

## Positions

A **position** is a barcode slot *within* a site. By default each site has a single position, labelled `A`. Calling `add_positions(n_per_site=k)` gives every site positions `A, B, …` (up to 8, labelled `A`–`H`).

Positions enable expanded **combinatorial** sites: instead of one barcode per site, you assemble several adjacent barcodes that ligate together in order.

## Overhangs

Adjacent positions within a site are joined by 4 bp **overhangs** (sticky ends), so the fragments assemble in a defined order. barcodeEZ assigns these automatically from a built-in set when you call `add_positions()`. Neighbouring positions share an overhang — the right overhang of position `A` is the left overhang of position `B`, and so on.

At the **outer** edges of a site (the left edge of the first position, the right edge of the last), there is no internal overhang; instead the boundary restriction enzyme's recognition sequence is used when building the oligo.

## Barcodes and the corpus

Raw barcode sequence is drawn from a **corpus** of ~20,000 orthogonal 60-mers bundled with the package. Each call to `generate_barcodes()` draws without replacement, so barcodes are unique across the whole library. Requesting a `bc_len` longer than 60 bp concatenates multiple corpus segments.

## Fixed sequences

A **fixed sequence** is constant flanking sequence — an adapter, spacer, or anchor — attached to the left or right side of a site with `add_fixed_sequence()`. It applies to every position in that site and sits immediately adjacent to the barcode, inside the enzyme/overhang ends.

The assembled barcode is therefore:

```text
[fixed_left] + barcode + [fixed_right]
```

## Oligos

For each barcode, barcodeEZ builds a **forward** and **reverse** oligo. These carry the assembled barcode plus the appropriate end sequence:

- at an **enzyme boundary**, the enzyme's recognition overhang;
- at an **internal position junction**, the shared 4 bp overhang.

The forward and reverse oligos are reverse complements designed to anneal into a double-stranded insert with the correct sticky ends. They are regenerated automatically whenever the design changes.

## Validation

`validate()` checks every assembled barcode (barcode + any fixed sequences) for unwanted **motifs** — restriction sites or arbitrary sequences that must not appear inside the insert. Any barcode containing a flagged motif is replaced with a fresh draw from the corpus, repeating until it is clean.

If a motif comes from a **fixed sequence or overhang** rather than the barcode itself, no replacement can ever clear it, and `validate()` raises a `RuntimeError` once the corpus is exhausted. See [`validate()`](api.md#validate) for the default motif panel and how to extend it.
