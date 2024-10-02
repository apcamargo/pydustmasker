# pydustmasker

`pydustmasker` is a Python library that provides an efficient implementation of the SDUST algorithm[^1], designed to identify and mask low-complexity regions in nucleotide sequences.

## Usage

`pydustmasker` provides a `DustMasker` class that enables identification of low-complexity regions in an input DNA sequence and mask these regions.

Here is a basic example of how to use `pydustmasker`:

```python
import pydustmasker

# Example nucleotide sequence
masker = pydustmasker.DustMasker("CGTATATATATAGTATGCGTACTGGGGGGGCT")

# Get the low-complexity regions in the sequence and the number of masked bases
>>> print(masker.intervals)
[(23, 30)]
>>> print(masker.n_masked_bases)
7

# The mask() method returns the sequence with low-complexity regions soft-masked
>>> print(masker.mask())
CGTATATATATAGTATGCGTACTgggggggCT

# Hard-masking can be enabled by setting the `hard` parameter to `True`
>>> print(masker.mask(hard=True))
CGTATATATATAGTATGCGTACTNNNNNNNCT

# The `window_size` and `score_threshold` parameters can be adjusted to tune the masking
>>> masker = pydustmasker.DustMasker(
...     "CGTATATATATAGTATGCGTACTGGGGGGGCT",
...     score_threshold=10
... )
>>> print(masker.intervals)
[(2, 12), (23, 30)]
>>> print(masker.mask())
CGtatatatataGTATGCGTACTgggggggCT
```

[^1]: Morgulis, Aleksandr, et al. "[A fast and symmetric DUST implementation to mask low-complexity DNA sequences](https://doi.org/10.1089/cmb.2006.13.1028)". *Journal of Computational Biology* **13.5** (2006): 1028-1040.