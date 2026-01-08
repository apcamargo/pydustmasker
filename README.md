# pydustmasker

`pydustmasker` is a Python library that enables efficient detection and masking of low-complexity regions in nucleotide sequences using the SDUST[^1] and Longdust[^2] algorithms.

## Usage

`pydustmasker` provides a `DustMasker` class that enables identification of low-complexity regions in an input DNA sequence and mask these regions.

Here is a basic example of how to use `pydustmasker`:

```py
>>> import pydustmasker

# Example nucleotide sequence
>>> seq = "CGTATATATATAGTATGCGTACTGGGGGGGCT"
# Create a DustMasker object to identify low-complexity regions with the SDUST algorithm
>>> sdust_masker = pydustmasker.DustMasker(seq)

# Get the low-complexity regions in the sequence and the number of masked bases
>>> print(sdust_masker.intervals)
((23, 30))
>>> print(sdust_masker.n_masked_bases)
7

# The mask() method returns the sequence with low-complexity regions soft-masked
>>> print(sdust_masker.mask())
CGTATATATATAGTATGCGTACTgggggggCT

# Hard-masking can be enabled by setting the `hard` parameter to `True`
>>> print(sdust_masker.mask(hard=True))
CGTATATATATAGTATGCGTACTNNNNNNNCT

# The `window_size` and `score_threshold` parameters can be adjusted to tune the masking
>>> masker = pydustmasker.DustMasker(seq, score_threshold=10)
>>> print(sdust_masker.intervals)
((2, 12), (23, 30))
>>> print(sdust_masker.mask())
CGtatatatataGTATGCGTACTgggggggCT
```

[^1]: Morgulis, Aleksandr, et al. **A Fast and Symmetric DUST Implementation to Mask Low-Complexity DNA Sequences**. *Journal of Computational Biology*, vol. 13, no. 5, June 2006, pp. 1028–40. <https://doi.org/10.1089/cmb.2006.13.1028>.

[^2]: Li, Heng, and Brian Li. **Finding Low-Complexity DNA Sequences with Longdust**. *arXiv*, 2025. <https://doi.org/10.48550/arxiv.2509.07357>.
