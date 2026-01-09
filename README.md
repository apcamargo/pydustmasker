# pydustmasker

`pydustmasker` is a Python library that enables efficient detection and masking of low-complexity regions in nucleotide sequences using the SDUST[^1] and Longdust[^2] algorithms.

## Documentation

The full documentation for `pydustmasker`, including installation instructions, theoretical background, and API reference, is available at <https://apcamargo.github.io/pydustmasker>.

## Installation

Using `pip`:

```sh
pip install pydustmasker
```

Using [`pixi`](https://pixi.sh/):

```sh
# Create a new Pixi workspace and navigate into the workspace directory
pixi init my_workspace && cd my_workspace
# Add Bioconda to the list of channels of your Pixi workspace
pixi workspace channel add bioconda
# Add pydustmasker to your Pixi workspace
pixi add pydustmasker
```

## Usage

To identify and mask low-complexity regions in a nucleotide sequence, create an instance of a masker class and provide your sequence to it. A masker class implements a specific low-complexity detection algorithm and provides methods to retrieve the detected regions and to generate a masked version of the sequence. `pydustmasker` provides two such classes, corresponding to different detection algorithms: [SDUST](https://apcamargo.github.io/pydustmasker/theory#sdust) and [Longdust](https://apcamargo.github.io/pydustmasker/theory#longdust). The SDUST algorithm is implemented in the [`DustMasker`](https://apcamargo.github.io/pydustmasker/api#pydustmasker.DustMasker) class, while the Longdust algorithm is implemented in the [`LongdustMasker`](https://apcamargo.github.io/pydustmasker/api#pydustmasker.LongdustMasker) class.

```py
>>> import pydustmasker
# Example nucleotide sequence
>>> seq = "CGTATATATATAGTATGCGTACTGGGGGGGCT"
# Create a DustMasker object to identify low-complexity regions with the SDUST algorithm
>>> masker = pydustmasker.DustMasker(seq)
# The len() function returns the number of low-complexity regions detected in the sequence
>>> len(masker)
1
# Get the number of bases within low-complexity regions and the intervals of these regions
>>> masker.n_masked_bases
7
>>> masker.intervals
((23, 30),)
# The masker object is iterable, yielding start and end positions of each low-complexity region
>>> for start, end in masker: # (4)!
...     print(f"{start}-{end}: {seq[start:end]}")
23-30: GGGGGGG
```

You can generate a masked version of the sequence using the `mask()` method. By default, low-complexity regions are soft-masked by converting bases to lowercase. Setting the `hard` parameter to `True` enables hard-masking, in which affected bases are replaced with the ambiguous nucleotide `N`.

```py
# The mask() method returns the sequence with low-complexity regions soft-masked
>>> masker.mask()
'CGTATATATATAGTATGCGTACTgggggggCT'
# Hard-masking can be enabled by setting the `hard` parameter to `True`
>>> masker.mask(hard=True)
'CGTATATATATAGTATGCGTACTNNNNNNNCT'
```

The identification of low-complexity regions can be tuned via algorithm-specific parameters. Both `DustMasker` and `LongdustMasker` provide multiple options, documented in the [API reference](https://apcamargo.github.io/pydustmasker/api), that control how low-complexity regions are determined. One shared parameter is `score_threshold`, which controls detection stringency: lowering this threshold results in more regions being classified as low-complexity, whereas increasing it restricts detection to the most clearly low-complexity regions.

```py
# Setting `score_threshold` to 10 results in more low-complexity regions being detected
>>> masker = pydustmasker.DustMasker(seq, score_threshold=10)
>>> len(masker)
2
>>> masker.intervals
((2, 12), (23, 30))
>>> masker.mask()
'CGtatatatataGTATGCGTACTgggggggCT'
```

## Processing sequences in parallel

When working with large numbers of sequences, you can run `pydustmasker` in parallel to process multiple sequences at the same time. This can substantially reduce the total time needed to process all sequences.

The example below uses [Biopython](https://biopython.org/) to parse a FASTA file containing multiple sequences, which are then processed in parallel using a pool of worker processes from the [`multiprocessing`](https://docs.python.org/3/library/multiprocessing.html) module. Each sequence record is submitted to the worker pool via `imap` and processed with `LongdustMasker` to identify low-complexity regions using the Longdust algorithm. The resulting intervals are written to the output file as they become available.

```py
#!/usr/bin/env python

import multiprocessing.pool

from Bio import SeqIO

import pydustmasker

input_file = "sequences.fna"
output_file = "lc_intervals.tsv"


def process_record(record):
    masker = pydustmasker.LongdustMasker(str(record.seq), score_threshold=12)
    return record.id, masker.intervals


if __name__ == "__main__":
    with open(output_file, "w") as f, multiprocessing.pool.Pool() as pool:
        records = SeqIO.parse(input_file, "fasta")
        for name, intervals in pool.imap(process_record, records):
            for start, end in intervals:
                f.write(f"{name}\t{start}\t{end}\n")
```

## References

[^1]: Morgulis, Aleksandr, et al. **A Fast and Symmetric DUST Implementation to Mask Low-Complexity DNA Sequences**. *Journal of Computational Biology*, vol. 13, no. 5, June 2006, pp. 1028–40. <https://doi.org/10.1089/cmb.2006.13.1028>.

[^2]: Li, Heng, and Brian Li. **Finding Low-Complexity DNA Sequences with Longdust**. *arXiv*, 2025. <https://doi.org/10.48550/arxiv.2509.07357>.
