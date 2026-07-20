# pydustmasker

`pydustmasker` provides a unified Python interface to the SDUST[^1], Longdust[^2], and tantan[^3] algorithms for detecting and masking low-complexity regions and tandem repeats in nucleotide and protein sequences.

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

To identify and mask low-complexity regions in a nucleotide sequence, create an instance of a masker class and provide your sequence to it. A masker class implements a specific low-complexity detection algorithm and provides methods to retrieve the detected regions and to generate a masked version of the sequence. `pydustmasker` provides three such classes, corresponding to different detection algorithms: [SDUST](https://apcamargo.github.io/pydustmasker/theory#sdust), [Longdust](https://apcamargo.github.io/pydustmasker/theory#longdust), and [tantan](https://apcamargo.github.io/pydustmasker/theory#tantan). These algorithms are implemented in three different classes: [`DustMasker`](https://apcamargo.github.io/pydustmasker/api#pydustmasker.DustMasker), [`LongdustMasker`](https://apcamargo.github.io/pydustmasker/api#pydustmasker.LongdustMasker), and [`TantanMasker`](https://apcamargo.github.io/pydustmasker/api#pydustmasker.TantanMasker).

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
>>> for start, end in masker:
...     print(f"{start}-{end}: {seq[start:end]}")
23-30: GGGGGGG
```

### Masking low-complexity regions

You can generate a masked version of the sequence using the `mask()` method. By default, low-complexity regions are soft-masked by converting bases to lowercase. Setting the `hard` parameter to `True` enables hard-masking, in which affected bases are replaced with the ambiguous nucleotide `N`.

```py
# The mask() method returns the sequence with low-complexity regions soft-masked
>>> masker.mask()
'CGTATATATATAGTATGCGTACTgggggggCT'
# Hard-masking can be enabled by setting the `hard` parameter to `True`
>>> masker.mask(hard=True)
'CGTATATATATAGTATGCGTACTNNNNNNNCT'
```

### Tuning the detection of low-complexity regions

The identification of low-complexity regions can be tuned via algorithm-specific parameters. All of the provided masker classes provide multiple parameters, documented in the [API reference](https://apcamargo.github.io/pydustmasker/api), that enable control how low-complexity regions are determined. One shared parameter is `score_threshold`, which controls detection stringency: lowering this threshold results in more regions being classified as low-complexity, whereas increasing it restricts detection to the most clearly low-complexity regions.

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

### Identifying tandem repeats in protein sequences

Although the SDUST and Longdust are specifically designed for nucleotide sequences, tandem repeats in proteins can be identified using the tantan algorithm via the `TantanMasker` class.

```py
# Example protein sequence with an imperfect tandem repeat
>>> prot_seq = "QAEMSTNPKPMSTNPKPMSTDPKPMSTNPKPMSTNPKPMSTNDEH"
# Set protein=True to identify tandem repeats in a protein sequence
>>> masker = pydustmasker.TantanMasker(prot_seq, protein=True)
# Get the intervals of the tandem repeats identified in the sequence
>>> masker.intervals
((9, 38),)
# Hard-mask the tandem repeats in the sequence with 'X' characters
>>> masker.mask(hard=True)
'QAEMSTNPKXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXEH'
```

In addition to masking, `TantanMasker` can determine tandem repeat units through the `repeat_units()` method.
Pass `period` to select candidates with a specific repeat period from the same
decoded repeat tracts.

```py
# Each repeat unit is represented by a (unit, start, end, copy_number) tuple
>>> masker.repeat_units()
(('MSTNPKP', 3, 42, 5.571428571428571),)
>>> masker.repeat_units(period=7)
(('MSTNPKP', 3, 42, 5.571428571428571),)
```

## Processing sequences in parallel

When working with large numbers of sequences, you can run `pydustmasker` in parallel to process multiple sequences at the same time. This can substantially reduce the total time needed to process all sequences.

The example below uses [Biopython](https://biopython.org/) to parse a FASTA file containing multiple sequences, which are then processed in parallel using a pool of worker processes from the [`multiprocessing`](https://docs.python.org/3/library/multiprocessing.html) module. Each sequence record is submitted to the worker pool via `imap` and processed with `LongdustMasker` to identify low-complexity regions using the Longdust algorithm. The resulting intervals are written to the output file as they become available.

```py
#!/usr/bin/env python

import multiprocessing.pool

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import pydustmasker

INPUT_FILE = "sequences.fna"
OUTPUT_FILE = "lc_intervals.tsv"

Intervals = tuple[tuple[int, int], ...]

# The `process_record()` function encapsulates the computation performed on a
# single sequence record. A wrapper like this is necessary because process pools
# distribute work by invoking a single function for each item
def process_record(record: SeqRecord) -> tuple[str, Intervals]:
    masker = pydustmasker.LongdustMasker(str(record.seq))
    return record.id, masker.intervals


if __name__ == "__main__":
    # A process pool is created using `multiprocessing.pool.Pool()`, which
    # automatically manages a set of worker processes. If you are using a
    # free-threaded version of Python, you can also use
    # `multiprocessing.pool.ThreadPool()` to create a pool of threads instead
    # of processes
    with open(OUTPUT_FILE, "w") as f, multiprocessing.pool.Pool() as pool:
        records = SeqIO.parse(INPUT_FILE, "fasta")
        # `SeqIO.parse()` is used to lazily read sequence records from the input
        # FASTA file, avoiding the need to load all sequences into memory at
        # once. Each record is submitted to the process pool via `pool.imap()`,
        # which returns an iterator that yields results in the order they were
        # submitted
        for name, intervals in pool.imap(process_record, records):
            for start, end in intervals:
                f.write(f"{name}\t{start}\t{end}\n")
```

## References

[^1]: Morgulis, Aleksandr, *et al*. **A Fast and Symmetric Dust Implementation to Mask Low-Complexity DNA Sequences**. *Journal of Computational Biology*, vol. 13, no. 5, June 2006, pp. 1028–40. <https://doi.org/10.1089/cmb.2006.13.1028>.

[^2]: Li, Heng, and Brian Li. **Finding Low-Complexity DNA Sequences with Longdust**. *Bioinformatics*, vol. 42, no. 3, Feb. 2026, p. btag112. <https://doi.org/10.1093/bioinformatics/btag112>.

[^3]: Frith, Martin C. **A New Repeat-Masking Method Enables Specific Detection of Homologous Sequences**. *Nucleic Acids Research*, vol. 39, no. 4, Mar. 2011, pp. e23–e23. <https://doi.org/10.1093/nar/gkq1212>.
