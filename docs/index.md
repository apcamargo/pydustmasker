---
icon: lucide/rocket
---

# Quick start

`pydustmasker` is a Python library for efficient identification and masking of [low-complexity](theory#low-complexity-sequences) regions in nucleotide sequences. Below, we describe the basic steps required to install and use the library. For a complete description of the available functionality, refer to the [API reference](/api).

## Installation

`pydustmasker` can be installed from PyPI with `pip` or [`uv`](https://docs.astral.sh/uv/), or from Bioconda with [`pixi`](https://pixi.sh/), [`conda`](https://docs.conda.io/projects/conda/), or [`mamba`](https://mamba.readthedocs.io/en/latest/).

=== "pip"

    ```sh
    pip install pydustmasker
    ```

=== "uv"

    ```sh
    uv init my_project && cd my_project # (1)!
    uv add pydustmasker # (2)!
    ```

    1. Create a new Python project with uv and navigate into the project directory.
    2. Add `pydustmasker` as a dependency to your project.

=== "Pixi"

    ```sh
    pixi init my_workspace && cd my_workspace # (1)!
    pixi workspace channel add bioconda # (2)!
    pixi add pydustmasker # (3)!
    ```

    1. Create a new Pixi workspace and navigate into the workspace directory.
    2. Add Bioconda to the list of channels of your Pixi workspace.
    3. Add `pydustmasker` to your Pixi workspace.

=== "Conda"

    ```sh
    conda create -n my_environment -c conda-forge -c bioconda pydustmasker # (1)!
    conda activate my_environment # (2)!
    ```

    1. Create a new Conda environment with `pydustmasker`.
    2. Activate the environment.

=== "Mamba"

    ```sh
    mamba create -n my_environment -c conda-forge -c bioconda pydustmasker # (1)!
    mamba activate my_environment # (2)!
    ```

    1. Create a new Mamba environment with `pydustmasker`.
    2. Activate the environment.

## Usage

To identify and mask low-complexity regions in a nucleotide sequence, create an instance of a masker class and provide your sequence to it. A masker class implements a specific low-complexity detection algorithm and provides methods to retrieve the detected regions and to generate a masked version of the sequence. `pydustmasker` provides two such classes, corresponding to different detection algorithms: [SDUST](theory#sdust) and [Longdust](theory#longdust). The SDUST algorithm is implemented in the [`DustMasker`](api#pydustmasker.DustMasker) class, while the Longdust algorithm is implemented in the [`LongdustMasker`](api#pydustmasker.LongdustMasker) class.

```pycon
>>> import pydustmasker
>>> seq = "CGTATATATATAGTATGCGTACTGGGGGGGCT"
>>> masker = pydustmasker.DustMasker(seq)
>>> len(masker) # (1)!
1
>>> masker.n_masked_bases # (2)!
7
>>> masker.intervals # (3)!
((23, 30),)
>>> for start, end in masker: # (4)!
...     print(f"{start}-{end}: {seq[start:end]}")
23-30: GGGGGGG
```

1. The `len()` function returns the number of low-complexity regions detected in the sequence.
2. The `n_masked_bases` attribute returns the total number of bases within low-complexity regions.
3. The `intervals` attribute returns a tuple of low-complexity regions detected in the sequence, represented as `(start, end)` index pairs.
4. The masker object is iterable, yielding `(start, end)` index pairs for each low-complexity region.

You can generate a masked version of the sequence using the `mask()` method. By default, low-complexity regions are soft-masked by converting bases to lowercase. Setting the `hard` parameter to `True` enables hard-masking, in which affected bases are replaced with the ambiguous nucleotide `N`.

```pycon
>>> masker.mask()
'CGTATATATATAGTATGCGTACTgggggggCT'
>>> masker.mask(hard=True)
'CGTATATATATAGTATGCGTACTNNNNNNNCT'
```

The identification of low-complexity regions can be tuned via algorithm-specific parameters. Both `DustMasker` and `LongdustMasker` provide multiple options, documented in the [API reference](/api), that control how low-complexity regions are determined. One shared parameter is `score_threshold`, which controls detection stringency: lowering this threshold results in more regions being classified as low-complexity, whereas increasing it restricts detection to the most clearly low-complexity regions.

```pycon
>>> masker = pydustmasker.DustMasker(seq, score_threshold=10) # (1)!
>>> len(masker)
2
>>> masker.intervals
((2, 12), (23, 30))
>>> masker.mask()
'CGtatatatataGTATGCGTACTgggggggCT'
```

1. The default `score_threshold` for `DustMasker` is `20`. Setting it to `10` results in more low-complexity regions being detected.

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
