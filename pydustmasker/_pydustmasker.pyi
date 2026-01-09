from __future__ import annotations

from collections.abc import Iterator
from typing import overload

__version__: str

class _BaseMasker:
    """
    Base class for sequence masking.

    Attributes
    ----------
    sequence : str
        The nucleotide sequence that was provided as input.
    intervals: tuple of tuples
        A tuple of tuples representing the start and end positions of the
        low-complexity regions identified in the sequence.

    Methods
    -------
    mask
        Returns the sequence with low-complexity regions masked.
    """

    sequence: str
    intervals: tuple[tuple[int, int]]

    @property
    def n_masked_bases(self) -> int: ...
    def mask(self, hard: bool = False) -> str:
        """
        Returns the sequence with low-complexity regions masked.

        Parameters
        ----------
        hard : bool, default: False
            If True, low-complexity regions will be masked with 'N' characters.
            By default, bases within low-complexity regions are converted to
            lowercase (i.e., soft-masking).
        """
        ...
    def __repr__(self) -> str: ...
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[tuple[int, int]]: ...
    @overload
    def __getitem__(self, index: int) -> tuple[int, int]: ...
    @overload
    def __getitem__(self, index: slice) -> tuple[tuple[int, int], ...]: ...

class DustMasker(_BaseMasker):
    """
    Identify and mask low-complexity regions in nucleotide sequences using the
    SDUST algorithm from DustMasker.

    Parameters
    ----------
    sequence : str
        The nucleotide sequence to be processed. Characters other than 'A', 'C',
        'G', 'T', 'a', 'c', 'g', 't' will be considered ambiguous bases.
        The minimum allowed sequence length is 4 bases.
    window_size : int, default: 64
        The length of the window used by symmetric DUST algorithm. The minimum
        allowed value is 4.
    score_threshold : int, default: 20
        Score threshold for identifying low-complexity regions. Higher values
        result in fewer regions being masked.

    Attributes
    ----------
    sequence : str
        The nucleotide sequence that was provided as input.
    window_size : int
        The length of the window used by symmetric DUST algorithm.
    score_threshold : int
        Score threshold for identifying low-complexity regions.
    intervals: tuple of tuples
        A tuple of tuples representing the start and end positions of the
        low-complexity regions identified in the sequence.
    n_masked_bases : int
        The total number of bases that were masked.

    Methods
    -------
    mask
        Returns the sequence with low-complexity regions masked.

    Raises
    ------
    ValueError
        If the input parameters violate the following constraints:

        * sequence length < 4
        * window_size < 4
    TypeError
        If the input parameters are not of the expected type.
    OverflowError
        If a negative integer is passed to `window_size` or `score_threshold`.
    """

    window_size: int
    score_threshold: int

    def __init__(
        self, sequence: str, window_size: int = 64, score_threshold: int = 20
    ) -> None: ...

class LongdustMasker(_BaseMasker):
    """
    Identify and mask low-complexity regions in nucleotide sequences using the
    Longdust algorithm.

    Parameters
    ----------
    sequence : str
        A string representing the nucleotide sequence to be processed. Characters
        other than 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't' will be considered
        ambiguous bases. The minimum allowed sequence length is 4 bases.
    window_size : int, default: 5000
        Maximum size of the sliding window used to scan for low-complexity regions.
        Larger windows can detect longer repeats but increase memory usage. For
        optimal performance, keep window_size < 4^kmer.
    score_threshold : float, default: 0.6
        Score threshold for identifying low-complexity regions. Higher values
        result in fewer regions being masked.
    kmer : int, default: 7
        The k-mer length used by the Longdust algorithm. Must be at least 1.
    gc : float | 'auto' | None, default: None
        GC content for bias correction. If None (default), assume a uniform base
        composition. If 'auto', compute GC from the input sequence. If a float
        between 0.0 and 1.0, use that value.
    xdrop : int | None, default: 50
        Maximum allowable score drop for X-drop extension termination. During
        backward scanning, extension continues as long as (max_score - current_score)
        remains below (score_threshold * xdrop). Once the score drops by more
        than this amount from the peak score observed during the scan, extension
        stops immediately. Lower values enforce stricter extensions and tighter
        boundaries, potentially truncating part of the low-complexity region, whereas
        higher values allow more permissive extensions and looser boundaries, which
        may include non-low-complexity regions. If set to None, X-drop is disabled.
    min_start_cnt : int, default: 3
        Minimum k-mer frequency in the window to trigger a backward scan.
        Only when a k-mer appears at least this many times does the algorithm
        attempt to identify a low-complexity region starting at that position.
        Must be at least 2. Lower values are more sensitive but slower, while
        higher values will result in faster processing but may miss shorter
        repeats.
    approx : bool, default: False
        If True, use approximate mode for guaranteed O(L*w) time complexity.
        In this mode, only the first candidate starting position is examined
        during backward scanning, rather than checking all candidates to find
        the optimal one.
    forward_only : bool, default: False
        If True, only process the forward strand. By default, both strands are processed.

    Attributes
    ----------
    sequence : str
        The nucleotide sequence that was provided as input.
    window_size : int
        The size of the sliding window used to scan for low-complexity regions.
    score_threshold : int
        Score threshold for determining low-complexity regions.
    kmer : int
        k-mer length.
    gc : float | 'auto' | None
        Option used for GC bias correction. Can be None (a uniform base composition
        was assumed), 'auto' (GC was computed from the input sequence), or a float
        between 0.0 and 1.0 (provided by the user).
    xdrop : int | None
        Extension X-drop length.
    min_start_cnt : int
        Minimum k-mer frequency to trigger backward scan.
    approx : bool
        Whether approximate mode was enabled.
    forward_only : bool
        Whether only the forward strand was processed.
    intervals: tuple of tuples
        A tuple of tuples representing the start and end positions of the
        low-complexity regions identified in the sequence.
    n_masked_bases : int
        The total number of bases that were masked.

    Methods
    -------
    mask
        Returns the sequence with low-complexity regions masked.

    Raises
    ------
    ValueError
        If the input parameters violate the following constraints:

        * sequence length < kmer + 1
        * window_size < kmer + 1
        * kmer is 0
        * score_threshold <= 0.0
        * min_start_cnt < 2
        * xdrop is 0
        * gc is invalid (not 'auto', None, or float between 0.0 and 1.0)
    TypeError
        If the input parameters are not of the expected type.
    OverflowError
        If a negative integer is passed to `window_size`, `kmer`, `xdrop`,
        or `min_start_cnt`.
    """

    window_size: int
    score_threshold: float
    kmer: int
    gc: float | str | None
    xdrop: int | None
    min_start_cnt: int
    approx: bool
    forward_only: bool

    def __init__(
        self,
        sequence: str,
        window_size: int = 5000,
        score_threshold: float = 0.6,
        kmer: int = 7,
        gc: float | str | None = None,
        xdrop: int | None = 50,
        min_start_cnt: int = 3,
        approx: bool = False,
        forward_only: bool = False,
    ) -> None: ...
