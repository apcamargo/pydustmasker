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
        The input sequence that was provided.
    intervals: tuple of tuples
        A tuple of tuples representing the start and end positions of the
        low-complexity regions identified in the sequence.

    Methods
    -------
    mask
        Returns the sequence with low-complexity regions masked.
    """

    sequence: str
    intervals: tuple[tuple[int, int], ...]

    @property
    def n_masked_bases(self) -> int: ...
    def mask(self, hard: bool = False) -> str:
        """
        Returns the sequence with low-complexity regions masked.

        Parameters
        ----------
        hard : bool, default: False
            If True, low-complexity regions will be masked with 'N' (for
            nucleotide sequences) or 'X' (for protein sequences). By default,
            bases within low-complexity regions are converted to lowercase
            (i.e., soft-masking).
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

        * sequence contains a non-ASCII character
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
        Nucleotide sequence to process. ASCII characters other than 'A', 'C',
        'G', 'T', 'a', 'c', 'g', and 't' are treated as ambiguous bases. Non-ASCII
        characters are rejected. Must contain at least `kmer + 1` bases.
    window_size : int, default: 5000
        Maximum sliding-window size. Larger windows can detect longer repeats
        but use more memory. Must be in the range [`kmer + 1`, 65535].
    score_threshold : float, default: 0.6
        Score threshold for identifying low-complexity regions. Higher values
        result in fewer regions being masked. Must be finite and greater than 0.0.
    kmer : int, default: 7
        The k-mer length used by the Longdust algorithm. Must be in the range
        [1, 12].
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
        Minimum k-mer frequency to start a backward scan. Must be in [2, 65535].
        Lower values are more sensitive but slower.
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
    score_threshold : float
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

        * sequence contains a non-ASCII character
        * sequence length < kmer + 1
        * window_size is not in [kmer + 1, 65535]
        * kmer is not in [1, 12]
        * score_threshold is non-finite or not greater than 0.0
        * min_start_cnt < 2
        * xdrop is 0
        * gc is invalid (not 'auto', None, or float between 0.0 and 1.0)
    TypeError
        If the input parameters are not of the expected type.
    OverflowError
        If an integer cannot be represented by its parameter type, including a
        negative value or `min_start_cnt` above 65535.
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

class TantanMasker(_BaseMasker):
    """
    Identify and mask low-complexity regions and short-period tandem repeats in
    nucleotide or protein sequences using the tantan algorithm.

    Parameters
    ----------
    sequence : str
        The nucleotide or protein sequence to be processed. When `protein` is
        False, characters other than 'A', 'C', 'G', 'T' (case-insensitive) are
        treated as ambiguous residues. When `protein` is True, characters other
        than the standard 20 amino acids (case-insensitive) are considered
        ambiguous.
    protein : bool, default: False
        If True, treat the sequence as a protein sequence; otherwise, treat it as
        a nucleotide sequence.
    repeat_start : float, default: 0.005
        Probability of transitioning from the background state to a repeat state.
        Must be in the range [0.0, 1.0).
    repeat_end : float, default: 0.05
        Probability of ending a repeat (transitioning from a repeat state back to
        the background). Must be in the range [0.0, 1.0].
    decay : float, default: 0.9
        Probability decay from one period offset to the next. Must be in the
        range (0.0, 1.0] and be a normal finite number.
    max_period : int | None, default: None
        Maximum repeat period (cycle length) to consider. If None (default), it
        resolves to 50 for protein sequences and 100 for nucleotide sequences.
        `probabilities` retains O(len(sequence)) data; `repeat_units` uses
        O(sqrt(len(sequence)) * max_period) DP workspace plus O(R) temporary
        candidates for an active tract.
    gap_open : int, default: 0
        Cost of opening a gap within a repeat.
    gap_extend : int | None, default: None
        Cost of extending a gap by one more letter. If None (default), gaps within
        repeats are disabled. If set, must be greater than 0.
    score_threshold : float, default: 0.5
        Posterior probability threshold above which a position is considered part
        of a repeat. Used by `intervals`/`mask`/`probabilities`; has no
        effect on `repeat_units`. Must be in the range [0.0, 1.0].
    min_copy_number : float, default: 2.0
        Minimum estimated copy number (tract length divided by consensus unit
        length) for a tandem repeat tract to be reported by `repeat_units`. Must
        be finite and non-negative. Has no effect on
        `intervals`/`mask`/`probabilities`.

    Attributes
    ----------
    sequence : str
        The sequence that was provided as input.
    protein : bool
        Whether the sequence was treated as a protein sequence.
    repeat_start : float
        Probability of a repeat starting per position.
    repeat_end : float
        Probability of a repeat ending per position.
    decay : float
        Probability decay per unit increase in repeat period.
    max_period : int
        Maximum tandem repeat period (in letters) considered. Resolved from
        `None` to 50 (protein) or 100 (nucleotide) when not provided explicitly.
    gap_open : int
        Cost of opening a gap within a repeat.
    gap_extend : int | None
        Cost of extending a gap by one more letter, or `None` if gaps within
        repeats are disabled.
    score_threshold : float
        Posterior probability threshold used to determine `intervals`/`mask`.
    min_copy_number : float
        Minimum copy number used to filter the tracts returned by `repeat_units`.
    intervals : tuple of tuples
        An immutable tuple of tuples representing the start and end positions of
        the tandem repeat regions identified in the sequence, based on
        per-position posterior probabilities.
    n_masked_bases : int
        The total number of bases/residues that were masked.
    probabilities : tuple of float
        The per-position posterior probability of being part of a tandem repeat,
        one value per character of `sequence`, in order. Computed once at
        construction time.

    Methods
    -------
    mask
        Returns the sequence with tandem repeat regions masked.
    repeat_units
        Returns the consensus tandem repeat unit(s) identified via a Viterbi
        decode, independent of `intervals`.

    Raises
    ------
    ValueError
        If the input parameters violate the following constraints:

        * sequence contains a non-ASCII character
        * repeat_start is not in [0.0, 1.0)
        * repeat_end is not in [0.0, 1.0]
        * decay is not a normal finite number in (0.0, 1.0]
        * score_threshold is not in [0.0, 1.0]
        * max_period is less than 1
        * gap_extend is 0
        * the combined repeat/gap probabilities form an invalid distribution
          (i.e. `repeat_end + 2 * first_gap_prob > 1.0`, where the gap
          probabilities are derived from `gap_open` and `gap_extend`). Raising
          `gap_open` or `gap_extend` lowers `first_gap_prob`, so this is
          resolved by increasing them or by lowering `repeat_end`.
        * min_copy_number is not finite, or is negative
    TypeError
        If the input parameters are not of the expected type.
    OverflowError
        If `max_period`, `gap_open` or `gap_extend` is supplied as a
        negative integer, since they are stored as unsigned values.
    """

    protein: bool
    repeat_start: float
    repeat_end: float
    decay: float
    max_period: int
    gap_open: int
    gap_extend: int | None
    score_threshold: float
    min_copy_number: float
    probabilities: tuple[float, ...]

    def __init__(
        self,
        sequence: str,
        protein: bool = False,
        repeat_start: float = 0.005,
        repeat_end: float = 0.05,
        decay: float = 0.9,
        max_period: int | None = None,
        gap_open: int = 0,
        gap_extend: int | None = None,
        score_threshold: float = 0.5,
        min_copy_number: float = 2.0,
    ) -> None: ...
    def repeat_units(self) -> tuple[tuple[str, int, int, float], ...]:
        """
        Returns consensus tandem repeat units from an independent Viterbi
        decode. Its tract boundaries can differ from `intervals`.

        Returns
        -------
        tuple of (str, int, int, float)
            One `(unit, start, end, copy_number)` tuple per tract, sorted by
            `start`. Gaps can make `copy_number` differ from tract length divided
            by unit length. Tracts below `min_copy_number` are omitted. Ambiguous
            letters in units are reported as 'N' ('X' for protein).
        """
