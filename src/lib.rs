mod common;
mod longdust;
mod sdust;
mod tantan;

use crate::common::{decode_sequence, Alphabet};
use crate::longdust::{GcOption, Longdust, LongdustOptions};
use crate::sdust::{SymmetricDust, SymmetricDustOptions};
use crate::tantan::{probabilities_to_intervals, RepeatTract, Tantan, TantanOptions};
use pyo3::{
    exceptions::{PyIndexError, PyTypeError, PyValueError},
    prelude::*,
    types::{PyAny, PySlice, PyTuple},
    IntoPyObjectExt,
};
use std::sync::OnceLock;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum InputError {
    #[error("sequence is too short, must be at least {1} characters long (got {0})")]
    SequenceLengthError(usize, usize),
    #[error("sequence must contain only ASCII characters")]
    NonAsciiSequenceError,
    #[error("invalid window size '{0}', must be at least '{1}'")]
    WindowSizeError(usize, usize),
    #[error("invalid GC content '{0}', must be between 0.0 and 1.0")]
    GcError(f64),
    #[error("invalid k-mer size '{0}', must be greater than 0")]
    KmerSizeError(usize),
    #[error("invalid score threshold '{0}', must be greater than 0.0")]
    LongdustScoreThresholdError(f64),
    #[error("invalid min_start_cnt '{0}', must be at least 2")]
    MinStartCntError(u16),
    #[error("invalid xdrop '{0}', must be at least 1 (or None to disable)")]
    XdropLenError(usize),
    #[error("invalid repeat_start '{0}', must be in the range [0.0, 1.0)")]
    RepeatStartError(f64),
    #[error("invalid repeat_end '{0}', must be in the range [0.0, 1.0]")]
    RepeatEndError(f64),
    #[error("invalid decay '{0}', must be in the range (0.0, 1.0]")]
    DecayError(f64),
    #[error("invalid max_period '{0}', must be in the range [1, {1}]")]
    MaxPeriodError(usize, usize),
    #[error("invalid gap_extend '{0}', must be greater than 0 (or None to disable gaps)")]
    GapExtendError(u32),
    #[error("invalid score_threshold '{0}', must be in the range [0.0, 1.0]")]
    TantanScoreThresholdError(f64),
    #[error(
        "repeat_end plus twice the derived first-gap probability exceeds 1.0; \
         reduce repeat_end, or increase gap_open or gap_extend"
    )]
    GapProbabilityError,
    #[error("invalid min_copy_number '{0}', must be finite and non-negative")]
    MinCopyNumberError(f64),
}

impl From<InputError> for PyErr {
    fn from(err: InputError) -> PyErr {
        PyValueError::new_err(err.to_string())
    }
}

trait Validate {
    fn validate_inputs(&self, sequence: &str) -> Result<(), InputError>;
}

fn validate_base_params(
    sequence: &str,
    window_size: usize,
    min_len: usize,
) -> Result<(), InputError> {
    if !sequence.is_ascii() {
        return Err(InputError::NonAsciiSequenceError);
    }
    if sequence.len() < min_len {
        return Err(InputError::SequenceLengthError(sequence.len(), min_len));
    }
    if window_size < min_len {
        return Err(InputError::WindowSizeError(window_size, min_len));
    }
    Ok(())
}

impl Validate for SymmetricDustOptions {
    fn validate_inputs(&self, sequence: &str) -> Result<(), InputError> {
        validate_base_params(sequence, self.window_size, 4)
    }
}

impl Validate for LongdustOptions {
    fn validate_inputs(&self, sequence: &str) -> Result<(), InputError> {
        let min_len = self.kmer + 1;
        validate_base_params(sequence, self.window_size, min_len)?;

        if self.kmer == 0 {
            return Err(InputError::KmerSizeError(self.kmer));
        }
        if self.score_threshold <= 0.0 {
            return Err(InputError::LongdustScoreThresholdError(
                self.score_threshold,
            ));
        }
        if let GcOption::Fixed(gc_val) = self.gc {
            if !(0.0..=1.0).contains(&gc_val) {
                return Err(InputError::GcError(gc_val));
            }
        }
        if self.min_start_cnt < 2 {
            return Err(InputError::MinStartCntError(self.min_start_cnt));
        }
        if let Some(len) = self.xdrop {
            if len == 0 {
                return Err(InputError::XdropLenError(len));
            }
        }
        Ok(())
    }
}

// Largest period representable by `f64::powi` in `repeat_offset_prob`.
const MAX_PERIOD: usize = i32::MAX as usize;

impl Validate for TantanOptions {
    fn validate_inputs(&self, sequence: &str) -> Result<(), InputError> {
        if !sequence.is_ascii() {
            return Err(InputError::NonAsciiSequenceError);
        }
        if sequence.is_empty() {
            return Err(InputError::SequenceLengthError(0, 1));
        }
        if !(0.0..1.0).contains(&self.repeat_start) {
            return Err(InputError::RepeatStartError(self.repeat_start));
        }
        if !(0.0..=1.0).contains(&self.repeat_end) {
            return Err(InputError::RepeatEndError(self.repeat_end));
        }
        if !(0.0..=1.0).contains(&self.score_threshold) {
            return Err(InputError::TantanScoreThresholdError(self.score_threshold));
        }
        // Keep the supported decay domain finite and well-conditioned.
        if !self.decay.is_normal() || !(0.0..=1.0).contains(&self.decay) {
            return Err(InputError::DecayError(self.decay));
        }
        // `repeat_offset_prob` passes this value to `f64::powi`.
        if self.max_period == 0 || self.max_period > MAX_PERIOD {
            return Err(InputError::MaxPeriodError(self.max_period, MAX_PERIOD));
        }
        if let Some(0) = self.gap_extend {
            return Err(InputError::GapExtendError(0));
        }
        // Foreground transition probabilities must form a valid distribution.
        let (first_gap_prob, _) = self.gap_probabilities();
        if self.repeat_end + first_gap_prob * 2.0 > 1.0 {
            return Err(InputError::GapProbabilityError);
        }
        if !self.min_copy_number.is_finite() || self.min_copy_number < 0.0 {
            return Err(InputError::MinCopyNumberError(self.min_copy_number));
        }
        Ok(())
    }
}

/// Helper to parse the GC parameter from the Python input
fn parse_gc_config(gc: Option<&Bound<'_, PyAny>>) -> PyResult<GcOption> {
    if let Some(obj) = gc {
        if let Ok(val) = obj.extract::<f64>() {
            return Ok(GcOption::Fixed(val));
        } else if let Ok(s) = obj.extract::<String>() {
            if s == "auto" {
                return Ok(GcOption::Auto);
            }
        }
        Err(PyValueError::new_err(
            "gc must be a float between 0.0 and 1.0, 'auto', or None",
        ))
    } else {
        Ok(GcOption::Uniform)
    }
}

#[pyclass]
struct BaseMaskerIter {
    masker: Py<BaseMasker>,
    index: usize,
}

#[pymethods]
impl BaseMaskerIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> Option<(usize, usize)> {
        let masker = self.masker.borrow(py);
        if self.index < masker.intervals.len() {
            let item = masker.intervals[self.index];
            self.index += 1;
            Some(item)
        } else {
            None
        }
    }
}

/// Base class for sequence masking.
///
/// Attributes
/// ----------
/// sequence : str
///     The input sequence that was provided.
/// intervals: tuple of tuples
///     A tuple of tuples representing the start and end positions of the
///     low-complexity regions identified in the sequence.
///
/// Methods
/// -------
/// mask
///     Returns the sequence with low-complexity regions masked.
#[pyclass(subclass, name = "_BaseMasker")]
struct BaseMasker {
    #[pyo3(get)]
    sequence: String,
    intervals: Vec<(usize, usize)>,
    // The letter used for hard-masking.
    mask_symbol: u8,
}

#[pymethods]
impl BaseMasker {
    #[getter]
    fn intervals(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let tuple = PyTuple::new(py, &self.intervals)?;
        Ok(tuple.into_any().unbind())
    }

    #[getter]
    fn n_masked_bases(&self) -> usize {
        self.intervals.iter().map(|(start, end)| end - start).sum()
    }

    fn __len__(&self) -> usize {
        self.intervals.len()
    }

    fn __iter__(slf: PyRef<'_, Self>) -> BaseMaskerIter {
        BaseMaskerIter {
            masker: slf.into(),
            index: 0,
        }
    }

    fn __getitem__(&self, item: Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        let py = item.py();
        let len = self.intervals.len();

        if let Ok(slice) = item.extract::<Bound<'_, PySlice>>() {
            let indices = slice.indices(len.try_into().unwrap())?;
            let mut result = Vec::with_capacity(indices.slicelength);
            let mut i = indices.start;
            for _ in 0..indices.slicelength {
                if i >= 0 && (i as usize) < len {
                    result.push(self.intervals[i as usize]);
                }
                i += indices.step;
            }
            let tuple = PyTuple::new(py, result)?;
            Ok(tuple.into_any().unbind())
        } else if let Ok(idx) = item.extract::<isize>() {
            let mut idx = idx;
            if idx < 0 {
                idx += len as isize;
            }
            if idx < 0 || idx >= len as isize {
                return Err(PyIndexError::new_err("list index out of range"));
            }
            Ok(self.intervals[idx as usize].into_py_any(py)?)
        } else {
            Err(PyTypeError::new_err("indices must be integers or slices"))
        }
    }

    /// Returns the sequence with low-complexity regions masked.
    ///
    /// Parameters
    /// ----------
    /// hard : bool, default: False
    ///     If True, low-complexity regions will be masked with 'N' (for
    ///     nucleotide sequences) or 'X' (for protein sequences). By default,
    ///     bases within low-complexity regions are converted to lowercase
    ///     (i.e., soft-masking).
    #[pyo3(signature = (hard=false))]
    fn mask(&self, hard: bool) -> String {
        let mut masked_sequence = self.sequence.clone();
        let hard_mask = hard.then(|| char::from(self.mask_symbol).to_string());
        for &(start, end) in &self.intervals {
            if let Some(symbol) = &hard_mask {
                let len = end - start;
                masked_sequence.replace_range(start..end, &symbol.repeat(len));
            } else {
                let lowercased = self.sequence[start..end].to_lowercase();
                masked_sequence.replace_range(start..end, &lowercased);
            }
        }
        masked_sequence
    }

    fn __repr__(slf: &Bound<'_, Self>) -> PyResult<String> {
        let class_name = slf.get_type().name()?;
        let inner = slf.borrow();
        let sequence_preview = if inner.sequence.len() > 8 {
            format!("{}…", &inner.sequence[..8])
        } else {
            inner.sequence.clone()
        };

        let mut intervals_repr = String::from("(");
        for (i, (start, end)) in inner.intervals.iter().take(3).enumerate() {
            if i > 0 {
                intervals_repr.push_str(", ");
            }
            intervals_repr.push_str(&format!("({}, {})", start, end));
        }

        if inner.intervals.len() > 3 {
            intervals_repr.push_str(", …");
        }
        intervals_repr.push(')');

        Ok(format!(
            "{}(sequence: '{}', intervals: {})",
            class_name, sequence_preview, intervals_repr
        ))
    }
}

/// Identify and mask low-complexity regions in nucleotide sequences using the
/// SDUST algorithm from DustMasker.
///
/// Parameters
/// ----------
/// sequence : str
///     The nucleotide sequence to be processed. Characters other than 'A', 'C',
///     'G', 'T', 'a', 'c', 'g', 't' will be considered ambiguous bases.
///     The minimum allowed sequence length is 4 bases.
/// window_size : int, default: 64
///     The length of the window used by symmetric DUST algorithm. The minimum
///     allowed value is 4.
/// score_threshold : int, default: 20
///     Score threshold for identifying low-complexity regions. Higher values
///     result in fewer regions being masked.
///
/// Attributes
/// ----------
/// sequence : str
///     The nucleotide sequence that was provided as input.
/// window_size : int
///     The length of the window used by symmetric DUST algorithm.
/// score_threshold : int
///     Score threshold for identifying low-complexity regions.
/// intervals: tuple of tuples
///     A tuple of tuples representing the start and end positions of the
///     low-complexity regions identified in the sequence.
/// n_masked_bases : int
///     The total number of bases that were masked.
///
/// Methods
/// -------
/// mask
///     Returns the sequence with low-complexity regions masked.
///
/// Raises
/// ------
/// ValueError
///    If the input parameters violate the following constraints:
///    * sequence contains a non-ASCII character
///    * sequence length < 4
///    * window_size < 4
/// TypeError
///    If the input parameters are not of the expected type.
/// OverflowError
///    If a negative integer is passed to `window_size` or `score_threshold`.
#[pyclass(extends=BaseMasker)]
struct DustMasker {
    #[pyo3(get)]
    window_size: usize,
    #[pyo3(get)]
    score_threshold: usize,
}

#[pymethods]
impl DustMasker {
    #[new]
    #[pyo3(signature = (sequence, window_size=64, score_threshold=20))]
    fn new(
        sequence: String,
        window_size: usize,
        score_threshold: usize,
    ) -> PyResult<(DustMasker, BaseMasker)> {
        let options = SymmetricDustOptions {
            window_size,
            score_threshold,
        };
        options.validate_inputs(&sequence)?;
        let intervals = SymmetricDust::process(sequence.as_bytes(), options);
        Ok((
            DustMasker {
                window_size,
                score_threshold,
            },
            BaseMasker {
                sequence,
                intervals,
                mask_symbol: b'N',
            },
        ))
    }
}

/// Identify and mask low-complexity regions in nucleotide sequences using the
/// Longdust algorithm.
///
/// Parameters
/// ----------
/// sequence : str
///     A string representing the nucleotide sequence to be processed. Characters
///     other than 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't' will be considered
///     ambiguous bases. The minimum allowed sequence length is 4 bases.
/// window_size : int, default: 5000
///     Maximum size of the sliding window used to scan for low-complexity regions.
///     Larger windows can detect longer repeats but increase memory usage. For
///     optimal performance, keep window_size < 4^kmer.
/// score_threshold : float, default: 0.6
///     Score threshold for identifying low-complexity regions. Higher values
///     result in fewer regions being masked.
/// kmer : int, default: 7
///     The k-mer length used by the Longdust algorithm. Must be at least 1.
/// xdrop : int | None, default: 50
///     Maximum allowable score drop for X-drop extension termination. During
///     backward scanning, extension continues as long as (max_score - current_score)
///     remains below (score_threshold * xdrop). Once the score drops by more
///     than this amount from the peak score observed during the scan, extension
///     stops immediately. Lower values enforce stricter extensions and tighter
///     boundaries, potentially truncating part of the low-complexity region, whereas
///     higher values allow more permissive extensions and looser boundaries, which
///     may include non-low-complexity regions. If set to None, X-drop is disabled.
/// min_start_cnt : int, default: 3
///     Minimum k-mer frequency in the window to trigger a backward scan.
///     Only when a k-mer appears at least this many times does the algorithm
///     attempt to identify a low-complexity region starting at that position.
///     Must be at least 2. Lower values are more sensitive but slower, while
///     higher values will result in faster processing but may miss shorter
///     repeats.
/// approx : bool, default: False
///     If True, use approximate mode for guaranteed O(L*w) time complexity.
///     In this mode, only the first candidate starting position is examined
///     during backward scanning, rather than checking all candidates to find
///     the optimal one.
/// gc : float | 'auto' | None, default: None
///     GC content for bias correction. If None (default), assume a uniform base
///     composition. If 'auto', compute GC from the input sequence. If a float
///     between 0.0 and 1.0, use that value.
/// forward_only : bool, default: False
///     If True, only process the forward strand. By default, both strands are processed.
///
/// Attributes
/// ----------
/// sequence : str
///     The nucleotide sequence that was provided as input.
/// window_size : int
///     The size of the sliding window used to scan for low-complexity regions.
/// score_threshold : int
///     Score threshold for determining low-complexity regions.
/// kmer : int
///     k-mer length.
/// gc : float | 'auto' | None
///     Option used for GC bias correction. Can be None (a uniform base composition
///     was assumed), 'auto' (GC was computed from the input sequence), or a float
///     between 0.0 and 1.0 (provided by the user).
/// xdrop : int | None
///     Extension X-drop length.
/// min_start_cnt : int
///     Minimum k-mer frequency to trigger backward scan.
/// approx : bool
///     Whether approximate mode was enabled.
/// forward_only : bool
///     Whether only the forward strand was processed.
/// intervals: tuple of tuples
///     A tuple of tuples representing the start and end positions of the
///     low-complexity regions identified in the sequence.
/// n_masked_bases : int
///     The total number of bases that were masked.
///
/// Methods
/// -------
/// mask
///     Returns the sequence with low-complexity regions masked.
///
/// Raises
/// ------
/// ValueError
///    If the input parameters violate the following constraints:
///    * sequence contains a non-ASCII character
///    * sequence length < kmer + 1
///    * window_size < kmer + 1
///    * kmer is 0
///    * score_threshold <= 0.0
///    * min_start_cnt < 2
///    * xdrop is 0
///    * gc is invalid (not 'auto', None, or float between 0.0 and 1.0)
/// TypeError
///    If the input parameters are not of the expected type.
/// OverflowError
///    If a negative integer is passed to `window_size`, `kmer`, `xdrop`,
///    or `min_start_cnt`.
#[pyclass(extends=BaseMasker)]
struct LongdustMasker {
    #[pyo3(get)]
    window_size: usize,
    #[pyo3(get)]
    score_threshold: f64,
    #[pyo3(get)]
    kmer: usize,
    gc: GcOption,
    #[pyo3(get)]
    xdrop: Option<usize>,
    #[pyo3(get)]
    min_start_cnt: u16,
    #[pyo3(get)]
    approx: bool,
    #[pyo3(get)]
    forward_only: bool,
}

#[pymethods]
impl LongdustMasker {
    #[new]
    #[pyo3(signature = (
        sequence,
        window_size=5000,
        score_threshold=0.6,
        kmer=7,
        gc=None,
        xdrop=Some(50),
        min_start_cnt=3,
        approx=false,
        forward_only=false
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        sequence: String,
        window_size: usize,
        score_threshold: f64,
        kmer: usize,
        gc: Option<&Bound<'_, PyAny>>,
        xdrop: Option<usize>,
        min_start_cnt: u16,
        approx: bool,
        forward_only: bool,
    ) -> PyResult<(LongdustMasker, BaseMasker)> {
        let gc_config = parse_gc_config(gc)?;
        let options = LongdustOptions {
            window_size,
            score_threshold,
            kmer,
            gc: gc_config,
            xdrop,
            min_start_cnt,
            approx,
            forward_only,
        };

        options.validate_inputs(&sequence)?;

        let intervals = Longdust::process(sequence.as_bytes(), options);

        Ok((
            LongdustMasker {
                window_size,
                score_threshold,
                kmer,
                gc: gc_config,
                xdrop,
                min_start_cnt,
                approx,
                forward_only,
            },
            BaseMasker {
                sequence,
                intervals,
                mask_symbol: b'N',
            },
        ))
    }

    /// Option used for GC bias correction. Can be None (a uniform base composition
    /// was assumed), 'auto' (GC was computed from the input sequence), or a float
    /// between 0.0 and 1.0 (provided by the user).
    #[getter]
    fn gc(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        match &self.gc {
            // Convert f64 to a Python float
            GcOption::Fixed(val) => Ok(val.into_pyobject(py)?.into_any().unbind()),
            // Convert string to a Python str
            GcOption::Auto => Ok("auto".into_pyobject(py)?.into_any().unbind()),
            // Return Python None
            GcOption::Uniform => Ok(py.None()),
        }
    }
}

/// Identify and mask low-complexity regions and short-period tandem repeats
/// in nucleotide or protein sequences using the tantan algorithm.
///
/// Parameters
/// ----------
/// sequence : str
///     The nucleotide or protein sequence to be processed. When `protein` is
///     False, characters other than 'A', 'C', 'G', 'T' (case-insensitive) are
///     treated as ambiguous residues. When `protein` is True, characters other
///     than the standard 20 amino acids (case-insensitive) are considered
///     ambiguous.
/// protein : bool, default: False
///     If True, treat the sequence as a protein sequence; otherwise, treat it as
///     a nucleotide sequence.
/// repeat_start : float, default: 0.005
///     Probability of transitioning from the background state to a repeat state.
///     Must be in the range [0.0, 1.0).
/// repeat_end : float, default: 0.05
///     Probability of ending a repeat (transitioning from a repeat state back to
///     the background). Must be in the range [0.0, 1.0].
/// decay : float, default: 0.9
///     Probability decay from one period offset to the next. Must be in the
///     range (0.0, 1.0] and be a normal finite number.
/// max_period : int | None, default: None
///     Maximum repeat period (cycle length) to consider. If None (default), it
///     resolves to 50 for protein sequences and 100 for nucleotide sequences.
///     `probabilities` retains O(len(sequence)) data; `repeat_units` uses
///     O(sqrt(len(sequence)) * max_period) DP workspace plus O(R) temporary
///     candidates for an active tract.
/// gap_open : int, default: 0
///     Cost of opening a gap within a repeat.
/// gap_extend : int | None, default: None
///     Cost of extending a gap by one more letter. If None (default), gaps within
///     repeats are disabled. If set, must be greater than 0.
/// score_threshold : float, default: 0.5
///     Posterior probability threshold above which a position is considered part
///     of a repeat. Used by `intervals`/`mask`/`probabilities`; has no
///     effect on `repeat_units`. Must be in the range [0.0, 1.0].
/// min_copy_number : float, default: 2.0
///     Minimum estimated copy number (tract length divided by consensus unit
///     length) for a tandem repeat tract to be reported by `repeat_units`. Must
///     be finite and non-negative. Has no effect on
///     `intervals`/`mask`/`probabilities`.
///
/// Attributes
/// ----------
/// sequence : str
///     The sequence that was provided as input.
/// protein : bool
///     Whether the sequence was treated as a protein sequence.
/// repeat_start : float
///     Probability of a repeat starting per position.
/// repeat_end : float
///     Probability of a repeat ending per position.
/// decay : float
///     Probability decay per unit increase in repeat period.
/// max_period : int
///     Maximum tandem repeat period (in letters) considered. Resolved from
///     `None` to 50 (protein) or 100 (nucleotide) when not provided explicitly.
/// gap_open : int
///     Cost of opening a gap within a repeat.
/// gap_extend : int | None
///     Cost of extending a gap by one more letter, or `None` if gaps within
///     repeats are disabled.
/// score_threshold : float
///     Posterior probability threshold used to determine `intervals`/`mask`.
/// min_copy_number : float
///     Minimum copy number used to filter the tracts returned by `repeat_units`.
/// intervals : tuple of tuples
///     An immutable tuple of tuples representing the start and end positions of
///     the tandem repeat regions identified in the sequence, based on
///     per-position posterior probabilities.
/// n_masked_bases : int
///     The total number of bases/residues that were masked.
/// probabilities : tuple of float
///     The per-position posterior probability of being part of a tandem repeat,
///     one value per character of `sequence`, in order. Computed once at
///     construction time.
///
/// Methods
/// -------
/// mask
///     Returns the sequence with tandem repeat regions masked.
/// repeat_units
///     Returns the consensus tandem repeat unit(s) identified via a
///     Viterbi decode, independent of `intervals`.
///
/// Raises
/// ------
/// ValueError
///    If the input parameters violate the following constraints:
///    * sequence contains a non-ASCII character
///    * repeat_start is not in [0.0, 1.0)
///    * repeat_end is not in [0.0, 1.0]
///    * decay is not a normal finite number in (0.0, 1.0]
///    * score_threshold is not in [0.0, 1.0]
///    * max_period is less than 1
///    * gap_extend is 0
///    * the combined repeat/gap probabilities form an invalid distribution
///      (i.e. `repeat_end + 2 * first_gap_prob > 1.0`, where the gap
///      probabilities are derived from `gap_open` and `gap_extend`). Raising
///      `gap_open` or `gap_extend` lowers `first_gap_prob`, so this is
///      resolved by increasing them or by lowering `repeat_end`.
///    * min_copy_number is not finite, or is negative
/// TypeError
///    If the input parameters are not of the expected type.
/// OverflowError
///    If `max_period`, `gap_open` or `gap_extend` is supplied as a negative
///    integer, since they are stored as unsigned values.
#[pyclass(extends=BaseMasker)]
struct TantanMasker {
    // Retained for the lazy repeat unit decode.
    options: TantanOptions,
    probabilities: Vec<f32>,
    // Memoizes the independent Viterbi decode.
    repeat_units: OnceLock<Vec<RepeatTract>>,
}

#[pymethods]
impl TantanMasker {
    #[new]
    #[pyo3(signature = (
        sequence,
        protein=false,
        repeat_start=0.005,
        repeat_end=0.05,
        decay=0.9,
        max_period=None,
        gap_open=0,
        gap_extend=None,
        score_threshold=0.5,
        min_copy_number=2.0
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        sequence: String,
        protein: bool,
        repeat_start: f64,
        repeat_end: f64,
        decay: f64,
        max_period: Option<usize>,
        gap_open: u32,
        gap_extend: Option<u32>,
        score_threshold: f64,
        min_copy_number: f64,
    ) -> PyResult<(TantanMasker, BaseMasker)> {
        let max_period = max_period.unwrap_or(if protein { 50 } else { 100 });
        let alphabet = if protein {
            Alphabet::Protein
        } else {
            Alphabet::Dna
        };
        let options = TantanOptions {
            alphabet,
            repeat_start,
            repeat_end,
            decay,
            max_period,
            gap_open,
            gap_extend,
            score_threshold,
            min_copy_number,
        };
        options.validate_inputs(&sequence)?;

        let probabilities = Tantan::probabilities(sequence.as_bytes(), options);
        let intervals = probabilities_to_intervals(&probabilities, options.score_threshold);

        Ok((
            TantanMasker {
                options,
                probabilities,
                repeat_units: OnceLock::new(),
            },
            BaseMasker {
                sequence,
                intervals,
                mask_symbol: match alphabet {
                    Alphabet::Dna => b'N',
                    Alphabet::Protein => b'X',
                },
            },
        ))
    }

    /// Whether the sequence was interpreted as protein.
    #[getter]
    fn protein(&self) -> bool {
        matches!(self.options.alphabet, Alphabet::Protein)
    }

    /// Probability of a repeat starting per position.
    #[getter]
    fn repeat_start(&self) -> f64 {
        self.options.repeat_start
    }

    /// Probability of ending a repeat per position.
    #[getter]
    fn repeat_end(&self) -> f64 {
        self.options.repeat_end
    }

    /// Probability decay per unit increase in repeat period.
    #[getter]
    fn decay(&self) -> f64 {
        self.options.decay
    }

    /// Maximum tandem repeat period (in letters) considered. Resolved from
    /// `None` to 50 (protein) or 100 (nucleotide) when not provided explicitly.
    #[getter]
    fn max_period(&self) -> usize {
        self.options.max_period
    }

    /// Cost of opening a gap within a repeat.
    #[getter]
    fn gap_open(&self) -> u32 {
        self.options.gap_open
    }

    /// Posterior probability threshold used to determine `intervals`/`mask`.
    #[getter]
    fn score_threshold(&self) -> f64 {
        self.options.score_threshold
    }

    /// Minimum copy number used to filter the tracts returned by `repeat_units`.
    #[getter]
    fn min_copy_number(&self) -> f64 {
        self.options.min_copy_number
    }

    /// Cost of extending a gap by one more letter, or `None` if gaps within
    /// repeats are disabled.
    #[getter]
    fn gap_extend(&self) -> Option<u32> {
        self.options.gap_extend
    }

    /// The per-position posterior probability of being part of a tandem repeat,
    /// one value per character of `sequence`, in order. Computed once at
    /// construction time.
    #[getter]
    fn probabilities(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let tuple = PyTuple::new(py, self.probabilities.iter().map(|&p| p as f64))?;
        Ok(tuple.into_any().unbind())
    }

    /// Returns consensus tandem repeat units from an independent Viterbi
    /// decode. Its tract boundaries can differ from `intervals`.
    ///
    /// Returns
    /// -------
    /// tuple of (str, int, int, float)
    ///     One `(unit, start, end, copy_number)` tuple per tract, sorted by
    ///     `start`. Gaps can make `copy_number` differ from tract length divided
    ///     by unit length. Tracts below `min_copy_number` are omitted. Ambiguous
    ///     letters in units are reported as 'N' ('X' for protein).
    fn repeat_units(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let tracts = slf
            .repeat_units
            .get_or_init(|| Tantan::repeat_units(slf.as_super().sequence.as_bytes(), slf.options));

        let alphabet = slf.options.alphabet;
        let tuple = PyTuple::new(
            py,
            tracts.iter().map(|t| {
                (
                    decode_sequence(&t.unit, alphabet),
                    t.start,
                    t.end,
                    t.copy_number,
                )
            }),
        )?;
        Ok(tuple.into_any().unbind())
    }
}

#[pymodule]
fn _pydustmasker(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<BaseMasker>()?;
    m.add_class::<BaseMaskerIter>()?;
    m.add_class::<DustMasker>()?;
    m.add_class::<LongdustMasker>()?;
    m.add_class::<TantanMasker>()?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
