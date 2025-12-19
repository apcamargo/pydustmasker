mod common;
mod longdust;
mod sdust;

use crate::longdust::{GcOption, Longdust, LongdustOptions};
use crate::sdust::{SymmetricDust, SymmetricDustOptions};
use pyo3::{
    exceptions::{PyIndexError, PyTypeError, PyValueError},
    prelude::*,
    types::{PyAny, PySlice, PyTuple},
    IntoPyObjectExt,
};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum InputError {
    #[error("sequence is too short, must be at least {1} characters long (got {0})")]
    SequenceLengthError(usize, usize),
    #[error("invalid window size '{0}', must be at least '{1}'")]
    WindowSizeError(usize, usize),
    #[error("invalid GC content '{0}', must be between 0.0 and 1.0")]
    GcError(f64),
    #[error("invalid k-mer size '{0}', must be greater than 0")]
    KmerSizeError(usize),
    #[error("invalid score threshold '{0}', must be greater than 0.0")]
    ScoreThresholdError(f64),
    #[error("invalid min_start_cnt '{0}', must be at least 2")]
    MinStartCntError(u16),
    #[error("invalid xdrop '{0}', must be at least 1 (or None to disable)")]
    XdropLenError(usize),
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
            return Err(InputError::ScoreThresholdError(self.score_threshold));
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
        // Borrow the masker to access intervals
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
#[pyclass(subclass, name = "_BaseMasker")]
struct BaseMasker {
    #[pyo3(get)]
    sequence: String,
    intervals: Vec<(usize, usize)>,
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
            // Converts Vec<(usize, usize)> to a Python tuple of tuples
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
    ///     If True, low-complexity regions will be masked with 'N' characters.
    ///     By default, bases within low-complexity regions are converted to
    ///     lowercase (i.e., soft-masking).
    ///
    /// Raises
    /// ------
    /// TypeError
    ///     If the input parameters are not of the expected type.
    #[pyo3(signature = (hard=false))]
    fn mask(&self, hard: bool) -> String {
        let mut masked_sequence = self.sequence.clone();
        for &(start, end) in &self.intervals {
            if hard {
                let len = end - start;
                masked_sequence.replace_range(start..end, &"N".repeat(len));
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
///     Score threshold for identifying low-complexity regions (10 times the
///     actual threshold value). Higher values result in fewer regions being masked.
///
/// Attributes
/// ----------
/// sequence : str
///     The nucleotide sequence that was provided as input.
/// window_size : int
///     The length of the window used by symmetric DUST algorithm.
/// score_threshold : int
///     Score threshold for identifying low-complexity regions.
/// intervals : tuple of tuples
///    An immutable tuple of tuples representing the start and end positions of
///    the low-complexity regions identified in the sequence.
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
///     Score threshold for identifying low-complexity regions.
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
///     may include non–low-complexity regions. If set to None, X-drop is disabled.
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
///    An immutable tuple of tuples representing the start and end positions of
///    the low-complexity regions identified in the sequence.
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
            },
        ))
    }
    /// Expose the resolved GC option as a Python attribute.
    /// Returns a float, 'auto', or None.
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

#[pymodule]
fn _pydustmasker(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<BaseMasker>()?;
    m.add_class::<BaseMaskerIter>()?;
    m.add_class::<DustMasker>()?;
    m.add_class::<LongdustMasker>()?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
