mod sdust;

use crate::sdust::SymmetricDust;
use pyo3::{exceptions::PyValueError, prelude::*};
use thiserror::Error;

const MIN_SEQUENCE_LENGTH: usize = 4;
const MIN_WINDOW_SIZE: usize = 3;

#[derive(Error, Debug)]
pub enum InputError {
    #[error("sequence is too short, it must be at least 4 characters long")]
    SequenceLengthError(usize),
    #[error("invalid window size '{0}', must be at least '3'")]
    WindowSizeError(usize),
}

fn validate_inputs(sequence: &str, window_size: usize) -> Result<(), InputError> {
    if sequence.len() < MIN_SEQUENCE_LENGTH {
        return Err(InputError::SequenceLengthError(sequence.len()));
    }
    if window_size < MIN_WINDOW_SIZE {
        return Err(InputError::WindowSizeError(window_size));
    }
    Ok(())
}

/// Identify and mask low-complexity regions in nucleotide sequences using the
/// symmetric DUST algorithm from DustMasker.
///
/// Parameters
/// ----------
/// sequence : str
///     A string representing the nucleotide sequence to be processed. Characters
///     other than 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't' will be considered
///     ambiguous bases. The minimum allowed sequence length is 4 bases.
/// window_size : int, default: 64
///     The length of the window used by symmetric DUST algorithm. The minimum
///     allowed value is 3.
/// score_threshold : int, default: 20
///     Score threshold for subwindows. The minimum allowed value is 0.
///
/// Attributes
/// ----------
/// sequence : str
///     A string representing the nucleotide sequence that was provided as input.
/// window_size : int
///     The length of the window used by symmetric DUST algorithm.
/// score_threshold : int
///     Score threshold for subwindows.
/// Intervals: list of tuples
///    A immutable list of tuples representing the start and end positions of
///    the low-complexity regions identified in the sequence.
/// n_masked_bases : int
///     The total number of bases that were masked.
///
/// Raises
/// ------
/// ValueError
///    If the input sequence is too short (less than 4 characters) or if the
///    window size is too small (less than 3).
/// TypeError
///    If the input parameters are not of the expected type.
/// OverflowError
///    If a negative integer is passed as the window size or score threshold.
#[pyclass]
struct DustMasker {
    #[pyo3(get)]
    sequence: String,
    #[pyo3(get)]
    window_size: usize,
    #[pyo3(get)]
    score_threshold: usize,
    #[pyo3(get)]
    intervals: Vec<(usize, usize)>,
}

#[pymethods]
impl DustMasker {
    #[new]
    #[pyo3(signature = (sequence, window_size=64, score_threshold=20))]
    fn new(sequence: String, window_size: usize, score_threshold: usize) -> PyResult<DustMasker> {
        validate_inputs(&sequence, window_size)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let intervals = SymmetricDust::process(sequence.as_bytes(), window_size, score_threshold);
        Ok(DustMasker {
            sequence,
            window_size,
            score_threshold,
            intervals,
        })
    }
    #[getter]
    fn n_masked_bases(&self) -> usize {
        self.intervals.iter().map(|(start, end)| end - start).sum()
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
    ///    If the input parameters are not of the expected type.
    #[pyo3(signature = (hard=false))]
    fn mask(&self, hard: bool) -> PyResult<String> {
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
        Ok(masked_sequence)
    }
    fn __repr__(slf: &Bound<'_, Self>) -> PyResult<String> {
        let sequence_preview = if slf.borrow().sequence.len() > 8 {
            format!("{}…", &slf.borrow().sequence[..8])
        } else {
            slf.borrow().sequence.clone()
        };
        Ok(format!(
            "DustMasker(sequence: '{}', intervals: {:?})",
            sequence_preview,
            slf.borrow().intervals
        ))
    }
}

#[pymodule]
fn _pydustmasker(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<DustMasker>()?;
    Ok(())
}
