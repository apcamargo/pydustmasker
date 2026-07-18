use crate::common::{encode_with_alphabet, Alphabet};
use std::collections::HashMap;
use std::sync::OnceLock;

const DNA_SIZE: usize = 4;
const PROTEIN_SIZE: usize = 20;

// Match/mismatch matrix for A, C, G, T (match +1, mismatch -1), indexed by the
// encoded letter (0..4).
const DNA_MATCH_MISMATCH: [[i32; DNA_SIZE]; DNA_SIZE] = [
    [1, -1, -1, -1],
    [-1, 1, -1, -1],
    [-1, -1, 1, -1],
    [-1, -1, -1, 1],
];

// BLOSUM62, reordered from NCBI's row/column order (ARNDCQEGHILKMFPSTWYV) into
// this crate's alphabet order (ACDEFGHIKLMNPQRSTVWY). Same values, different
// order.
#[rustfmt::skip]
const BLOSUM62: [[i32; PROTEIN_SIZE]; PROTEIN_SIZE] = [
    [ 4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2], // A
    [ 0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2], // C
    [-2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3], // D
    [-1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2], // E
    [-2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3], // F
    [ 0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3], // G
    [-2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2], // H
    [-1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1], // I
    [-1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2], // K
    [-1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1], // L
    [-1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1], // M
    [-2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2], // N
    [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3], // P
    [-1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1], // Q
    [-1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2], // R
    [ 1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2], // S
    [ 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2], // T
    [ 0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1], // V
    [-3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2], // W
    [-2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7], // Y
];

// Gaussian elimination with partial pivoting, solving `m * v = rhs` in place
// (`v` is overwritten with the solution). Returns false if the system is
// singular.
fn linalg_solve(m: &mut [Vec<f64>], v: &mut [f64]) -> bool {
    let n = v.len();
    for k in 0..n {
        // Partial pivoting: swap in the row with the largest |value| in column k.
        let pivot_row = (k..n)
            .max_by(|&a, &b| m[a][k].abs().total_cmp(&m[b][k].abs()))
            .unwrap();
        if pivot_row != k {
            m.swap(k, pivot_row);
            v.swap(k, pivot_row);
        }
        if m[k][k] == 0.0 {
            return false;
        }
        // Cold path (n <= 21): clone the pivot row to read it while the others
        // are mutated.
        let row_k = m[k].clone();
        for i in 0..n {
            if i == k {
                continue;
            }
            let q = m[i][k] / row_k[k];
            for (mi_j, &mk_j) in m[i][k..n].iter_mut().zip(&row_k[k..n]) {
                *mi_j -= q * mk_j;
            }
            v[i] -= q * v[k];
        }
    }
    for k in 0..n {
        v[k] /= m[k][k];
    }
    true
}

// Largest and smallest value yielded by `values`, in a single pass.
fn max_min(values: impl Iterator<Item = f64>) -> (f64, f64) {
    values.fold((f64::NEG_INFINITY, f64::INFINITY), |(max, min), v| {
        (max.max(v), min.min(v))
    })
}

// Finds a valid upper bound for the lambda search by examining per-row and
// per-column score extremes.
fn find_upper_bound(matrix: &[Vec<f64>]) -> Option<f64> {
    let n = matrix.len();
    let (mut r_max_min, mut c_max_min) = (f64::MAX, f64::MAX);
    let mut l_r = 0usize;
    let mut l_c = 0usize;

    for row in matrix {
        let (r_max, r_min) = max_min(row.iter().copied());
        if r_max == 0.0 && r_min == 0.0 {
            l_r += 1;
        } else if r_max <= 0.0 || r_min >= 0.0 {
            return None;
        } else {
            r_max_min = r_max_min.min(r_max);
        }
    }
    for j in 0..n {
        let (c_max, c_min) = max_min(matrix.iter().map(|row| row[j]));
        if c_max == 0.0 && c_min == 0.0 {
            l_c += 1;
        } else if c_max <= 0.0 || c_min >= 0.0 {
            return None;
        } else {
            c_max_min = c_max_min.min(c_max);
        }
    }
    if l_r == n {
        return None;
    }
    // Empirical 1.1 safety margin to avoid too-tight a bound.
    Some(if r_max_min > c_max_min {
        1.1 * ((n - l_r) as f64).ln() / r_max_min
    } else {
        1.1 * ((n - l_c) as f64).ln() / c_max_min
    })
}

// exp(tau * matrix), solve A^T x = 1, return sum(x).
fn inv_sum(matrix: &[Vec<f64>], tau: f64) -> Option<f64> {
    let n = matrix.len();
    let mut a: Vec<Vec<f64>> = matrix
        .iter()
        .map(|row| row.iter().map(|&s| (tau * s).exp()).collect())
        .collect();
    let mut x = vec![1.0; n];
    linalg_solve(&mut a, &mut x).then(|| x.iter().sum())
}

// Bisects `[ub * 1e-6, ub]` until the bracket contains adjacent doubles.
// Early tolerance-based exits can change near-tie outcomes in both decoders.
fn calculate_lambda(matrix: &[Vec<f64>]) -> f64 {
    let solve = |tau| inv_sum(matrix, tau).expect("well-conditioned matrix");
    let ub = find_upper_bound(matrix).expect("built-in score matrix must yield a valid bound");
    let (mut lo, mut hi) = (ub * 1e-6, ub);
    let (mut lo_sum, mut hi_sum) = (solve(lo), solve(hi));

    // inv_sum is > 1 at lo and < 1 at hi, so the single crossing of
    // inv_sum(tau) - 1 == 0 is always bracketed.
    debug_assert!(lo_sum > 1.0);
    debug_assert!(hi_sum < 1.0);

    while lo_sum != 1.0 && hi_sum != 1.0 {
        let mid = f64::midpoint(lo, hi);
        // The bracket is down to two adjacent doubles, nothing left to test.
        if mid == lo || mid == hi {
            break;
        }
        let mid_sum = solve(mid);
        if mid_sum > 1.0 {
            (lo, lo_sum) = (mid, mid_sum);
        } else {
            (hi, hi_sum) = (mid, mid_sum);
        }
    }

    if (lo_sum - 1.0).abs() < (hi_sum - 1.0).abs() {
        lo
    } else {
        hi
    }
}

// A flat square scoring matrix with one row per encoded letter.
struct Matrix {
    stride: usize,
    values: Vec<f64>,
}

impl Matrix {
    // The scoring row for `letter`, to be indexed by a second encoded letter.
    fn row(&self, letter: u8) -> &[f64] {
        let start = letter as usize * self.stride;
        &self.values[start..start + self.stride]
    }
}

// The scoring tables for one alphabet, all derived from its built-in score
// matrix. Computed at most once per alphabet per process.
struct Tables {
    lambda: f64,
    // exp(lambda * score), used by the forward-backward HMM.
    likelihood: Matrix,
    // lambda * score, used by the Viterbi decoder.
    log_odds: Matrix,
}

fn tables(alphabet: Alphabet) -> &'static Tables {
    static DNA: OnceLock<Tables> = OnceLock::new();
    static PROTEIN: OnceLock<Tables> = OnceLock::new();
    match alphabet {
        Alphabet::Dna => &DNA,
        Alphabet::Protein => &PROTEIN,
    }
    .get_or_init(|| {
        let scores = match alphabet {
            Alphabet::Dna => to_f64_rows(&DNA_MATCH_MISMATCH),
            Alphabet::Protein => to_f64_rows(&BLOSUM62),
        };
        let lambda = calculate_lambda(&scores);
        let log_odds = expand_matrix(&scores, |score| lambda * score);
        Tables {
            lambda,
            likelihood: Matrix {
                stride: log_odds.stride,
                values: log_odds.values.iter().copied().map(f64::exp).collect(),
            },
            log_odds,
        }
    })
}

// Widens a fixed-size integer score matrix into the Vec<Vec<f64>> shape the
// lambda solver works in.
fn to_f64_rows<const N: usize>(matrix: &[[i32; N]; N]) -> Vec<Vec<f64>> {
    matrix
        .iter()
        .map(|row| row.iter().map(|&v| f64::from(v)).collect())
        .collect()
}

// Builds the (size+1) x (size+1) matrix for these raw scores: the last
// row/column is an "ambiguous letter" bucket, filled with `transform` applied
// to the lowest score in the matrix, so an ambiguous letter scores no better
// than the worst real pairing. For DNA that worst pairing is an ordinary
// mismatch (-1), so ambiguous bases are penalized but not excluded.
fn expand_matrix(scores: &[Vec<f64>], transform: impl Fn(f64) -> f64) -> Matrix {
    let size = scores.len();
    let stride = size + 1;
    let min_score = scores.iter().flatten().copied().fold(f64::MAX, f64::min);

    let mut values = Vec::with_capacity(stride * stride);
    for row in scores {
        values.extend(row.iter().chain(&[min_score]).map(|&s| transform(s)));
    }
    values.extend(std::iter::repeat_n(transform(min_score), stride));
    Matrix { stride, values }
}

#[derive(Debug, Clone, Copy)]
pub struct TantanOptions {
    pub alphabet: Alphabet,
    // Background probability of starting a repeat, per position.
    pub repeat_start: f64,
    // Probability of ending a repeat, per position.
    pub repeat_end: f64,
    // Geometric decay of the repeat-start probability with period (1.0 = none).
    pub decay: f64,
    // Maximum tandem repeat period (repeat unit length) to consider.
    pub max_period: usize,
    // Gap-open (existence) cost.
    pub gap_open: u32,
    // Gap-extend (per-letter) cost. None disables gaps.
    pub gap_extend: Option<u32>,
    // Posterior-repeat probability threshold for masking.
    pub score_threshold: f64,
    // Minimum copy number for a tract to be reported by repeat_units().
    pub min_copy_number: f64,
}

impl TantanOptions {
    pub fn gap_probabilities(&self) -> (f64, f64) {
        let Some(extend_cost) = self.gap_extend else {
            return (0.0, 0.0);
        };
        let lambda = tables(self.alphabet).lambda;
        let first_gap_cost = (u64::from(self.gap_open) + u64::from(extend_cost)) as f64;
        let other_gap_prob = (-lambda * extend_cost as f64).exp();
        // The gap-existence cost includes the gap-ending cost.
        let first_gap_prob = (-lambda * first_gap_cost).exp() / (1.0 - other_gap_prob);
        (first_gap_prob, other_gap_prob)
    }

    // A gap needs a nonzero opening probability and a period greater than one.
    fn has_gaps(&self, first_gap_prob: f64) -> bool {
        first_gap_prob > 0.0 && self.max_period > 1
    }
}

// Background-to-foreground transition probability for offset `max_period`,
// given the per-step decay.
fn repeat_offset_prob(prob_mult: f64, max_period: usize) -> f64 {
    if prob_mult != 1.0 {
        (1.0 - prob_mult) / (1.0 - prob_mult.powi(max_period as i32))
    } else {
        1.0 / max_period as f64
    }
}

pub fn probabilities_to_intervals(probabilities: &[f32], threshold: f64) -> Vec<(usize, usize)> {
    let mut intervals = Vec::new();
    let mut start: Option<usize> = None;
    for (i, &p) in probabilities.iter().enumerate() {
        if p as f64 >= threshold {
            start.get_or_insert(i);
        } else if let Some(s) = start.take() {
            intervals.push((s, i));
        }
    }
    if let Some(s) = start {
        intervals.push((s, probabilities.len()));
    }
    intervals
}

#[derive(Debug)]
pub struct RepeatTract {
    pub unit: Vec<u8>,
    pub start: usize,
    pub end: usize,
    pub copy_number: f64,
}

const SCALE_STEP_SIZE: usize = 16;

// How many repeat periods `pos` can look back over: every period up to
// `max_period`, until the start of the sequence cuts it short.
fn max_offset_in_sequence(pos: usize, max_period: usize) -> usize {
    pos.min(max_period)
}

// Transition probabilities shared by both decoders (the Viterbi decoder takes
// their logs), kept in one place so the formulas are not duplicated.
struct TransitionProbs {
    b2b: f64,
    f2b: f64,
    g2g: f64,
    one_gap: f64,
    end_gap: f64,
    f2f0: f64,
    f2f1: f64,
    f2f2: f64,
    // Period-1 background-to-foreground transition. Later periods decay by `decay`.
    b2f_first: f64,
}

impl TransitionProbs {
    fn new(options: &TantanOptions, first_gap_prob: f64, other_gap_prob: f64) -> Self {
        let repeat_end = options.repeat_end;
        Self {
            b2b: 1.0 - options.repeat_start,
            f2b: repeat_end,
            g2g: other_gap_prob,
            one_gap: first_gap_prob * (1.0 - other_gap_prob),
            end_gap: first_gap_prob,
            f2f0: 1.0 - repeat_end,
            f2f1: 1.0 - repeat_end - first_gap_prob,
            f2f2: 1.0 - repeat_end - first_gap_prob * 2.0,
            b2f_first: options.repeat_start * repeat_offset_prob(options.decay, options.max_period),
        }
    }

    // Same probabilities in log space. `b2f_first` keeps only its own log, and
    // the decoder adds the per-period `log(decay)` term itself.
    fn ln(&self) -> Self {
        Self {
            b2b: log_safe(self.b2b),
            f2b: log_safe(self.f2b),
            g2g: log_safe(self.g2g),
            one_gap: log_safe(self.one_gap),
            end_gap: log_safe(self.end_gap),
            f2f0: log_safe(self.f2f0),
            f2f1: log_safe(self.f2f1),
            f2f2: log_safe(self.f2f2),
            b2f_first: log_safe(self.b2f_first),
        }
    }
}

// Forward-backward decoder for the tantan HMM.
struct TantanHmm<'a> {
    likelihood: &'a Matrix,
    max_period: usize,
    has_gaps: bool,

    b2b: f64,
    f2b: f64,
    g2g: f64,
    one_gap_prob: f64,
    end_gap_prob: f64,
    f2f0: f64,
    f2f1: f64,
    f2f2: f64,
    // Background-to-foreground probabilities by increasing period.
    b2f_probs: Vec<f64>,
    foreground_probs: Vec<f64>, // len == max_period
    // Read by the gapped recursion only. len == max_period - 1, or 0 when
    // !has_gaps (which includes max_period == 1, since a gap needs a period to
    // insert into).
    insertion_probs: Vec<f64>,

    background_prob: f64,
    scale_factors: Vec<f64>, // len == seq_len / SCALE_STEP_SIZE
}

impl<'a> TantanHmm<'a> {
    fn new(
        seq_len: usize,
        likelihood: &'a Matrix,
        options: &TantanOptions,
        first_gap_prob: f64,
        other_gap_prob: f64,
    ) -> Self {
        let max_period = options.max_period;
        let tp = TransitionProbs::new(options, first_gap_prob, other_gap_prob);
        let has_gaps = options.has_gaps(first_gap_prob);
        let b2f_probs = std::iter::successors(Some(tp.b2f_first), |p| Some(p * options.decay))
            .take(max_period)
            .collect();

        Self {
            likelihood,
            max_period,
            has_gaps,
            b2b: tp.b2b,
            f2b: tp.f2b,
            g2g: tp.g2g,
            one_gap_prob: tp.one_gap,
            // Only read on the gapped path, which implies max_period > 1.
            end_gap_prob: tp.end_gap,
            f2f0: tp.f2f0,
            f2f1: tp.f2f1,
            f2f2: tp.f2f2,
            b2f_probs,
            foreground_probs: vec![0.0; max_period],
            insertion_probs: vec![0.0; if has_gaps { max_period - 1 } else { 0 }],
            background_prob: 1.0,
            scale_factors: vec![0.0; seq_len / SCALE_STEP_SIZE],
        }
    }

    fn forward_total(&self) -> f64 {
        let from_foreground: f64 = self.foreground_probs.iter().sum();
        self.background_prob * self.b2b + from_foreground * self.f2b
    }

    fn initialize_backward(&mut self) {
        self.background_prob = self.b2b;
        self.foreground_probs.fill(self.f2b);
        self.insertion_probs.fill(0.0);
    }

    fn calc_forward_transition_probs_with_gaps(&mut self) {
        // Only reached when `has_gaps`, which implies max_period > 1. The
        // indexing below relies on it.
        debug_assert!(
            self.max_period > 1,
            "gapped recursion indexes foreground_probs[mp-1] / insertion_probs[mp-2]"
        );
        let mp = self.max_period;
        let b2f_probs = &self.b2f_probs;
        let from_background = self.background_prob * b2f_probs[mp - 1];
        let f_last = self.foreground_probs[mp - 1];
        self.foreground_probs[mp - 1] =
            from_background + f_last * self.f2f1 + self.insertion_probs[mp - 2] * self.end_gap_prob;
        let mut from_foreground = f_last;
        let mut d = f_last;
        for i in (1..mp - 1).rev() {
            let f = self.foreground_probs[i];
            from_foreground += f;
            let prev_ins = self.insertion_probs[i - 1];
            self.foreground_probs[i] = self.background_prob * b2f_probs[i]
                + f * self.f2f2
                + (prev_ins + d) * self.one_gap_prob;
            self.insertion_probs[i] = f + prev_ins * self.g2g;
            d = f + d * self.g2g;
        }
        let f0 = self.foreground_probs[0];
        from_foreground += f0;
        self.foreground_probs[0] =
            self.background_prob * b2f_probs[0] + f0 * self.f2f1 + d * self.end_gap_prob;
        self.insertion_probs[0] = f0;
        self.background_prob = self.background_prob * self.b2b + from_foreground * self.f2b;
    }

    fn calc_backward_transition_probs_with_gaps(&mut self) {
        // Same precondition as the forward direction (see above).
        debug_assert!(
            self.max_period > 1,
            "gapped recursion indexes foreground_probs[mp-1] / insertion_probs[mp-2]"
        );
        let mp = self.max_period;
        let to_background = self.f2b * self.background_prob;
        let f0 = self.foreground_probs[0];
        let mut to_foreground = self.b2f_probs[0] * f0;
        let ins0 = self.insertion_probs[0];
        self.foreground_probs[0] = to_background + self.f2f1 * f0 + ins0;
        let mut d = self.end_gap_prob * f0;
        for i in 1..mp - 1 {
            let f = self.foreground_probs[i];
            to_foreground += self.b2f_probs[i] * f;
            let ins = self.insertion_probs[i];
            self.foreground_probs[i] = to_background + self.f2f2 * f + (ins + d);
            let one_gap_f = self.one_gap_prob * f;
            self.insertion_probs[i - 1] = one_gap_f + self.g2g * ins;
            d = one_gap_f + self.g2g * d;
        }
        let f_last = self.foreground_probs[mp - 1];
        to_foreground += self.b2f_probs[mp - 1] * f_last;
        self.foreground_probs[mp - 1] = to_background + self.f2f1 * f_last + d;
        self.insertion_probs[mp - 2] = self.end_gap_prob * f_last;
        self.background_prob = self.b2b * self.background_prob + to_foreground;
    }

    fn calc_emission_probs(&mut self, sequence: &[u8], pos: usize) {
        let lr_row = self.likelihood.row(sequence[pos]);
        let max_offset = max_offset_in_sequence(pos, self.max_period);
        let letters = &sequence[pos - max_offset..pos];
        let foreground = &mut self.foreground_probs[..max_offset];

        for i in 0..max_offset {
            let letter = letters[max_offset - 1 - i];
            foreground[i] *= lr_row[letter as usize];
        }
        // Periods reaching back past the start of the sequence have no
        // emission to contribute.
        self.foreground_probs[max_offset..].fill(0.0);
    }

    fn calc_forward_transition_and_emission(&mut self, sequence: &[u8], pos: usize) {
        if self.has_gaps {
            self.calc_forward_transition_probs_with_gaps();
            self.calc_emission_probs(sequence, pos);
            return;
        }

        let b = self.background_prob;
        let f2f0 = self.f2f0;
        let lr_row = self.likelihood.row(sequence[pos]);
        let max_offset = max_offset_in_sequence(pos, self.max_period);

        let letters = &sequence[pos - max_offset..pos];
        let foreground = &mut self.foreground_probs[..max_offset];
        let b2f_probs = &self.b2f_probs[..max_offset];

        let mut from_foreground = 0.0;
        for i in 0..max_offset {
            let f = foreground[i];
            from_foreground += f;
            let letter = letters[max_offset - 1 - i];
            foreground[i] = (b * b2f_probs[i] + f * f2f0) * lr_row[letter as usize];
        }
        self.background_prob = b * self.b2b + from_foreground * self.f2b;
    }

    fn calc_emission_and_backward_transition(&mut self, sequence: &[u8], pos: usize) {
        if self.has_gaps {
            self.calc_emission_probs(sequence, pos);
            self.calc_backward_transition_probs_with_gaps();
            return;
        }

        let to_background = self.f2b * self.background_prob;
        let f2f0 = self.f2f0;
        let lr_row = self.likelihood.row(sequence[pos]);
        let max_offset = max_offset_in_sequence(pos, self.max_period);

        let letters = &sequence[pos - max_offset..pos];
        let foreground = &mut self.foreground_probs[..max_offset];
        let b2f_probs = &self.b2f_probs[..max_offset];

        let mut to_foreground = 0.0;
        for i in 0..max_offset {
            let letter = letters[max_offset - 1 - i];
            let f = foreground[i] * lr_row[letter as usize];
            to_foreground += b2f_probs[i] * f;
            foreground[i] = to_background + f2f0 * f;
        }
        self.background_prob = self.b2b * self.background_prob + to_foreground;
    }

    fn rescale(&mut self, scale: f64) {
        self.background_prob *= scale;
        for x in &mut self.foreground_probs {
            *x *= scale;
        }
        for x in &mut self.insertion_probs {
            *x *= scale;
        }
    }

    fn rescale_forward(&mut self, pos: usize) {
        if pos % SCALE_STEP_SIZE == SCALE_STEP_SIZE - 1 {
            let scale = 1.0 / self.background_prob;
            self.scale_factors[pos / SCALE_STEP_SIZE] = scale;
            self.rescale(scale);
        }
    }

    fn rescale_backward(&mut self, pos: usize) {
        if pos % SCALE_STEP_SIZE == SCALE_STEP_SIZE - 1 {
            let scale = self.scale_factors[pos / SCALE_STEP_SIZE];
            self.rescale(scale);
        }
    }

    // Forward pass stores the background probability into the output buffer,
    // then a backward pass computes `1 - (forward_bg[t] * backward_bg[t] / Z)`
    // in place (standard forward-backward posterior decoding).
    fn calc_repeat_probs(&mut self, sequence: &[u8]) -> Vec<f32> {
        let seq_len = sequence.len();
        let mut probs = vec![0.0f32; seq_len];

        for (pos, slot) in probs.iter_mut().enumerate() {
            self.calc_forward_transition_and_emission(sequence, pos);
            self.rescale_forward(pos);
            *slot = self.background_prob as f32;
        }

        let z = self.forward_total();

        self.initialize_backward();
        for pos in (0..seq_len).rev() {
            let non_repeat = probs[pos] as f64 * self.background_prob / z;
            // Round to f32 before subtracting: at a true-zero posterior (e.g.
            // position 0, with no offset to look back to) `non_repeat` sits a
            // few ulps off 1.0, and rounding first lands the result on exactly
            // 0.0 rather than ~1e-16. The clamp guards the other side, where the
            // same ulps would dip below zero and a `score_threshold` of 0.0
            // would miss the position.
            probs[pos] = (1.0 - non_repeat as f32).clamp(0.0, 1.0);
            self.rescale_backward(pos);
            self.calc_emission_and_backward_transition(sequence, pos);
        }

        probs
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ViterbiState {
    Background,
    Repeat(usize),    // period, 1..=max_period
    Insertion(usize), // period, 1..=max_period-1 (only reachable if has_gaps)
}

impl ViterbiState {
    // Maps a state to its physical index into a DP column. This is the column
    // index used to read a score, with no caller-side offset to add.
    fn column_index(self, max_period: usize) -> usize {
        match self {
            ViterbiState::Background => 0,
            ViterbiState::Repeat(period) => period,
            ViterbiState::Insertion(period) => max_period + period,
        }
    }
}

fn log_safe(x: f64) -> f64 {
    if x > 0.0 {
        x.ln()
    } else {
        f64::NEG_INFINITY
    }
}

fn max3(x: f64, y: f64, z: f64) -> f64 {
    x.max(y).max(z)
}

struct ViterbiHmm<'a> {
    log_odds: &'a Matrix,
    max_period: usize,
    has_gaps: bool,

    b2b: f64,
    f2b: f64,
    g2g: f64,
    one_gap_score: f64,
    end_gap_score: f64,
    f2f0: f64,
    f2f1: f64,
    f2f2: f64,
    // Indexed by (period - 1): b2f_scores[period-1] is the background-to-
    // foreground transition score for that period.
    b2f_scores: Vec<f64>,
}

impl<'a> ViterbiHmm<'a> {
    fn new(
        log_odds: &'a Matrix,
        options: &TantanOptions,
        first_gap_prob: f64,
        other_gap_prob: f64,
    ) -> Self {
        let max_period = options.max_period;
        let ln = TransitionProbs::new(options, first_gap_prob, other_gap_prob).ln();

        let log_decay = log_safe(options.decay);
        let b2f_scores: Vec<f64> = (0..max_period)
            .map(|i| ln.b2f_first + i as f64 * log_decay)
            .collect();

        Self {
            log_odds,
            max_period,
            has_gaps: options.has_gaps(first_gap_prob),
            b2b: ln.b2b,
            f2b: ln.f2b,
            g2g: ln.g2g,
            one_gap_score: ln.one_gap,
            // Only read on the gapped path, which implies max_period > 1.
            end_gap_score: ln.end_gap,
            f2f0: ln.f2f0,
            f2f1: ln.f2f1,
            f2f2: ln.f2f2,
            b2f_scores,
        }
    }

    fn column_len(&self) -> usize {
        if self.has_gaps {
            2 * self.max_period
        } else {
            self.max_period + 1
        }
    }

    fn initialize_backward(&self, column: &mut [f64]) {
        let mp = self.max_period;
        column[0] = self.b2b;
        column[1..=mp].fill(self.f2b);
        if self.has_gaps {
            column[mp + 1..2 * mp].fill(f64::NEG_INFINITY);
        }
    }

    // Computes the DP column for position `pos` (best scores-to-go from
    // `pos` to the end), given `old` = the already-computed column for
    // position `pos + 1`. Writes into `new`.
    fn calc_scores(&self, old: &[f64], new: &mut [f64], seq: &[u8], pos: usize) {
        if self.has_gaps {
            self.calc_scores_with_gaps(old, new, seq, pos);
        } else {
            self.calc_scores_no_gaps(old, new, seq, pos);
        }
    }

    fn calc_scores_no_gaps(&self, old: &[f64], new: &mut [f64], seq: &[u8], pos: usize) {
        let max_offset = max_offset_in_sequence(pos, self.max_period);
        let log_row = self.log_odds.row(seq[pos]);
        let to_background = self.f2b + old[0];
        let mut to_foreground = f64::NEG_INFINITY;

        // Index `i` is period `i + 1`, whose letter sits `i + 1` places before
        // `pos`. Binding the three rows as equal-length slices up front keeps
        // the per-period bounds checks out of the loop.
        let letters = &seq[pos - max_offset..pos];
        let old_repeat = &old[1..=max_offset];
        let new_repeat = &mut new[1..=max_offset];
        let b2f_scores = &self.b2f_scores[..max_offset];
        let f2f0 = self.f2f0;

        for i in 0..max_offset {
            let letter = letters[max_offset - 1 - i];
            let f = old_repeat[i] + log_row[letter as usize];
            to_foreground = to_foreground.max(f + b2f_scores[i]);
            new_repeat[i] = to_background.max(f2f0 + f);
        }
        new[(max_offset + 1)..=self.max_period].fill(to_background);
        new[0] = (self.b2b + old[0]).max(to_foreground);
    }

    fn calc_scores_with_gaps(&self, old: &[f64], new: &mut [f64], seq: &[u8], pos: usize) {
        let mp = self.max_period;
        let max_offset = max_offset_in_sequence(pos, self.max_period);
        let log_row = self.log_odds.row(seq[pos]);
        let to_background = self.f2b + old[0];

        // Repeat(i) score-with-emission, for i = 1..max_offset.
        for i in 1..=max_offset {
            new[i] = old[i] + log_row[seq[pos - i] as usize];
        }
        // Periods beyond what the sequence supports have no valid emission.
        new[(max_offset + 1)..=mp].fill(f64::NEG_INFINITY);

        // Rolling insertion/deletion accumulators (same O(1)-extra-work
        // trick as the forward-backward HMM's gapped recursion).
        let mut f = new[1];
        let mut to_foreground = f + self.b2f_scores[0];
        let ins = old[mp + 1]; // old Insertion(1)
        new[1] = max3(to_background, self.f2f1 + f, ins);
        let mut del = self.end_gap_score + f;

        for i in 2..mp {
            f = new[i];
            to_foreground = to_foreground.max(f + self.b2f_scores[i - 1]);
            let prev_ins = old[mp + i]; // old Insertion(i)
            new[i] = max3(to_background, self.f2f2 + f, prev_ins.max(del));
            let one_gap_f = self.one_gap_score + f;
            // new Insertion(i - 1): physical index = (mp - 1) + i = mp + (i - 1)
            new[mp - 1 + i] = one_gap_f.max(self.g2g + prev_ins);
            del = one_gap_f.max(self.g2g + del);
        }

        f = new[mp];
        to_foreground = to_foreground.max(f + self.b2f_scores[mp - 1]);
        new[mp] = max3(to_background, self.f2f1 + f, del);
        new[mp - 1 + mp] = self.end_gap_score + f; // new Insertion(max_period - 1)

        new[0] = (self.b2b + old[0]).max(to_foreground);
    }

    // Given the current `state` and the DP columns at `pos` (`current`)
    // and `pos + 1` (`next`), decides the state at `pos + 1` along the
    // best path.
    fn next_state(
        &self,
        state: ViterbiState,
        current: &[f64],
        next: &[f64],
        seq: &[u8],
        pos: usize,
    ) -> ViterbiState {
        let max_score = current[state.column_index(self.max_period)];
        match state {
            ViterbiState::Background => {
                if self.b2b + next[0] < max_score {
                    self.best_new_repeat_period(next, seq, pos)
                } else {
                    ViterbiState::Background
                }
            }
            ViterbiState::Repeat(period) => {
                if self.f2b + next[0] >= max_score {
                    ViterbiState::Background
                } else if !self.has_gaps {
                    ViterbiState::Repeat(period)
                } else {
                    self.repeat_or_gap_transition(period, max_score, next, seq, pos)
                }
            }
            ViterbiState::Insertion(period) => {
                let next_period = period + 1;
                // Insertion(max_period) does not exist (its column index would
                // be one past the end), so the last period always retires into
                // a repeat. The guard below keeps the column lookup in bounds.
                let stays_inserted = next_period < self.max_period
                    && self.g2g + next[self.max_period + next_period] >= max_score;
                if stays_inserted {
                    ViterbiState::Insertion(next_period)
                } else {
                    ViterbiState::Repeat(next_period)
                }
            }
        }
    }

    // Picks the period maximizing the background-to-foreground transition
    // score at this position.
    fn best_new_repeat_period(&self, next: &[f64], seq: &[u8], pos: usize) -> ViterbiState {
        let log_row = self.log_odds.row(seq[pos]);
        let max_offset = max_offset_in_sequence(pos, self.max_period);

        // Start from -inf and update only on a *strict* increase, so the
        // earliest period wins ties, and an all-`-inf` column (or the empty
        // range at pos == 0) leaves the result as Background.
        let mut best: Option<(usize, f64)> = None;
        for period in 1..=max_offset {
            let score =
                next[period] + log_row[seq[pos - period] as usize] + self.b2f_scores[period - 1];
            if score > f64::NEG_INFINITY && best.is_none_or(|(_, b)| score > b) {
                best = Some((period, score));
            }
        }
        best.map_or(ViterbiState::Background, |(period, _)| {
            ViterbiState::Repeat(period)
        })
    }

    // Decides whether staying in Repeat(period) or diverting into an Insertion
    // state best explains the best path.
    fn repeat_or_gap_transition(
        &self,
        period: usize,
        max_score: f64,
        next: &[f64],
        seq: &[u8],
        pos: usize,
    ) -> ViterbiState {
        let log_row = self.log_odds.row(seq[pos]);
        let f = |p: usize| next[p] + log_row[seq[pos - p] as usize];

        if period == 1 {
            if self.f2f1 + f(1) < max_score {
                ViterbiState::Insertion(1)
            } else {
                ViterbiState::Repeat(1)
            }
        } else if period == self.max_period {
            if self.f2f1 + f(period) < max_score {
                self.best_deletion_period(period, next, seq, pos)
            } else {
                ViterbiState::Repeat(period)
            }
        } else if self.f2f2 + f(period) < max_score {
            if next[self.max_period + period] >= max_score {
                ViterbiState::Insertion(period)
            } else {
                self.best_deletion_period(period, next, seq, pos)
            }
        } else {
            ViterbiState::Repeat(period)
        }
    }

    // Finds which shorter period, reached by "deleting" through periods
    // 2..period, best explains the best path.
    fn best_deletion_period(
        &self,
        period: usize,
        next: &[f64],
        seq: &[u8],
        pos: usize,
    ) -> ViterbiState {
        let log_row = self.log_odds.row(seq[pos]);
        let f = |p: usize| next[p] + log_row[seq[pos - p] as usize];

        let mut best_period = 1;
        // `running` accumulates `g2g` on each step and is overwritten whenever
        // a longer-offset deletion scores higher.
        let mut running = self.end_gap_score + f(1);
        for p in 2..period {
            running += self.g2g;
            let candidate = self.one_gap_score + f(p);
            if candidate > running {
                running = candidate;
                best_period = p;
            }
        }
        ViterbiState::Repeat(best_period)
    }
}

// Runs the backward Viterbi recursion, retaining columns at block boundaries.
fn compute_checkpoints(hmm: &ViterbiHmm, seq: &[u8], block_size: usize) -> Vec<f64> {
    let n = seq.len();
    let col_len = hmm.column_len();
    let num_blocks = n.div_ceil(block_size);
    let mut checkpoints = vec![0.0; (num_blocks + 1) * col_len];

    let last = num_blocks * col_len;
    hmm.initialize_backward(&mut checkpoints[last..last + col_len]);

    let mut cur = checkpoints[last..last + col_len].to_vec();
    let mut prev = vec![0.0; col_len];
    for pos in (0..n).rev() {
        hmm.calc_scores(&cur, &mut prev, seq, pos);
        std::mem::swap(&mut cur, &mut prev);
        if pos % block_size == 0 {
            let off = (pos / block_size) * col_len;
            checkpoints[off..off + col_len].copy_from_slice(&cur);
        }
    }
    checkpoints
}

// Recomputes one block at a time, then traces it forward.
fn traceback(
    hmm: &ViterbiHmm,
    seq: &[u8],
    checkpoints: &[f64],
    block_size: usize,
) -> Vec<RepeatTract> {
    let mut state = ViterbiState::Background;
    let mut accumulator = TractAccumulator::new(hmm.max_period);

    // Reused scratch columns for one block.
    let col_len = hmm.column_len();
    let mut columns = vec![0.0; (block_size + 1) * col_len];
    let num_blocks = checkpoints.len() / col_len - 1;

    for block_idx in 0..num_blocks {
        let lo = block_idx * block_size;
        let hi = ((block_idx + 1) * block_size).min(seq.len());

        let boundary = (block_idx + 1) * col_len;
        let seed = (hi - lo) * col_len;
        columns[seed..seed + col_len].copy_from_slice(&checkpoints[boundary..boundary + col_len]);
        for pos in (lo..hi).rev() {
            let ci = pos - lo;
            let (before, after) = columns.split_at_mut((ci + 1) * col_len);
            hmm.calc_scores(
                &after[..col_len],
                &mut before[ci * col_len..(ci + 1) * col_len],
                seq,
                pos,
            );
        }

        for pos in lo..hi {
            let ci = pos - lo;
            let current = &columns[ci * col_len..(ci + 1) * col_len];
            let next = &columns[(ci + 1) * col_len..(ci + 2) * col_len];
            let new_state = hmm.next_state(state, current, next, seq, pos);
            accumulator.push(state, new_state, pos, seq);
            state = new_state;
        }
    }
    accumulator.finish(state, seq.len(), seq)
}

// Accumulates maximal non-background Viterbi tracts.
struct TractAccumulator {
    max_period: usize,
    tract_start: usize,
    // Whole-unit count and the most recent unit boundary.
    completed_units: usize,
    last_boundary: usize,
    // Repeat unit candidates from the current tract.
    units_seen: Vec<(usize, usize)>,
    // Period frequencies indexed by period.
    period_counts: Vec<usize>,
    tracts: Vec<RepeatTract>,
}

impl TractAccumulator {
    fn new(max_period: usize) -> Self {
        Self {
            max_period,
            tract_start: 0,
            completed_units: 0,
            last_boundary: 0,
            units_seen: Vec::new(),
            period_counts: vec![0; max_period + 1],
            tracts: Vec::new(),
        }
    }

    fn push(&mut self, old_state: ViterbiState, new_state: ViterbiState, pos: usize, seq: &[u8]) {
        match new_state {
            ViterbiState::Background => {
                if old_state != ViterbiState::Background {
                    self.close_tract(pos, old_state, seq);
                }
            }
            ViterbiState::Repeat(period) => {
                if old_state == ViterbiState::Background {
                    self.tract_start = pos - period;
                    self.units_seen.clear();
                    self.completed_units = 0;
                    self.last_boundary = pos - period;
                } else if let ViterbiState::Repeat(old_period) = old_state {
                    // A deletion can complete skipped unit boundaries.
                    for skipped in (period + 1..=old_period).rev() {
                        self.mark_unit_boundary(pos, skipped);
                    }
                }
                self.units_seen.push((pos - period, period));
                self.mark_unit_boundary(pos, period);
            }
            ViterbiState::Insertion(_) => {
                // Insertions extend a tract without adding a repeat unit.
            }
        }
    }

    // Mark each completed unit boundary.
    fn mark_unit_boundary(&mut self, pos: usize, n: usize) {
        if pos - self.last_boundary >= n {
            self.completed_units += 1;
            self.last_boundary = pos;
        }
    }

    fn close_tract(&mut self, end: usize, final_state: ViterbiState, seq: &[u8]) {
        let final_offset = final_state.column_index(self.max_period);
        if let Some(tract) = self.build_tract(end, final_offset, seq) {
            self.tracts.push(tract);
        }
    }

    fn build_tract(&mut self, end: usize, final_offset: usize, seq: &[u8]) -> Option<RepeatTract> {
        self.period_counts.fill(0);
        for &(_, period) in &self.units_seen {
            self.period_counts[period] += 1;
        }
        let mut best_len = 0;
        let mut best_count = 0;
        for (period, &count) in self.period_counts.iter().enumerate() {
            if count > best_count {
                best_count = count;
                best_len = period;
            }
        }
        if best_count == 0 {
            return None;
        }

        // Resolve unit ties by earliest sequence position.
        let mut unit_counts: HashMap<&[u8], (usize, usize)> = HashMap::new();
        for &(unit_start, period) in &self.units_seen {
            if period != best_len {
                continue;
            }
            let unit = &seq[unit_start..unit_start + period];
            unit_counts.entry(unit).or_insert((0, unit_start)).0 += 1;
        }
        let (&best_unit, _) = unit_counts
            .iter()
            .max_by_key(|(_, &(count, first_pos))| (count, std::cmp::Reverse(first_pos)))?;

        Some(RepeatTract {
            unit: best_unit.to_vec(),
            start: self.tract_start,
            end,
            copy_number: self.completed_units as f64
                + (end - self.last_boundary) as f64 / final_offset as f64,
        })
    }

    fn finish(mut self, final_state: ViterbiState, seq_len: usize, seq: &[u8]) -> Vec<RepeatTract> {
        if final_state != ViterbiState::Background {
            self.close_tract(seq_len, final_state, seq);
        }
        self.tracts
    }
}

#[derive(Debug)]
pub struct Tantan {
    options: TantanOptions,
    first_gap_prob: f64,
    other_gap_prob: f64,
}

impl Tantan {
    fn new(options: TantanOptions) -> Self {
        let (first_gap_prob, other_gap_prob) = options.gap_probabilities();
        Self {
            options,
            first_gap_prob,
            other_gap_prob,
        }
    }

    pub fn probabilities(sequence: &[u8], options: TantanOptions) -> Vec<f32> {
        Self::new(options).inner_probabilities(sequence)
    }

    pub fn repeat_units(sequence: &[u8], options: TantanOptions) -> Vec<RepeatTract> {
        Self::new(options).inner_repeat_units(sequence)
    }

    fn inner_probabilities(&self, sequence: &[u8]) -> Vec<f32> {
        let encoded = encode_with_alphabet(sequence, self.options.alphabet);
        let mut hmm = TantanHmm::new(
            encoded.len(),
            &tables(self.options.alphabet).likelihood,
            &self.options,
            self.first_gap_prob,
            self.other_gap_prob,
        );
        hmm.calc_repeat_probs(&encoded)
    }

    fn inner_repeat_units(&self, sequence: &[u8]) -> Vec<RepeatTract> {
        let encoded = encode_with_alphabet(sequence, self.options.alphabet);
        let hmm = ViterbiHmm::new(
            &tables(self.options.alphabet).log_odds,
            &self.options,
            self.first_gap_prob,
            self.other_gap_prob,
        );
        let block_size = (encoded.len() as f64).sqrt().ceil().max(1.0) as usize;
        let checkpoints = compute_checkpoints(&hmm, &encoded, block_size);
        let mut tracts = traceback(&hmm, &encoded, &checkpoints, block_size);
        tracts.retain(|t| t.copy_number >= self.options.min_copy_number);
        tracts
    }
}
