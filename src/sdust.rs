// Code adapted from: https://crates.io/crates/sdust

use crate::common::encode_sequence;
use std::collections::VecDeque;
use std::ops::Range;

#[derive(Debug)]
struct PerfectInterval {
    start: usize,
    finish: usize,
    score: usize,
    l: usize,
}

// Options for SymmetricDust
#[derive(Debug, Clone, Copy)]
pub struct SymmetricDustOptions {
    // The length of the window used by symmetric DUST algorithm.
    // `W` in the paper.
    pub window_size: usize,
    // 10 times the score threshold used by symmetric DUST algorithm.
    // `T` in the paper.
    pub score_threshold: usize,
}

impl Default for SymmetricDustOptions {
    // Provides common default values for the DUST algorithm.
    fn default() -> Self {
        SymmetricDustOptions {
            // Default window size based on typical implementations (e.g., 64)
            window_size: 64,
            // Default score threshold (e.g., 20)
            score_threshold: 20,
        }
    }
}

// The main structure for the Symmetric DUST algorithm execution.
#[derive(Debug)]
pub struct SymmetricDust {
    // Parameters struct
    options: SymmetricDustOptions,
    // `P` in the paper - stores detected intervals that meet the criteria
    perfect_intervals: VecDeque<PerfectInterval>,
    // `res` in the paper - the final, merged results
    results: Vec<Range<usize>>,
    // `w` in the paper - the sliding window of triplets
    window: VecDeque<usize>,
    // counts in the current window
    cv: [usize; 64],
    cw: [usize; 64],
    // runnings counts
    rv: usize,
    rw: usize,
    // `L` in the paper - The biggest number of triplets whose count is <= 2*T/10
    biggest_num_triplets: usize,
}

impl SymmetricDust {
    // Initialize and run the SymmetricDust algorithm on the input sequence
    pub fn process(sequence: &[u8], options: SymmetricDustOptions) -> Vec<(usize, usize)> {
        let mut obj = SymmetricDust {
            options,
            perfect_intervals: VecDeque::new(),
            results: Vec::new(),
            window: VecDeque::new(),
            cv: [0; 64],
            cw: [0; 64],
            rv: 0,
            rw: 0,
            biggest_num_triplets: 0,
        };

        let encoded_seq = encode_sequence(sequence);
        obj.inner_process(&encoded_seq);

        // The algorithm can sometimes give end ranges outside of the sequence
        // https://github.com/lh3/sdust/issues/2
        obj.results
            .into_iter()
            .map(|r| (r.start, r.end.min(sequence.len())))
            .collect()
    }

    fn inner_process(&mut self, encoded_seq: &[u8]) {
        // Mask for triplet extraction (63 = 0b111111)
        const MASK: u8 = 63;
        let mut triplet: u8 = 0;
        let mut l: usize = 0;

        // The chain ensures the loop executes one last time with '4' (non-ACGT) for cleanup
        for (i, &b) in encoded_seq.iter().chain([4].iter()).enumerate() {
            // A/T/C/G
            if b < 4 {
                l += 1;
                triplet = (triplet << 2 | b) & MASK;

                // We have at least 3 chars, we can look at them
                if l >= 3 {
                    // Calculate the starting position in the original sequence
                    let window_start = l.saturating_sub(self.options.window_size) + i + 1 - l;

                    self.save_masked_regions(window_start);
                    self.shift_window(triplet as usize);

                    if self.rw * 10 > self.biggest_num_triplets * self.options.score_threshold {
                        self.find_perfect(window_start);
                    }
                }
            } else {
                // Suggested fix for Ambiguous nucleotides causing end ranges
                // falling outside of the sequence
                // https://github.com/lh3/sdust/issues/2
                // A `N` (or end‐of‐seq) resets the sequence:
                // 1) flush any pending perfect intervals
                let mut window_start = if l > self.options.window_size - 1 {
                    l - self.options.window_size + 1
                } else {
                    0
                };
                window_start += i + 1 - l;

                while !self.perfect_intervals.is_empty() {
                    window_start += 1;
                    self.save_masked_regions(window_start);
                }
                // 2) reset the local context
                l = 0;
                triplet = 0;
                // 3) clear the sliding window and zero out all counts
                self.window.clear();
                self.cw.fill(0);
                self.cv.fill(0);
                self.rw = 0;
                self.rv = 0;
                self.biggest_num_triplets = 0;
            }
        }
    }

    // Save all the intervals that are before the `window_start`
    // This can only insert one result at a time
    fn save_masked_regions(&mut self, window_start: usize) {
        if self.perfect_intervals.is_empty() {
            return;
        }

        let back = self.perfect_intervals.back().unwrap();
        if back.start >= window_start {
            return;
        }

        let num_results = self.results.len();
        // If we already have a result, see if we can merge the last perfect interval with it
        // if they are overlapping
        if num_results > 0 {
            let last_res = &self.results[num_results - 1];
            if back.start <= last_res.end {
                self.results[num_results - 1] =
                    last_res.start..std::cmp::max(last_res.end, back.finish);
            } else {
                self.results.push(back.start..back.finish);
            }
        } else {
            self.results.push(back.start..back.finish);
        }

        while let Some(b) = self.perfect_intervals.back() {
            if b.start < window_start {
                self.perfect_intervals.pop_back();
            } else {
                break;
            }
        }
    }

    // Add a triplet to the window, shifting all the data to represent the new window
    fn shift_window(&mut self, triplet: usize) {
        let mut s;

        if self.window.len() >= self.options.window_size - 2 {
            s = self.window.pop_front().unwrap();
            self.cw[s] -= 1;
            self.rw -= self.cw[s];
            if self.biggest_num_triplets > self.window.len() {
                self.biggest_num_triplets -= 1;
                self.cv[s] -= 1;
                self.rv -= self.cv[s];
            }
        }

        self.window.push_back(triplet);
        self.biggest_num_triplets += 1;

        self.rw += self.cw[triplet];
        self.cw[triplet] += 1;
        self.rv += self.cv[triplet];
        self.cv[triplet] += 1;

        if self.cv[triplet] * 10 > 2 * self.options.score_threshold {
            loop {
                s = self.window[self.window.len() - self.biggest_num_triplets];
                self.biggest_num_triplets -= 1;
                self.cv[s] -= 1;
                self.rv -= self.cv[s];

                if s == triplet {
                    break;
                }
            }
        }
    }

    // Find all the perfect intervals in the window
    fn find_perfect(&mut self, window_start: usize) {
        let mut c = self.cv;
        let mut r = self.rv;
        let mut max_score = 0;
        let mut max_l = 0;

        for i in (0..=self.window.len() - self.biggest_num_triplets - 1).rev() {
            let triplet = self.window[i];
            r += c[triplet];
            c[triplet] += 1;
            let new_score = r;
            let new_l = self.window.len() - i - 1;
            if new_score * 10 > self.options.score_threshold * new_l {
                let mut insertion_position = 0;
                // Figure out where to insert the new interval
                for (j, interval) in self.perfect_intervals.iter().enumerate() {
                    if interval.start < i + window_start {
                        break;
                    }
                    insertion_position = j + 1;
                    if max_score == 0 || interval.score * max_l > max_score * interval.l {
                        max_score = interval.score;
                        max_l = interval.l;
                    }
                }

                // And insert it
                if max_score == 0 || new_score * max_l >= max_score * new_l {
                    max_score = new_score;
                    max_l = new_l;
                    let new_perf = PerfectInterval {
                        start: i + window_start,
                        // +2 => triplet size (3) - 1
                        finish: self.window.len() + 2 + window_start,
                        score: new_score,
                        l: new_l,
                    };

                    self.perfect_intervals.insert(insertion_position, new_perf);
                }
            }
        }
    }
}
