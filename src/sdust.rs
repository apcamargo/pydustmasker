// Code adapted from: https://crates.io/crates/sdust

use std::collections::VecDeque;
use std::ops::Range;

const MASK: u8 = 63;
const ENCODING_LOOKUP: [u8; 256] = {
    let mut lookup = [4; 256];
    lookup[b'A' as usize] = 0;
    lookup[b'C' as usize] = 1;
    lookup[b'G' as usize] = 2;
    lookup[b'T' as usize] = 3;
    lookup[b'a' as usize] = 0;
    lookup[b'c' as usize] = 1;
    lookup[b'g' as usize] = 2;
    lookup[b't' as usize] = 3;
    lookup
};

#[derive(Debug)]
struct PerfectInterval {
    start: usize,
    finish: usize,
    score: usize,
    l: usize,
}

#[derive(Debug)]
pub struct SymmetricDust<'a> {
    /// `q` in the paper
    sequence: &'a [u8],
    /// The length of the window used by symmetric DUST algorithm
    /// `W` in the paper
    window_size: usize,
    /// 10 times the score threshold used by symmetric DUST algorithm.
    /// `T` in the paper
    score_threshold: usize,
    /// `P` in the paper
    perfect_intervals: VecDeque<PerfectInterval>,
    /// `res` in the paper
    results: Vec<Range<usize>>,
    /// `w` in the paper
    window: VecDeque<usize>,
    // counts in the current window
    cv: [usize; 64],
    cw: [usize; 64],
    // runnings counts
    rv: usize,
    rw: usize,
    /// `L` in the paper
    biggest_num_triplets: usize,
}

impl<'a> SymmetricDust<'a> {
    pub fn process(
        sequence: &'a [u8],
        window_size: usize,
        score_threshold: usize,
    ) -> Vec<(usize, usize)> {
        let mut obj = SymmetricDust {
            sequence,
            window_size,
            score_threshold,
            perfect_intervals: VecDeque::new(),
            results: Vec::new(),
            window: VecDeque::new(),
            cv: [0; 64],
            cw: [0; 64],
            rv: 0,
            rw: 0,
            biggest_num_triplets: 0,
        };

        obj.inner_process();
        let mut res = Vec::with_capacity(obj.results.len());

        // The algorithm can sometimes give end ranges outside of the sequence
        // https://github.com/lh3/sdust/issues/2
        for mut range in obj.results {
            range.end = std::cmp::min(range.end, sequence.len());
            res.push((range.start, range.end));
        }
        res
    }

    fn inner_process(&mut self) {
        // We're going to represent 3 chars in that u8
        let mut triplet: u8 = 0;
        let mut l = 0;
        for i in 0..=self.sequence.len() {
            let b = if i < self.sequence.len() {
                ENCODING_LOOKUP[self.sequence[i] as usize]
            } else {
                4
            };

            // A/T/C/G
            if b < 4 {
                l += 1;
                triplet = (triplet << 2 | b) & MASK;

                // We have at least 3 chars, we can look at them
                if l >= 3 {
                    let mut window_start = if l > self.window_size {
                        l - self.window_size
                    } else {
                        0
                    };
                    window_start += i + 1 - l;

                    self.save_masked_regions(window_start);
                    self.shift_window(triplet as usize);
                    if self.rw * 10 > self.biggest_num_triplets * self.score_threshold {
                        self.find_perfect(window_start);
                    }
                }
            } else {
                // A `N` resets the sequence
                // When we are there (N or end of seq), we empty the intervals found so far
                let mut window_start = if l > self.window_size - 1 {
                    l - self.window_size + 1
                } else {
                    0
                };
                window_start += i + 1 - l;
                while !self.perfect_intervals.is_empty() {
                    window_start += 1;
                    self.save_masked_regions(window_start);
                }

                l = 0;
                triplet = 0;
            }
        }
    }

    /// Save all the intervals that are before the `window_start`
    /// This can only insert one result at a time
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

    /// Add a triplet to the window, shifting all the data to represent the new window
    fn shift_window(&mut self, triplet: usize) {
        let mut s;
        if self.window.len() >= self.window_size - 2 {
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

        if self.cv[triplet] * 10 > 2 * self.score_threshold {
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

    /// Find all the perfect intervals in the window
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
            if new_score * 10 > self.score_threshold * new_l {
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
