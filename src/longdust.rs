use std::collections::VecDeque;
use std::f64::consts::{E, PI};
use std::ops::Range;

const MAX_N: i32 = 10000;

/// Lookup to encode ASCII DNA letters into 0..4
/// A -> 0, C -> 1, G -> 2, T -> 3, others -> 4
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

pub fn encode_sequence(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .map(|&b| ENCODING_LOOKUP[b as usize])
        .collect()
}

/// Reverse complement for encoded sequence:
/// 0(A) <-> 3(T), 1(C) <-> 2(G), 4(N) -> 4(N)
fn reverse_complement_encoded_sequence(encoded_seq: &[u8]) -> Vec<u8> {
    encoded_seq
        .iter()
        .rev()
        .map(|&b| match b {
            0 => 3,
            1 => 2,
            2 => 1,
            3 => 0,
            x => x,
        })
        .collect()
}

/// Candidate forward positions (when backward scan suggests forward check)
#[derive(Clone, Debug)]
struct ForwardPosition {
    pos: usize,
    max_score: f64,
}

/// Options for Longdust
#[derive(Debug, Clone)]
pub struct LongdustOptions {
    pub kmer: usize,
    pub window_size: usize,
    pub threshold: f64,
    pub xdrop_len: usize,
    pub min_start_cnt: usize,
    pub approx: bool,
    pub gc: f64,
}

impl Default for LongdustOptions {
    fn default() -> Self {
        Self {
            kmer: 3,
            window_size: 64,
            threshold: 2.0,
            xdrop_len: 0,
            min_start_cnt: 1,
            approx: true,
            gc: 0.5,
        }
    }
}

#[derive(Debug)]
pub struct Longdust {
    // Input & parameters
    encoded_seq: Vec<u8>,
    opts: LongdustOptions,
    // Precomputed tables
    f: Vec<f64>,
    c: Vec<f64>,
    // Sliding queue: store u32 packed value ((kmer << 1) | ambi)
    q: VecDeque<u32>,
    // Counters used across passes (preallocated)
    ht: Vec<u16>,        // counts used in backward pass (like ld->ht)
    ht_for: Vec<u16>,    // counts used in forward pass
    window_ht: Vec<u16>, // counts of k-mers in the current sliding window (allocated once and reused)
    max_test: usize,
    // Temp for forward candidate positions
    for_pos: Vec<ForwardPosition>,
    // Output intervals
    results: Vec<Range<usize>>,
}

impl Longdust {
    /// Main entry point using the Options struct to resolve the "too many arguments" lint.
    pub fn process(sequence: &[u8], opts: LongdustOptions) -> Vec<(usize, usize)> {
        // Encode sequence
        let encoded_seq = encode_sequence(sequence);
        let table_size = 1usize << (2 * opts.kmer);

        // Calculate f table (probability correction values)
        let f = if opts.gc > 0.0 && opts.gc < 1.0 {
            Self::calculate_f_gc(opts.kmer, opts.window_size + 1, opts.gc)
        } else {
            Self::calculate_f(opts.kmer, opts.window_size + 1)
        };

        // Precompute log values (c)
        let mut c = vec![0.0f64; opts.window_size + 1];
        // Fix for "loop variable i is used to index c": use iter_mut().enumerate()
        for (i, val) in c.iter_mut().enumerate().skip(2) {
            *val = (i as f64).ln();
        }

        // Calculate max_test
        let mut max_test = 0usize;
        let mut s = 0.0f64;
        for (i, &c_val) in c.iter().enumerate().skip(1).take(opts.window_size) {
            s += c_val - opts.threshold;
            let sl = s - f[i];
            if sl > 0.0 {
                max_test = ((i as f64) * (i as f64).ln() / opts.threshold) as usize;
                break;
            }
        }

        // Construct the struct, preallocating tables once (Phase A)
        let mut obj = Self {
            encoded_seq,
            opts, // Store options directly
            f,
            c,
            q: VecDeque::new(),
            ht: vec![0u16; table_size],
            ht_for: vec![0u16; table_size],
            window_ht: vec![0u16; table_size],
            max_test,
            for_pos: Vec::new(),
            results: Vec::new(),
        };

        // Process forward and reverse using same tables (avoid cloning f/c)
        obj.inner_process_both_strands();

        // Convert results into Vec<(usize, usize)>
        let mut out = Vec::with_capacity(obj.results.len());
        for r in obj.results {
            // clamp end to sequence length as in original
            let end = if r.end > sequence.len() {
                sequence.len()
            } else {
                r.end
            };
            out.push((r.start, end));
        }
        out
    }

    /// Process forward and reverse strands, reusing precomputed tables
    fn inner_process_both_strands(&mut self) {
        // Move encoded sequence out to avoid overlapping borrows when calling inner_process
        let seq = std::mem::take(&mut self.encoded_seq);
        // Forward
        self.inner_process(&seq);
        let fwd_intervals = self.results.clone();
        // Reverse
        let len = seq.len();
        let rev_seq = reverse_complement_encoded_sequence(&seq);
        self.results.clear();
        self.inner_process(&rev_seq);
        // Transform reverse intervals back into forward coordinates
        let rev_intervals: Vec<Range<usize>> = self
            .results
            .iter()
            .rev()
            .map(|intv| (len - intv.end)..(len - intv.start))
            .collect();
        // Restore original encoded sequence
        self.encoded_seq = seq;
        // Merge forward and reverse intervals
        self.merge_intervals(fwd_intervals, rev_intervals);
    }

    /// Process a single encoded strand (operates on a slice to avoid cloning tables)
    fn inner_process(&mut self, encoded_seq: &[u8]) {
        let mask = (1u32 << (2 * self.opts.kmer)) - 1;

        self.results.clear();
        self.q.clear();

        // Ensure window_ht is sized appropriately and zero it (reuse allocation across calls)
        let expected = ((mask + 1) as usize).max(1);
        if self.window_ht.len() != expected {
            self.window_ht.resize(expected, 0);
        }
        self.window_ht.fill(0);
        let mut ht_sum = 0.0f64;

        let mut x: u32 = 0;
        let mut l: usize = 0;
        let mut st: i64 = -1;
        let mut en: i64 = -1;
        let mut last_q: i64 = -1;

        // Fix for "loop variable i used to index encoded_seq":
        // We chain the sequence with a single sentinel value (4) to handle the final iteration.
        // This removes the check `if i < len` inside the loop.
        for (i, &b) in encoded_seq.iter().chain(std::iter::once(&4)).enumerate() {
            // Update current k-mer and ambi flag
            let ambi = if b < 4 {
                x = ((x << 2) | (b as u32)) & mask;
                l += 1;
                l < self.opts.kmer
            } else {
                l = 0;
                true
            };

            // Pop front if window is full
            if self.q.len() >= self.opts.window_size {
                let p = self.q.pop_front().unwrap();
                if (p & 1) == 0 {
                    let k = (p >> 1) as usize;
                    if self.window_ht[k] > 0 {
                        // subtract c[count] and decrement
                        ht_sum -= self.c[self.window_ht[k] as usize];
                        self.window_ht[k] -= 1;
                    }
                }

                if last_q == 0 {
                    if (p & 1) == 0 {
                        let k = (p >> 1) as usize;
                        if self.ht[k] > 0 {
                            self.ht[k] -= 1;
                        }
                    }
                } else if last_q > 0 {
                    last_q -= 1;
                }
            }

            let packed = (x << 1) | (if ambi { 1 } else { 0 });
            self.q.push_back(packed);

            if ambi {
                continue;
            }

            let kmer_idx = x as usize;
            self.window_ht[kmer_idx] += 1;
            ht_sum += self.c[self.window_ht[kmer_idx] as usize];

            let mut j: i32 = -1;

            if self.window_ht[kmer_idx] >= self.opts.min_start_cnt {
                let qlen = self.q.len();
                let swin = ht_sum - self.f[qlen] - (qlen as f64) * self.opts.threshold;

                // Attempt extend (only when end matches and some conditions)
                if (i as i64) == en && (last_q == 0 || (i as i64) - st >= qlen as i64) && swin > 0.0
                {
                    j = self.extend();
                }

                // If no extend, check backward possibility and do backward if plausible
                if j < 0 && self.if_backward(self.max_test) {
                    j = self.dust_backward(ht_sum);
                }
            }

            if j >= 0 {
                // Found LCR; compute start of LCR range
                let st2 =
                    (i as i64) - (self.q.len() as i64 - 1 - j as i64) - (self.opts.kmer as i64 - 1);

                if st2 < en {
                    // overlap with active interval
                    if st < 0 || st2 < st {
                        st = st2;
                    }
                } else {
                    // save previous interval and start a new one
                    if st >= 0 {
                        self.save_interval(st as usize, en as usize);
                    }
                    st = st2;
                }
                en = (i + 1) as i64;
                last_q = j as i64;
            }
        }

        if st >= 0 {
            self.save_interval(st as usize, en as usize);
        }
    }

    /// Backward scan to find candidate start positions; returns queue index of start or -1 if none
    fn dust_backward(&mut self, _win_sum: f64) -> i32 {
        let xdrop = self.opts.threshold
            * if self.opts.xdrop_len > 0 {
                self.opts.xdrop_len as f64
            } else {
                self.opts.window_size as f64
            };

        self.ht.fill(0);
        self.for_pos.clear();

        let mut max_i: i32 = -1;
        let mut max_sb: f64 = 0.0;
        let mut last_sl: f64 = -1.0;
        let mut s: f64 = 0.0;
        let mut sw: f64 = 0.0;

        let q_size = self.q.len();
        let mut l: usize = 1;

        // iterate backwards over the queue
        for i in (0..q_size).rev() {
            let x = self.q[i];
            // compute backward score s
            let score_val = if (x & 1) == 0 {
                let k = (x >> 1) as usize;
                self.ht[k] += 1;
                self.c[self.ht[k] as usize]
            } else {
                0.0
            };
            s += score_val - self.opts.threshold;
            let sl = s - self.f[l];
            // compute forward feasibility score sw using the preallocated window_ht to avoid borrowing self immutably while calling other &mut methods
            let sw_val = if (x & 1) == 0 {
                let k = (x >> 1) as usize;
                let idx = self.window_ht[k] + 1 - self.ht[k];
                self.c[idx as usize]
            } else {
                0.0
            };
            sw += sw_val - self.opts.threshold;
            // if forward can't reach, break
            if sw - self.f[l] < 0.0 {
                break;
            }
            // record candidate forward positions where forward pass may be needed
            if sl < last_sl && last_sl > 0.0 && (last_sl - max_sb).abs() < 1e-9 {
                self.for_pos.push(ForwardPosition {
                    pos: i + 1,
                    max_score: max_sb,
                });
            }
            if sl >= max_sb {
                max_sb = sl;
                max_i = i as i32;
            } else if max_i >= 0 && max_sb - sl > xdrop {
                break;
            }
            last_sl = sl;
            l += 1;
        }
        if max_i < 0 {
            return -1;
        }

        // ensure max_i is present in for_pos
        if self.for_pos.is_empty() || max_i < self.for_pos.last().unwrap().pos as i32 {
            self.for_pos.push(ForwardPosition {
                pos: max_i as usize,
                max_score: max_sb,
            });
        }

        // forward examine candidate positions
        let mut max_end: i32 = -1;
        let n_for = self.for_pos.len();
        for idx in (0..n_for).rev() {
            let pos = self.for_pos[idx].pos;
            let max_score = self.for_pos[idx].max_score;
            if (pos as i32) < max_end {
                continue;
            }
            let k = self.dust_forward(pos, max_score);
            if k == (q_size - 1) as i32 {
                return pos as i32;
            }
            if self.opts.approx {
                break;
            }
            max_end = max_end.max(k);
        }
        -1
    }

    /// Forward scan starting at i0; returns index achieving max score or -1
    fn dust_forward(&mut self, i0: usize, max_back: f64) -> i32 {
        self.ht_for.fill(0);
        let mut max_i: i32 = -1;
        let mut max_sf: f64 = 0.0;
        let mut s: f64 = 0.0;
        let mut l: usize = 1;
        for i in i0..self.q.len() {
            let x = self.q[i];
            let score_val = if (x & 1) == 0 {
                let k = (x >> 1) as usize;
                self.ht_for[k] += 1;
                self.c[self.ht_for[k] as usize]
            } else {
                0.0
            };
            s += score_val - self.opts.threshold;

            let sl = s - self.f[l];
            if sl >= max_sf {
                max_sf = sl;
                max_i = i as i32;
            }
            if sl > max_back + 1e-6 {
                break;
            }
            l += 1;
        }

        max_i
    }

    /// Quick backward-check heuristic: returns true if backward scan is worth trying
    fn if_backward(&self, max_step: usize) -> bool {
        let mut s = 0.0;
        for i in (0..self.q.len()).rev().take(max_step) {
            let x = self.q[i];
            let val = if (x & 1) == 0 {
                let k = (x >> 1) as usize;
                self.c[self.window_ht[k] as usize]
            } else {
                0.0
            };
            s += val - self.opts.threshold;
            if s < 0.0 {
                return false;
            }
        }
        true
    }

    /// Try to extend at the last position in queue; returns 0 on success, -1 on fail
    fn extend(&mut self) -> i32 {
        if self.q.is_empty() {
            return -1;
        }
        let x = *self.q.back().unwrap();
        if (x & 1) != 0 {
            return -1;
        }
        let k = (x >> 1) as usize;
        let l = self.q.len().saturating_sub(1);
        let ht_k = *self.ht.get(k).unwrap_or(&0u16);
        let idx = (ht_k as usize).saturating_add(1);
        if idx >= self.c.len() || (l + 1) >= self.f.len() {
            return -1;
        }
        let diff = self.c[idx] - (self.f[l + 1] - self.f[l]);
        if diff < self.opts.threshold {
            return -1;
        }
        if let Some(v) = self.ht.get_mut(k) {
            *v += 1;
        }
        0
    }

    /// Merge two sorted lists of intervals (forward and reverse) into results
    fn merge_intervals(&mut self, fwd: Vec<Range<usize>>, rev: Vec<Range<usize>>) {
        self.results.clear();
        let mut i = 0usize;
        let mut j = 0usize;
        let mut st = 0usize;
        let mut en = 0usize;

        while i < fwd.len() || j < rev.len() {
            let intv: &Range<usize> = if j >= rev.len() {
                i += 1;
                &fwd[i - 1]
            } else if i >= fwd.len() {
                j += 1;
                &rev[j - 1]
            } else if fwd[i].start < rev[j].start {
                i += 1;
                &fwd[i - 1]
            } else {
                j += 1;
                &rev[j - 1]
            };

            if intv.start <= en {
                en = en.max(intv.end);
            } else {
                if en > st {
                    self.results.push(st..en);
                }
                st = intv.start;
                en = intv.end;
            }
        }
        if en > st {
            self.results.push(st..en);
        }
    }

    /// Save an interval into results keeping the sorted/merged invariant
    fn save_interval(&mut self, st: usize, en: usize) {
        let mut k = self.results.len();
        while k > 0 && st <= self.results[k - 1].end {
            k -= 1;
        }
        if k < self.results.len() {
            if st < self.results[k].start {
                self.results[k].start = st;
            }
            if en > self.results[k].end {
                self.results[k].end = en;
            }
            self.results.truncate(k + 1);
        } else {
            self.results.push(st..en);
        }
    }

    // Math helpers for f() table computation (same as reference)
    fn f_large(lambda: f64) -> f64 {
        let x = 0.5 * (2.0 * PI * E * lambda).ln()
            - 1.0 / (12.0 * lambda) * (1.0 + 0.5 / lambda + 19.0 / (30.0 * lambda * lambda));
        x + lambda * (lambda.ln() - 1.0)
    }

    fn calculate_f(k: usize, max_l: usize) -> Vec<f64> {
        let n_kmer = 1i32 << (2 * k);
        let dr = 1.0;
        Self::calculate_f_internal(k, max_l, 1, &[n_kmer], &[dr])
    }

    fn calculate_f_gc(k: usize, max_l: usize, gc: f64) -> Vec<f64> {
        let n_kmer = 1usize << (2 * k);
        let mut dr = vec![0.0f64; k + 1];
        for (i, val) in dr.iter_mut().enumerate().take(k + 1) {
            *val = (gc / 0.5).powi(i as i32) * ((1.0 - gc) / 0.5).powi((k - i) as i32);
        }
        let mut n_dr = vec![0i32; k + 1];
        for x in 0..n_kmer {
            let mut n_gc = 0;
            for i in 0..k {
                let nt = (x >> (2 * i)) & 3;
                if nt == 1 || nt == 2 {
                    n_gc += 1;
                }
            }
            n_dr[n_gc] += 1;
        }
        Self::calculate_f_internal(k, max_l, k + 1, &n_dr, &dr)
    }

    fn calculate_f_internal(
        k: usize,
        max_l: usize,
        nn_dr: usize,
        n_dr: &[i32],
        dr: &[f64],
    ) -> Vec<f64> {
        let n_kmer = 1usize << (2 * k);
        let mut f = vec![0.0f64; max_l + 1];
        // Fix: Use enumerate skip(1) for l indexing
        for (l, val) in f.iter_mut().enumerate().skip(1).take(max_l) {
            let mut accum = 0.0f64;
            for (dr_i, n_dr_i) in dr.iter().zip(n_dr.iter()).take(nn_dr) {
                let lambda = (l as f64) / (n_kmer as f64) * dr_i;
                let fli = if lambda < 30.0 {
                    let mut x = 0.0f64;
                    let mut sn = 0.0f64;
                    let mut y = lambda;
                    for n in 2..=MAX_N {
                        sn += (n as f64).ln();
                        y *= lambda / (n as f64);
                        let z = y * sn;
                        if z < x * f64::EPSILON {
                            break;
                        }
                        x += z;
                    }
                    x * (-lambda).exp()
                } else {
                    Self::f_large(lambda)
                };
                accum += fli * (*n_dr_i as f64);
            }
            *val = accum;
        }
        f
    }
}
