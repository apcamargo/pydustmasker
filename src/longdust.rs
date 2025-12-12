use crate::common::{compute_gc_content, encode_sequence, reverse_complement_encoded_sequence};
use std::collections::VecDeque;
use std::f64::consts::{E, PI};
use std::ops::Range;

// Candidate forward positions (when backward scan suggests forward check)
#[derive(Clone, Debug)]
struct ForwardPosition {
    pos: i32,
    max_score: f64,
}

// Handling the GC parameter
#[derive(Debug, Clone, Copy)]
pub enum GcOption {
    // Use uniform distribution (no GC bias correction)
    Uniform,
    // Use fixed GC value for bias correction
    Fixed(f64),
    // Compute GC from the input sequence
    Auto,
}

// Options for Longdust
#[derive(Debug, Clone, Copy)]
pub struct LongdustOptions {
    // The size of the sliding window
    pub window_size: usize,
    // Score threshold for identifying low-complexity regions
    pub score_threshold: f64,
    // The k-mer size used by the Longdust algorithm
    pub kmer: usize,
    // GC content handling mode
    pub gc: GcOption,
    // X-drop extension length: None=disabled (use full window), Some(n)=limit to n
    pub xdrop: Option<usize>,
    // Minimum k-mer count to trigger backward scan
    pub min_start_cnt: u16,
    // Use approximate mode (faster but less accurate)
    pub approx: bool,
    // Only scan forward strand (skip reverse complement)
    pub forward_only: bool,
}

impl Default for LongdustOptions {
    fn default() -> Self {
        Self {
            window_size: 5000,
            score_threshold: 0.6,
            kmer: 7,
            gc: GcOption::Uniform,
            xdrop: Some(50),
            min_start_cnt: 3,
            approx: false,
            forward_only: false,
        }
    }
}

#[derive(Debug)]
pub struct Longdust {
    // Parameters struct
    options: LongdustOptions,
    // f[L]: precomputed expected score for random sequence of length L,
    // subtracted from raw score to normalize for statistical significance
    f: Vec<f64>,
    // c[i] = ln(i): precomputed log values for efficient ln(count!) computation
    c: Vec<f64>,
    // Sliding window queue of encoded k-mers (VecDeque for O(1) push_back/pop_front)
    // Each element: (kmer_value << 1) | ambiguous_flag
    window: VecDeque<u32>,
    // K-mer count hash tables (size = 4^kmer, preallocated for reuse)
    // ht: temporary count table used during backward/forward scans to track
    // k-mer occurrences within the current scan window
    ht: Vec<u16>,
    // window_ht: count table tracking k-mer occurrences in the sliding window,
    // maintained incrementally as the window advances
    window_ht: Vec<u16>,
    // max_test: maximum number of steps to check in `if_backward()` heuristic,
    // limiting how far back we look before deciding to trigger a full backward scan
    max_test: i32,
    // Temp for forward candidate positions
    for_pos: Vec<ForwardPosition>,
    // Output intervals
    results: Vec<Range<usize>>,
}

impl Longdust {
    // Initialize and run the Longdust algorithm on the input sequence.
    pub fn process(sequence: &[u8], options: LongdustOptions) -> Vec<(usize, usize)> {
        // Encode sequence
        let encoded_seq = encode_sequence(sequence);
        let table_size = 1usize << (2 * options.kmer);

        // Resolve GC mode and calculate the f table
        let f = match options.gc {
            GcOption::Uniform => Self::calculate_f(options.kmer, options.window_size + 1),
            GcOption::Fixed(gc_val) => {
                Self::calculate_f_gc(options.kmer, options.window_size + 1, gc_val)
            }
            GcOption::Auto => {
                // Compute GC from encoded sequence
                let gc_val = compute_gc_content(&encoded_seq);
                Self::calculate_f_gc(options.kmer, options.window_size + 1, gc_val)
            }
        };

        // Precompute log values
        let mut c = vec![0.0f64; options.window_size + 1];
        for (i, val) in c.iter_mut().enumerate().skip(2) {
            *val = (i as f64).ln();
        }

        // Calculate max_test
        let mut max_test = 0i32;
        let mut s = 0.0f64;
        for (i, &c_val) in c.iter().enumerate().skip(1).take(options.window_size) {
            s += c_val - options.score_threshold;
            let sl = s - f[i];
            if sl > 0.0 {
                max_test = ((i as f64) * (i as f64).ln() / options.score_threshold) as i32;
                break;
            }
        }

        // Construct the struct, preallocating tables once
        let mut obj = Self {
            options,
            f,
            c,
            window: VecDeque::new(),
            ht: vec![0u16; table_size],
            window_ht: vec![0u16; table_size],
            max_test,
            for_pos: Vec::new(),
            results: Vec::new(),
        };

        // Process the sequence (either forward-only or both strands)
        if obj.options.forward_only {
            obj.inner_process(&encoded_seq);
        } else {
            obj.inner_process_both_strands(&encoded_seq);
        }

        // Convert results into Vec<(usize, usize)>
        obj.results
            .into_iter()
            .map(|r| (r.start, r.end.min(sequence.len())))
            .collect()
    }

    // Process forward and reverse strands, reusing precomputed tables
    fn inner_process_both_strands(&mut self, encoded_seq: &[u8]) {
        // Forward
        self.inner_process(encoded_seq);
        let fwd_intervals = std::mem::take(&mut self.results);
        // Reverse
        let encoded_seq_rc = reverse_complement_encoded_sequence(encoded_seq);
        self.inner_process(&encoded_seq_rc);
        // Transform reverse intervals back into forward coordinates
        let rev_intervals = self
            .results
            .iter()
            .rev()
            .map(|intv| (encoded_seq.len() - intv.end)..(encoded_seq.len() - intv.start))
            .collect();
        // Merge forward and reverse intervals
        self.merge_intervals(fwd_intervals, rev_intervals);
    }

    // Process a single encoded strand
    fn inner_process(&mut self, encoded_seq: &[u8]) {
        let mask = (1u32 << (2 * self.options.kmer)) - 1;

        self.results.clear();
        self.window.clear();

        // Ensure window_ht is sized appropriately and zero it
        let expected = ((mask + 1) as usize).max(1);
        if self.window_ht.len() != expected {
            self.window_ht.resize(expected, 0);
        }
        self.window_ht.fill(0);
        let mut ht_sum = 0.0f64;
        let mut kmer_val: u32 = 0;
        let mut l: usize = 0;
        let mut start: i64 = -1;
        let mut end: i64 = -1;
        let mut last_q: i64 = -1;

        let len = encoded_seq.len();
        // Main hot loop - process each position plus one sentinel
        for i in 0..=len {
            // Get base, using sentinel value 4 at the end
            let b = if i < len {
                // SAFETY: i < len, so this is in bounds
                unsafe { *encoded_seq.get_unchecked(i) }
            } else {
                4
            };

            // Update current k-mer and ambi flag
            let ambi = if b < 4 {
                kmer_val = ((kmer_val << 2) | (b as u32)) & mask;
                l += 1;
                l < self.options.kmer
            } else {
                l = 0;
                true
            };

            // Pop front if window is full
            if self.window.len() >= self.options.window_size {
                let p = self.window.pop_front().unwrap();
                if (p & 1) == 0 {
                    let k = (p >> 1) as usize;
                    // SAFETY: k is extracted from a packed value (p >> 1) where
                    // p was stored as (x << 1) | ambi_bit. x was masked by
                    // (1 << (2*kmer)) - 1, so x < 2^(2*kmer). Therefore
                    // k < 2^(2*kmer) = window_ht.len()
                    let wht_k = unsafe { *self.window_ht.get_unchecked(k) };
                    if wht_k > 0 {
                        // SAFETY: wht_k is a count in the sliding window, so wht_k <= window_size
                        // and c is sized to window_size + 1
                        ht_sum -= unsafe { *self.c.get_unchecked(wht_k as usize) };
                        // SAFETY: Same bounds as above
                        unsafe {
                            *self.window_ht.get_unchecked_mut(k) = wht_k - 1;
                        }
                    }
                }

                if last_q == 0 {
                    if (p & 1) == 0 {
                        let k = (p >> 1) as usize;
                        // SAFETY: k < 2^(2*kmer) = ht.len()
                        let ht_k = unsafe { *self.ht.get_unchecked(k) };
                        if ht_k > 0 {
                            unsafe {
                                *self.ht.get_unchecked_mut(k) = ht_k - 1;
                            }
                        }
                    }
                } else if last_q > 0 {
                    last_q -= 1;
                }
            }

            let packed = (kmer_val << 1) | (if ambi { 1 } else { 0 });
            self.window.push_back(packed);

            if ambi {
                continue;
            }

            let kmer_idx = kmer_val as usize;
            // SAFETY: kmer_idx = kmer_val, where kmer_val is masked by (1 << (2*kmer)) - 1
            // So kmer_idx < 2^(2*kmer) = window_ht.len()
            let wht_kmer = unsafe { *self.window_ht.get_unchecked(kmer_idx) };
            unsafe {
                *self.window_ht.get_unchecked_mut(kmer_idx) = wht_kmer + 1;
            }

            // SAFETY: wht_kmer + 1 is the new count, which is at most window_size
            // (since we pop elements when queue reaches window_size)
            // and c is sized to window_size + 1
            ht_sum += unsafe { *self.c.get_unchecked((wht_kmer + 1) as usize) };

            let mut j: i32 = -1;

            if wht_kmer + 1 >= self.options.min_start_cnt {
                let qlen = self.window.len();
                // SAFETY: qlen <= window_size (due to pop_front above), and f
                // is sized to window_size + 1
                let f_qlen = unsafe { *self.f.get_unchecked(qlen) };
                let swin = ht_sum - f_qlen - (qlen as f64) * self.options.score_threshold;

                // Attempt extend (only when end matches and some conditions)
                if (i as i64) == end
                    && (last_q == 0 || (i as i64) - start >= qlen as i64)
                    && swin > 0.0
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
                let st2 = (i as i64)
                    - (self.window.len() as i64 - 1 - j as i64)
                    - (self.options.kmer as i64 - 1);

                if st2 < end {
                    // overlap with active interval
                    if start < 0 || st2 < start {
                        start = st2;
                    }
                } else {
                    // save previous interval and start a new one
                    if start >= 0 {
                        self.save_interval(start as usize, end as usize);
                    }
                    start = st2;
                }
                end = (i + 1) as i64;
                last_q = j as i64;
            }
        }

        if start >= 0 {
            self.save_interval(start as usize, end as usize);
        }
    }

    // Backward scan to find candidate start positions; returns queue index of start or -1 if none
    fn dust_backward(&mut self, _win_sum: f64) -> i32 {
        let xdrop = self.options.score_threshold
            * match self.options.xdrop {
                Some(len) => len as f64,
                None => self.options.window_size as f64,
            };

        self.ht.fill(0);
        self.for_pos.clear();

        let mut max_i: i32 = -1;
        let mut max_sb: f64 = 0.0;
        let mut last_sl: f64 = -1.0;
        let mut s: f64 = 0.0;
        let mut sw: f64 = 0.0;

        let q_size = self.window.len() as i32;
        let mut l: usize = 1;

        // Iterate backwards over the queue
        for i in (0..q_size).rev() {
            // SAFETY: i is in range [0, q_size), and q.len() = q_size
            let kmer_val = unsafe { *self.window.get(i as usize).unwrap_unchecked() };

            // Compute backward score s
            let score_val = if (kmer_val & 1) == 0 {
                let k = (kmer_val >> 1) as usize;
                // SAFETY: k < 2^(2*kmer) = ht.len()
                let ht_k = unsafe { *self.ht.get_unchecked(k) };
                let new_ht_k = ht_k + 1;
                unsafe {
                    *self.ht.get_unchecked_mut(k) = new_ht_k;
                }
                // SAFETY: new_ht_k is a count, bounded by queue length <= window_size
                // c is sized to window_size + 1
                unsafe { *self.c.get_unchecked(new_ht_k as usize) }
            } else {
                0.0
            };
            s += score_val - self.options.score_threshold;

            // SAFETY: l is bounded by loop iterations, starting at 1 and incrementing
            // l <= q_size <= window_size, and f is sized to window_size + 1
            let f_l = unsafe { *self.f.get_unchecked(l) };
            let sl = s - f_l;

            // Compute forward feasibility score sw
            let sw_val = if (kmer_val & 1) == 0 {
                let k = (kmer_val >> 1) as usize;
                // SAFETY: k < 2^(2*kmer) = window_ht.len() and ht.len()
                let wht_k = unsafe { *self.window_ht.get_unchecked(k) };
                let ht_k = unsafe { *self.ht.get_unchecked(k) };
                // SAFETY: idx = (wht_k + 1) - ht_k where wht_k <= window_size and ht_k <= wht_k + 1
                // So idx <= window_size, and c is sized to window_size + 1
                unsafe { *self.c.get_unchecked((wht_k + 1 - ht_k) as usize) }
            } else {
                0.0
            };
            sw += sw_val - self.options.score_threshold;

            // If forward can't reach, break
            if sw - f_l < 0.0 {
                break;
            }

            // Record candidate forward positions where forward pass may be needed
            if sl < last_sl && last_sl > 0.0 && (last_sl - max_sb).abs() < 1e-9 {
                self.for_pos.push(ForwardPosition {
                    pos: i + 1,
                    max_score: max_sb,
                });
            }
            if sl >= max_sb {
                max_sb = sl;
                max_i = i;
            } else if max_i >= 0 && max_sb - sl > xdrop {
                break;
            }
            last_sl = sl;
            l += 1;
        }

        if max_i < 0 {
            return -1;
        }

        // Ensure max_i is present in for_pos
        if self.for_pos.is_empty() || max_i < self.for_pos.last().unwrap().pos {
            self.for_pos.push(ForwardPosition {
                pos: max_i,
                max_score: max_sb,
            });
        }

        // Forward examine candidate positions
        let mut max_end: i32 = -1;
        for idx in (0..self.for_pos.len()).rev() {
            // SAFETY: idx < for_pos.len()
            let (pos, max_score) = unsafe {
                let for_pos = self.for_pos.get_unchecked(idx);
                (for_pos.pos, for_pos.max_score)
            };
            if pos < max_end {
                continue;
            }
            let k = self.dust_forward(pos, max_score);
            if k == (q_size - 1) {
                return pos;
            }
            // If `self.options.approx` is enabled, we act greedily and accept the
            // result of this first candidate immediately. We skip checking the
            // remaining candidate positions in 'for_pos'. This guarantees O(L*w)
            // performance by limiting the inner loop to 1 iteration, but we might
            // miss a theoretically higher score that started at an earlier position.
            if self.options.approx {
                break;
            }
            max_end = max_end.max(k);
        }
        -1
    }

    // Forward scan starting at i0; returns index achieving max score or -1
    fn dust_forward(&mut self, i0: i32, max_back: f64) -> i32 {
        self.ht.fill(0);
        let mut max_i: i32 = -1;
        let mut max_sf: f64 = 0.0;
        let mut s: f64 = 0.0;
        let mut l: usize = 1;
        let q_len = self.window.len();
        for i in (i0 as usize)..q_len {
            // SAFETY: i is in range [i0, q_len), verified by loop bounds
            let kmer_val = unsafe { *self.window.get(i).unwrap_unchecked() };
            let score_val = if (kmer_val & 1) == 0 {
                let k = (kmer_val >> 1) as usize;
                // SAFETY: k < 2^(2*kmer) = ht.len()
                let htf_k = unsafe { *self.ht.get_unchecked(k) };
                let new_htf_k = htf_k + 1;
                unsafe {
                    *self.ht.get_unchecked_mut(k) = new_htf_k;
                }
                // SAFETY: new_htf_k <= window_size (bounded by queue length)
                // c is sized to window_size + 1
                unsafe { *self.c.get_unchecked(new_htf_k as usize) }
            } else {
                0.0
            };
            s += score_val - self.options.score_threshold;

            // SAFETY: l starts at 1, increments each iteration
            // l <= (q_len - i0) <= window_size, f is sized to window_size + 1
            let sl = s - unsafe { *self.f.get_unchecked(l) };
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

    // Quick backward-check heuristic: returns true if backward scan is worth trying
    fn if_backward(&self, max_step: i32) -> bool {
        let mut s = 0.0;
        for i in (0..self.window.len())
            .rev()
            .take((max_step as usize).min(self.window.len()))
        {
            // SAFETY: i is from the reverse iterator over 0..q.len(), so it's in bounds
            let kmer_val = unsafe { *self.window.get(i).unwrap_unchecked() };
            let val = if (kmer_val & 1) == 0 {
                let k = (kmer_val >> 1) as usize;
                // SAFETY: k < 2^(2*kmer) = window_ht.len()
                let wht_k = unsafe { *self.window_ht.get_unchecked(k) };
                // SAFETY: wht_k <= window_size, c is sized to window_size + 1
                unsafe { *self.c.get_unchecked(wht_k as usize) }
            } else {
                0.0
            };
            s += val - self.options.score_threshold;
            if s < 0.0 {
                return false;
            }
        }
        true
    }

    // Try to extend at the last position in queue. Returns 0 on success, -1 on fail
    fn extend(&mut self) -> i32 {
        if self.window.is_empty() {
            return -1;
        }
        let kmer_val = *self.window.back().unwrap();
        if (kmer_val & 1) != 0 {
            return -1;
        }
        let k = (kmer_val >> 1) as usize;
        let l = self.window.len() - 1;
        // SAFETY: k < 2^(2*kmer) = ht.len()
        let ht_k = unsafe { *self.ht.get_unchecked(k) };
        let idx = ht_k as usize + 1;
        if idx >= self.c.len() || (l + 1) >= self.f.len() {
            return -1;
        }
        // SAFETY: We just checked idx < c.len() and l+1 < f.len()
        if unsafe { *self.c.get_unchecked(idx) }
            - (unsafe { *self.f.get_unchecked(l + 1) } - unsafe { *self.f.get_unchecked(l) })
            < self.options.score_threshold
        {
            return -1;
        }
        // SAFETY: k < ht.len() as established above
        unsafe {
            *self.ht.get_unchecked_mut(k) = ht_k + 1;
        }
        0
    }

    // Merge two sorted lists of intervals (forward and reverse) into results
    fn merge_intervals(&mut self, fwd: Vec<Range<usize>>, rev: Vec<Range<usize>>) {
        self.results.clear();
        let mut i = 0;
        let mut j = 0;
        let mut start = 0;
        let mut end = 0;

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

            if intv.start <= end {
                end = end.max(intv.end);
            } else {
                if end > start {
                    self.results.push(start..end);
                }
                start = intv.start;
                end = intv.end;
            }
        }
        if end > start {
            self.results.push(start..end);
        }
    }

    // Save an interval into results keeping the sorted/merged invariant
    fn save_interval(&mut self, start: usize, end: usize) {
        let mut k = self.results.len();
        while k > 0 && start <= self.results[k - 1].end {
            k -= 1;
        }
        if k < self.results.len() {
            if start < self.results[k].start {
                self.results[k].start = start;
            }
            if end > self.results[k].end {
                self.results[k].end = end;
            }
            self.results.truncate(k + 1);
        } else {
            self.results.push(start..end);
        }
    }

    // Math helpers for f() table computation:

    // Stirling's approximation for log(n!)
    fn stirlings_approx(lambda: f64) -> f64 {
        lambda * (lambda.ln() - 1.0) + 0.5 * (2.0 * PI * E * lambda).ln()
    }

    // Calculates a high-precision approximation for the statistical score
    // adjustment term f[L] when the expected k-mer frequency (lambda) is large (>= 30.0).
    // It uses a specialized form of Stirling's approximation to ensure numerical stability
    // and high performance, avoiding computationally expensive summations and overflows.
    fn f_large(lambda: f64) -> f64 {
        Self::stirlings_approx(lambda)
            - 1.0 / (12.0 * lambda) * (1.0 + 0.5 / lambda + 19.0 / (30.0 * lambda * lambda))
    }

    // Calculates the default DUST score adjustment table (f[L]) for all lengths L <= max_l.
    // This version assumes a uniform background (50% GC content) as it passes a
    // density ratio (dr) of 1.0, making it an optimized shortcut for ld_cal_f2 when
    // GC correction is disabled.
    fn calculate_f(k: usize, max_l: usize) -> Vec<f64> {
        let n_kmer = 1i32 << (2 * k);
        let dr = 1.0;
        Self::calculate_f_internal(k, max_l, 1, &[n_kmer], &[dr])
    }

    // Calculates the DUST score adjustment table (f[L]) for all lengths L <= max_l.
    // This version uses k-mer density ratios (dr) based on genome-wide GC content
    // to adjust the expected scores, preventing over-masking of GC/AT-rich sequences.
    // It relies on a summation for small lambda and f_large() for large lambda.
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
        const MAX_N: usize = 10000;
        let n_kmer = 1usize << (2 * k);
        let mut f = vec![0.0f64; max_l + 1];
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
