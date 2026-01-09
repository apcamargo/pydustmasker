/// Lookup to encode ASCII DNA letters into 0..4
/// A -> 0, C -> 1, G -> 2, T -> 3, others -> 4
pub const ENCODING_LOOKUP: [u8; 256] = {
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

/// Encode ASCII DNA sequence into numeric values (A→0, C→1, G→2, T→3, others→4)
pub fn encode_sequence(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .map(|&b| ENCODING_LOOKUP[b as usize])
        .collect()
}

/// Reverse complement for encoded sequence:
/// 0(A) <-> 3(T), 1(C) <-> 2(G), 4(N) -> 4(N)
pub fn reverse_complement_encoded_sequence(encoded_seq: &[u8]) -> Vec<u8> {
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

// Compute GC content from an encoded sequence.
// Returns a value between 0.0 and 1.0. Returns 0.0 if the sequence is empty.
pub fn compute_gc_content(encoded_seq: &[u8]) -> f64 {
    if encoded_seq.is_empty() {
        return 0.0;
    }
    let gc_count = encoded_seq.iter().filter(|&&b| b == 1 || b == 2).count();
    gc_count as f64 / encoded_seq.len() as f64
}
