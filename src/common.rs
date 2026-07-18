/// The alphabet a sequence is encoded and scored against.
#[derive(Debug, Clone, Copy)]
pub enum Alphabet {
    Dna,
    Protein,
}

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

/// Lookup to encode ASCII protein letters, in the order
/// ACDEFGHIKLMNPQRSTVWY -> 0..20, everything else -> 20 (ambiguous)
const PROTEIN_ENCODING_LOOKUP: [u8; 256] = {
    let mut lookup = [20; 256];
    let letters = b"ACDEFGHIKLMNPQRSTVWY";
    let mut i = 0;
    while i < letters.len() {
        lookup[letters[i] as usize] = i as u8;
        lookup[letters[i].to_ascii_lowercase() as usize] = i as u8;
        i += 1;
    }
    lookup
};

/// Encode ASCII protein sequence into numeric values using the protein
/// alphabet (ACDEFGHIKLMNPQRSTVWY -> 0..20, others -> 20).
fn encode_protein_sequence(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .map(|&b| PROTEIN_ENCODING_LOOKUP[b as usize])
        .collect()
}

/// Encode an ASCII sequence using the given alphabet's lookup table.
pub fn encode_with_alphabet(sequence: &[u8], alphabet: Alphabet) -> Vec<u8> {
    match alphabet {
        Alphabet::Dna => encode_sequence(sequence),
        Alphabet::Protein => encode_protein_sequence(sequence),
    }
}

/// Inverse tables for decoding a consensus repeat unit. The trailing entry is
/// the ambiguous bucket ('N' for DNA, 'X' for protein) both encoders fall back
/// to, since a repeat unit may span an ambiguous letter.
const DNA_DECODING: [u8; 5] = *b"ACGTN";
const PROTEIN_DECODING: [u8; 21] = *b"ACDEFGHIKLMNPQRSTVWYX";

/// Decode an encoded consensus unit back into letters. Every value the two
/// encoders can produce (0..=4 for DNA, 0..=20 for protein) has an entry, so
/// the result is always ASCII.
pub fn decode_sequence(encoded: &[u8], alphabet: Alphabet) -> String {
    let table: &[u8] = match alphabet {
        Alphabet::Dna => &DNA_DECODING,
        Alphabet::Protein => &PROTEIN_DECODING,
    };
    encoded.iter().map(|&b| table[b as usize] as char).collect()
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
