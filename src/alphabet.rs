use crate::error::{Error, Result};

pub(crate) const GAP_MASK: u8 = 1 << 4;
pub(crate) const PROTEIN_SYMBOLS: &str = "ACDEFGHIKLMNPQRSTVWYBJZXUO*-.";

pub(crate) fn normalize_symbols(sequence: &str) -> String {
    sequence
        .chars()
        .filter(|symbol| !symbol.is_whitespace())
        .map(|symbol| symbol.to_ascii_uppercase())
        .collect()
}

fn invalid_symbol(alphabet: &'static str, symbol: char, position: usize) -> Error {
    Error::InvalidSymbol {
        alphabet,
        symbol,
        position,
    }
}

fn encode_with(
    sequence: &str,
    alphabet: &'static str,
    code: fn(char) -> Option<u8>,
) -> Result<Vec<u8>> {
    let mut encoded = Vec::with_capacity(sequence.len());
    for (index, symbol) in sequence.chars().enumerate() {
        encoded.push(code(symbol).ok_or_else(|| invalid_symbol(alphabet, symbol, index + 1))?);
    }
    Ok(encoded)
}

pub(crate) fn nucleotide_2bit_code(symbol: char) -> Option<u8> {
    match symbol {
        'A' => Some(0b00),
        'C' => Some(0b01),
        'G' => Some(0b10),
        'T' | 'U' => Some(0b11),
        _ => None,
    }
}

pub(crate) fn nucleotide_exact_code(symbol: char) -> Option<u8> {
    match symbol {
        'A' => Some(0),
        'C' => Some(1),
        'G' => Some(2),
        'T' | 'U' => Some(3),
        'R' => Some(4),
        'Y' => Some(5),
        'S' => Some(6),
        'W' => Some(7),
        'K' => Some(8),
        'M' => Some(9),
        'B' => Some(10),
        'D' => Some(11),
        'H' => Some(12),
        'V' => Some(13),
        'N' => Some(14),
        '-' => Some(15),
        '.' => Some(16),
        _ => None,
    }
}

pub(crate) fn nucleotide_iupac_query_mask(symbol: char) -> Option<u8> {
    const A: u8 = 1 << 0;
    const C: u8 = 1 << 1;
    const G: u8 = 1 << 2;
    const T: u8 = 1 << 3;

    match symbol {
        'A' => Some(A),
        'C' => Some(C),
        'G' => Some(G),
        'T' | 'U' => Some(T),
        'R' => Some(A | G),
        'Y' => Some(C | T),
        'S' => Some(G | C),
        'W' => Some(A | T),
        'K' => Some(G | T),
        'M' => Some(A | C),
        'B' => Some(C | G | T),
        'D' => Some(A | G | T),
        'H' => Some(A | C | T),
        'V' => Some(A | C | G),
        'N' => Some(A | C | G | T),
        '-' | '.' => Some(GAP_MASK),
        _ => None,
    }
}

pub(crate) fn nucleotide_iupac_target_mask(symbol: char) -> Option<u8> {
    match symbol {
        'A' | 'C' | 'G' | 'T' | 'U' => nucleotide_iupac_query_mask(symbol),
        'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'B' | 'D' | 'H' | 'V' | 'N' => Some(0),
        '-' | '.' => Some(GAP_MASK),
        _ => None,
    }
}

pub(crate) fn protein_code(symbol: char) -> Option<u8> {
    if !symbol.is_ascii() {
        return None;
    }
    PROTEIN_SYMBOLS
        .bytes()
        .position(|candidate| candidate == symbol as u8)
        .map(|index| index as u8)
}

pub(crate) fn encode_nucleotide_exact(sequence: &str) -> Result<Vec<u8>> {
    encode_with(sequence, "nucleotide", nucleotide_exact_code)
}

pub(crate) fn encode_nucleotide_iupac_query(sequence: &str) -> Result<Vec<u8>> {
    encode_with(
        sequence,
        "IUPAC nucleotide query",
        nucleotide_iupac_query_mask,
    )
}

pub(crate) fn encode_protein_exact(sequence: &str) -> Result<Vec<u8>> {
    encode_with(sequence, "amino-acid", protein_code)
}

pub(crate) fn nucleotide_comparison_bytes(sequence: &str) -> Result<Vec<u8>> {
    encode_nucleotide_exact(sequence)
}

pub(crate) fn protein_comparison_bytes(sequence: &str) -> Result<Vec<u8>> {
    encode_protein_exact(sequence)
}

pub(crate) fn reverse_complement(sequence: &str) -> Result<String> {
    let mut result = String::with_capacity(sequence.len());
    let symbol_count = sequence.chars().count();

    for (reverse_index, symbol) in sequence.chars().rev().enumerate() {
        let complement = match symbol {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' | 'U' => 'A',
            'R' => 'Y',
            'Y' => 'R',
            'S' => 'S',
            'W' => 'W',
            'K' => 'M',
            'M' => 'K',
            'B' => 'V',
            'D' => 'H',
            'H' => 'D',
            'V' => 'B',
            'N' => 'N',
            '-' => '-',
            '.' => '.',
            _ => {
                let position = symbol_count - reverse_index;
                return Err(invalid_symbol("nucleotide", symbol, position));
            }
        };
        result.push(complement);
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalizes_case_and_whitespace() {
        assert_eq!(normalize_symbols("a c\ng\tt"), "ACGT");
    }

    #[test]
    fn exact_nucleotide_codes_t_and_u_equally() {
        assert_eq!(encode_nucleotide_exact("TU").unwrap(), vec![3, 3]);
    }

    #[test]
    fn target_ambiguity_is_asymmetric() {
        assert_eq!(nucleotide_iupac_target_mask('N'), Some(0));
        assert_eq!(nucleotide_iupac_query_mask('N'), Some(0b1111));
    }

    #[test]
    fn reverse_complement_handles_iupac_and_gaps() {
        assert_eq!(reverse_complement("ARYN-.").unwrap(), ".-NRYT");
    }

    #[test]
    fn protein_alphabet_fits_five_bits() {
        assert!(PROTEIN_SYMBOLS.len() <= 32);
    }
}
