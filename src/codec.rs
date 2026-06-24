use crate::alphabet::{
    encode_nucleotide_exact, encode_nucleotide_iupac_query, encode_protein_exact,
    nucleotide_comparison_bytes, nucleotide_exact_code, nucleotide_iupac_query_mask, protein_code,
    protein_comparison_bytes, reverse_complement,
};
use crate::error::Result;
use crate::packed::{Packed2, Packed5, QueryTarget};

#[derive(Debug, Clone)]
pub(crate) enum EncodedTarget {
    Nucleotide2Code(Packed2),
    NucleotideExact5(Packed5),
    Nucleotide2Mask(Packed2),
    NucleotideQueryTarget(QueryTarget),
    NucleotideIupac5(Packed5),
    Protein5(Packed5),
}

impl EncodedTarget {
    pub(crate) fn len(&self) -> usize {
        match self {
            Self::Nucleotide2Code(target) | Self::Nucleotide2Mask(target) => target.len(),
            Self::NucleotideExact5(target)
            | Self::NucleotideIupac5(target)
            | Self::Protein5(target) => target.len(),
            Self::NucleotideQueryTarget(target) => target.len(),
        }
    }

    pub(crate) fn symbol_at(&self, index: usize) -> u8 {
        debug_assert!(index < self.len());
        match self {
            Self::Nucleotide2Code(target) => target.get(index),
            Self::NucleotideExact5(target) | Self::Protein5(target) => target.get(index),
            Self::Nucleotide2Mask(target) => 1 << target.get(index),
            Self::NucleotideQueryTarget(target) => target.get(index),
            Self::NucleotideIupac5(target) => target.get(index),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Codec {
    NucleotideExact,
    NucleotideIupac { target_ambiguity: bool },
    ProteinExact,
}

impl Codec {
    pub(crate) fn encode_query(self, sequence: &str) -> Result<Vec<u8>> {
        match self {
            Self::NucleotideExact => encode_nucleotide_exact(sequence),
            Self::NucleotideIupac { .. } => encode_nucleotide_iupac_query(sequence),
            Self::ProteinExact => encode_protein_exact(sequence),
        }
    }

    pub(crate) fn encode_target(self, sequence: &str) -> Result<EncodedTarget> {
        match self {
            Self::NucleotideExact => match Packed2::encode(sequence) {
                Ok(target) => Ok(EncodedTarget::Nucleotide2Code(target)),
                Err(_) => {
                    let codes = sequence
                        .chars()
                        .enumerate()
                        .map(|(index, symbol)| {
                            nucleotide_exact_code(symbol).ok_or_else(|| {
                                crate::Error::InvalidSymbol {
                                    alphabet: "nucleotide",
                                    symbol,
                                    position: index + 1,
                                }
                            })
                        })
                        .collect::<Result<Vec<_>>>()?;
                    Ok(EncodedTarget::NucleotideExact5(Packed5::from_codes(&codes)))
                }
            },
            Self::NucleotideIupac {
                target_ambiguity: false,
            } => match Packed2::encode(sequence) {
                Ok(target) => Ok(EncodedTarget::Nucleotide2Mask(target)),
                Err(_) => Ok(EncodedTarget::NucleotideQueryTarget(QueryTarget::encode(
                    sequence,
                )?)),
            },
            Self::NucleotideIupac {
                target_ambiguity: true,
            } => match Packed2::encode(sequence) {
                Ok(target) => Ok(EncodedTarget::Nucleotide2Mask(target)),
                Err(_) => {
                    let masks = sequence
                        .chars()
                        .enumerate()
                        .map(|(index, symbol)| {
                            nucleotide_iupac_query_mask(symbol).ok_or_else(|| {
                                crate::Error::InvalidSymbol {
                                    alphabet: "IUPAC nucleotide",
                                    symbol,
                                    position: index + 1,
                                }
                            })
                        })
                        .collect::<Result<Vec<_>>>()?;
                    Ok(EncodedTarget::NucleotideIupac5(Packed5::from_codes(&masks)))
                }
            },
            Self::ProteinExact => {
                let codes = sequence
                    .chars()
                    .enumerate()
                    .map(|(index, symbol)| {
                        protein_code(symbol).ok_or_else(|| crate::Error::InvalidSymbol {
                            alphabet: "amino-acid",
                            symbol,
                            position: index + 1,
                        })
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(EncodedTarget::Protein5(Packed5::from_codes(&codes)))
            }
        }
    }

    pub(crate) fn comparison_bytes(self, sequence: &str) -> Result<Vec<u8>> {
        match self {
            Self::NucleotideExact => nucleotide_comparison_bytes(sequence),
            Self::ProteinExact => protein_comparison_bytes(sequence),
            Self::NucleotideIupac { .. } => unreachable!("IUPAC search is not exact"),
        }
    }

    pub(crate) fn compatible(self, query_symbol: u8, target_symbol: u8) -> bool {
        match self {
            Self::NucleotideExact | Self::ProteinExact => query_symbol == target_symbol,
            Self::NucleotideIupac { .. } => query_symbol & target_symbol != 0,
        }
    }

    pub(crate) fn reverse_complement(self, pattern: &str) -> Result<String> {
        match self {
            Self::NucleotideExact | Self::NucleotideIupac { .. } => reverse_complement(pattern),
            Self::ProteinExact => Err(crate::Error::InvalidOption(
                "reverse-complement search is not valid for amino acids".to_owned(),
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_nucleotide_uses_two_bits_when_canonical() {
        let target = Codec::NucleotideExact.encode_target("ACGTU").unwrap();
        assert!(matches!(target, EncodedTarget::Nucleotide2Code(_)));
        assert_eq!(
            (0..target.len())
                .map(|index| target.symbol_at(index))
                .collect::<Vec<_>>(),
            vec![0, 1, 2, 3, 3]
        );
    }

    #[test]
    fn exact_nucleotide_uses_five_bits_when_mixed() {
        let target = Codec::NucleotideExact.encode_target("ANRY-.").unwrap();
        assert!(matches!(target, EncodedTarget::NucleotideExact5(_)));
    }

    #[test]
    fn query_mode_unknown_targets_decode_to_zero() {
        let target = Codec::NucleotideIupac {
            target_ambiguity: false,
        }
        .encode_target("ANRY")
        .unwrap();
        assert_eq!(
            (0..target.len())
                .map(|index| target.symbol_at(index))
                .collect::<Vec<_>>(),
            vec![1, 0, 0, 0]
        );
    }

    #[test]
    fn both_mode_preserves_target_masks() {
        let target = Codec::NucleotideIupac {
            target_ambiguity: true,
        }
        .encode_target("ANRY")
        .unwrap();
        assert_eq!(
            (0..target.len())
                .map(|index| target.symbol_at(index))
                .collect::<Vec<_>>(),
            vec![1, 15, 5, 10]
        );
    }
}
