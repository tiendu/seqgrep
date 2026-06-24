use crate::alphabet::{nucleotide_2bit_code, nucleotide_iupac_target_mask, GAP_MASK};
use crate::error::{Error, Result};

#[derive(Debug, Clone)]
pub(crate) struct Packed2 {
    data: Vec<u8>,
    len: usize,
}

impl Packed2 {
    pub(crate) fn encode(sequence: &str) -> Result<Self> {
        let len = sequence.chars().count();
        let mut data = vec![0_u8; len.div_ceil(4)];

        for (index, symbol) in sequence.chars().enumerate() {
            let code = nucleotide_2bit_code(symbol).ok_or_else(|| Error::InvalidSymbol {
                alphabet: "canonical nucleotide",
                symbol,
                position: index + 1,
            })?;
            let shift = 6 - (index % 4) * 2;
            data[index / 4] |= code << shift;
        }

        Ok(Self { data, len })
    }

    pub(crate) fn from_packed_bytes(data: Vec<u8>, len: usize) -> Self {
        Self { data, len }
    }

    pub(crate) fn get(&self, index: usize) -> u8 {
        debug_assert!(index < self.len);
        let shift = 6 - (index % 4) * 2;
        (self.data[index / 4] >> shift) & 0b11
    }

    pub(crate) fn len(&self) -> usize {
        self.len
    }

    #[cfg(test)]
    pub(crate) fn storage_len(&self) -> usize {
        self.data.len()
    }
}

#[derive(Debug, Clone)]
pub(crate) struct Packed5 {
    data: Vec<u8>,
    len: usize,
}

impl Packed5 {
    pub(crate) fn from_codes(codes: &[u8]) -> Self {
        let len = codes.len();
        let mut data = vec![0_u8; (len * 5).div_ceil(8)];

        for (index, code) in codes.iter().copied().enumerate() {
            debug_assert!(code <= 0b1_1111);
            let bit_offset = index * 5;
            let byte_index = bit_offset / 8;
            let offset_in_byte = bit_offset % 8;
            let first_bits = 5.min(8 - offset_in_byte);
            let remaining_bits = 5 - first_bits;

            if remaining_bits == 0 {
                data[byte_index] |= code << (8 - offset_in_byte - 5);
            } else {
                data[byte_index] |= code >> remaining_bits;
                let remainder_mask = (1 << remaining_bits) - 1;
                data[byte_index + 1] |= (code & remainder_mask) << (8 - remaining_bits);
            }
        }

        Self { data, len }
    }

    pub(crate) fn get(&self, index: usize) -> u8 {
        debug_assert!(index < self.len);
        let bit_offset = index * 5;
        let byte_index = bit_offset / 8;
        let offset_in_byte = bit_offset % 8;

        let mut value = u16::from(self.data[byte_index]) << 8;
        if byte_index + 1 < self.data.len() {
            value |= u16::from(self.data[byte_index + 1]);
        }

        let shift = 16 - offset_in_byte - 5;
        ((value >> shift) & 0b1_1111) as u8
    }

    pub(crate) fn len(&self) -> usize {
        self.len
    }

    #[cfg(test)]
    pub(crate) fn storage_len(&self) -> usize {
        self.data.len()
    }
}

#[derive(Debug, Clone)]
pub(crate) struct BitSet {
    data: Vec<u8>,
    len: usize,
}

impl BitSet {
    pub(crate) fn new(len: usize) -> Self {
        Self {
            data: vec![0; len.div_ceil(8)],
            len,
        }
    }

    pub(crate) fn set(&mut self, index: usize) {
        debug_assert!(index < self.len);
        self.data[index / 8] |= 1 << (7 - index % 8);
    }

    pub(crate) fn contains(&self, index: usize) -> bool {
        debug_assert!(index < self.len);
        self.data[index / 8] & (1 << (7 - index % 8)) != 0
    }
}

#[derive(Debug, Clone)]
pub(crate) struct QueryTarget {
    codes: Packed2,
    valid: BitSet,
    gaps: Option<BitSet>,
}

impl QueryTarget {
    pub(crate) fn encode(sequence: &str) -> Result<Self> {
        let len = sequence.chars().count();
        let mut code_bytes = vec![0_u8; len.div_ceil(4)];
        let mut valid = BitSet::new(len);
        let mut gaps: Option<BitSet> = None;

        for (index, symbol) in sequence.chars().enumerate() {
            let mask = nucleotide_iupac_target_mask(symbol).ok_or(Error::InvalidSymbol {
                alphabet: "IUPAC nucleotide target",
                symbol,
                position: index + 1,
            })?;

            match mask {
                0 => {}
                GAP_MASK => {
                    gaps.get_or_insert_with(|| BitSet::new(len)).set(index);
                }
                singleton if singleton.is_power_of_two() => {
                    let code = singleton.trailing_zeros() as u8;
                    let shift = 6 - (index % 4) * 2;
                    code_bytes[index / 4] |= code << shift;
                    valid.set(index);
                }
                _ => unreachable!("target masks are canonical, unknown, or gap"),
            }
        }

        Ok(Self {
            codes: Packed2::from_packed_bytes(code_bytes, len),
            valid,
            gaps,
        })
    }

    pub(crate) fn get(&self, index: usize) -> u8 {
        if self.valid.contains(index) {
            return 1 << self.codes.get(index);
        }

        if self.gaps.as_ref().is_some_and(|gaps| gaps.contains(index)) {
            return GAP_MASK;
        }

        0
    }

    pub(crate) fn len(&self) -> usize {
        self.codes.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_bit_round_trip() {
        let packed = Packed2::encode("ACGTUACGT").unwrap();
        assert_eq!(packed.storage_len(), 3);
        assert_eq!(
            (0..packed.len())
                .map(|index| packed.get(index))
                .collect::<Vec<_>>(),
            vec![0, 1, 2, 3, 3, 0, 1, 2, 3]
        );
    }

    #[test]
    fn five_bit_round_trip_across_boundaries() {
        for len in 1_usize..65 {
            let codes = (0..len).map(|index| (index % 29) as u8).collect::<Vec<_>>();
            let packed = Packed5::from_codes(&codes);
            assert_eq!(packed.storage_len(), (len * 5).div_ceil(8));
            assert_eq!(
                (0..len).map(|index| packed.get(index)).collect::<Vec<_>>(),
                codes
            );
        }
    }

    #[test]
    fn query_target_marks_unknowns_and_gaps() {
        let target = QueryTarget::encode("ANRY-.").unwrap();
        assert_eq!(
            (0..target.len())
                .map(|index| target.get(index))
                .collect::<Vec<_>>(),
            vec![1, 0, 0, 0, GAP_MASK, GAP_MASK]
        );
    }
}
