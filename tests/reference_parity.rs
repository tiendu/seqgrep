use seqgrep::{
    AmbigMode, ExecutionConfig, SearchEngine, SearchQuery, SequenceRecord, SequenceType,
};

const GAP: u8 = 1 << 4;

fn next_u64(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    *state
}

fn generated_sequence(state: &mut u64, alphabet: &[u8], length: usize) -> String {
    (0..length)
        .map(|_| alphabet[(next_u64(state) as usize) % alphabet.len()] as char)
        .collect()
}

fn nucleotide_mask(symbol: u8) -> u8 {
    const A: u8 = 1;
    const C: u8 = 2;
    const G: u8 = 4;
    const T: u8 = 8;

    match symbol {
        b'A' => A,
        b'C' => C,
        b'G' => G,
        b'T' | b'U' => T,
        b'R' => A | G,
        b'Y' => C | T,
        b'S' => G | C,
        b'W' => A | T,
        b'K' => G | T,
        b'M' => A | C,
        b'B' => C | G | T,
        b'D' => A | G | T,
        b'H' => A | C | T,
        b'V' => A | C | G,
        b'N' => A | C | G | T,
        b'-' | b'.' => GAP,
        _ => panic!("test generated an unsupported nucleotide"),
    }
}

fn nucleotide_matches(query: u8, target: u8, mode: AmbigMode) -> bool {
    match mode {
        AmbigMode::None => {
            let query = if query == b'U' { b'T' } else { query };
            let target = if target == b'U' { b'T' } else { target };
            query == target
        }
        AmbigMode::Query => {
            let target_mask = match target {
                b'A' | b'C' | b'G' | b'T' | b'U' => nucleotide_mask(target),
                b'-' | b'.' => GAP,
                _ => 0,
            };
            nucleotide_mask(query) & target_mask != 0
        }
        AmbigMode::Both => nucleotide_mask(query) & nucleotide_mask(target) != 0,
        _ => unreachable!("test only uses published ambiguity modes"),
    }
}

fn brute_nucleotide_starts(
    target: &str,
    pattern: &str,
    mode: AmbigMode,
    circular: bool,
) -> Vec<usize> {
    let target = target.as_bytes();
    let pattern = pattern.as_bytes();
    if target.is_empty() || (!circular && pattern.len() > target.len()) {
        return Vec::new();
    }

    let starts = if circular {
        target.len()
    } else {
        target.len() - pattern.len() + 1
    };

    (0..starts)
        .filter(|&start| {
            pattern.iter().copied().enumerate().all(|(offset, query)| {
                let target_index = if circular {
                    (start + offset) % target.len()
                } else {
                    start + offset
                };
                nucleotide_matches(query, target[target_index], mode)
            })
        })
        .map(|start| start + 1)
        .collect()
}

#[test]
fn nucleotide_search_matches_independent_reference() {
    let alphabet = b"ACGTUNRY-.";
    let serial = SearchEngine::default();
    let parallel = SearchEngine::new(ExecutionConfig::new(4, 3).unwrap());
    let mut state = 0x5e_ed_5e_ed_d1_ce_ba_5e;

    for case_index in 0..2_000 {
        let target_len = (next_u64(&mut state) % 10 + 1) as usize;
        let pattern_len = (next_u64(&mut state) % 13 + 1) as usize;
        let target = generated_sequence(&mut state, alphabet, target_len);
        let pattern = generated_sequence(&mut state, alphabet, pattern_len);
        let circular = next_u64(&mut state) & 1 == 1;
        let mode = match next_u64(&mut state) % 3 {
            0 => AmbigMode::None,
            1 => AmbigMode::Query,
            _ => AmbigMode::Both,
        };

        let record = SequenceRecord::new("generated", &target);
        let query = SearchQuery::new(&pattern)
            .with_ambig_mode(mode)
            .with_circular(circular);
        let expected = brute_nucleotide_starts(&target, &pattern, mode, circular);
        let actual = serial
            .search_record(&record, &query)
            .unwrap()
            .into_iter()
            .map(|matched| matched.start)
            .collect::<Vec<_>>();

        assert_eq!(
            actual, expected,
            "case {case_index}: target={target:?}, pattern={pattern:?}, mode={mode:?}, circular={circular}"
        );

        if case_index % 20 == 0 {
            assert_eq!(
                parallel.search_record(&record, &query).unwrap(),
                serial.search_record(&record, &query).unwrap(),
                "parallel mismatch in case {case_index}"
            );
        }
    }
}

#[test]
fn protein_search_matches_independent_reference() {
    let alphabet = b"ACDEFGHIKLMNPQRSTVWYBJZXUO*-.";
    let serial = SearchEngine::default();
    let parallel = SearchEngine::new(ExecutionConfig::new(3, 2).unwrap());
    let mut state = 0xc0_ff_ee_12_34_56_78_90;

    for case_index in 0..500 {
        let target_len = (next_u64(&mut state) % 12 + 1) as usize;
        let pattern_len = (next_u64(&mut state) % 15 + 1) as usize;
        let target = generated_sequence(&mut state, alphabet, target_len);
        let pattern = generated_sequence(&mut state, alphabet, pattern_len);
        let circular = next_u64(&mut state) & 1 == 1;

        let expected = if !circular && pattern.len() > target.len() {
            Vec::new()
        } else {
            let starts = if circular {
                target.len()
            } else {
                target.len() - pattern.len() + 1
            };
            (0..starts)
                .filter(|&start| {
                    pattern.bytes().enumerate().all(|(offset, query)| {
                        let index = if circular {
                            (start + offset) % target.len()
                        } else {
                            start + offset
                        };
                        query == target.as_bytes()[index]
                    })
                })
                .map(|start| start + 1)
                .collect::<Vec<_>>()
        };

        let record = SequenceRecord::new("protein", &target);
        let query = SearchQuery::new(&pattern)
            .with_sequence_type(SequenceType::AminoAcid)
            .with_circular(circular);
        let actual = serial
            .search_record(&record, &query)
            .unwrap()
            .into_iter()
            .map(|matched| matched.start)
            .collect::<Vec<_>>();

        assert_eq!(
            actual, expected,
            "case {case_index}: target={target:?}, pattern={pattern:?}, circular={circular}"
        );

        if case_index % 10 == 0 {
            assert_eq!(
                parallel.search_record(&record, &query).unwrap(),
                serial.search_record(&record, &query).unwrap(),
                "parallel mismatch in case {case_index}"
            );
        }
    }
}
