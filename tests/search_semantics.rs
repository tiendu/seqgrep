use seqgrep::{
    AmbigMode, ExecutionConfig, Match, SearchEngine, SearchQuery, SequenceRecord, SequenceType,
    Strand,
};

fn search(record: &str, pattern: &str) -> Vec<Match> {
    SearchEngine::default()
        .search_record(
            &SequenceRecord::new("seq", record),
            &SearchQuery::new(pattern),
        )
        .unwrap()
}

fn values(matches: &[Match]) -> Vec<(usize, usize, &str)> {
    matches
        .iter()
        .map(|item| (item.start, item.end, item.matched.as_str()))
        .collect()
}

#[test]
fn default_mode_is_exact_nucleotide() {
    assert_eq!(values(&search("ATGNATGA", "ATGN")), vec![(1, 4, "ATGN")]);
}

#[test]
fn exact_nucleotide_treats_t_and_u_equally() {
    assert_eq!(
        values(&search("AUGATG", "ATG")),
        vec![(1, 3, "AUG"), (4, 6, "ATG")]
    );
}

#[test]
fn exact_gaps_are_literal() {
    assert_eq!(values(&search("A-C.A.C", "A-C")), vec![(1, 3, "A-C")]);
    assert_eq!(values(&search("A-C.A.C", "A.C")), vec![(5, 7, "A.C")]);
}

#[test]
fn iupac_gaps_are_equivalent() {
    let query = SearchQuery::new("A-C").with_ambig_mode(AmbigMode::Query);
    let matches = SearchEngine::default()
        .search_record(&SequenceRecord::new("aligned", "A-CA.C"), &query)
        .unwrap();
    assert_eq!(values(&matches), vec![(1, 3, "A-C"), (4, 6, "A.C")]);
}

#[test]
fn query_n_matches_canonical_bases() {
    let query = SearchQuery::new("N").with_ambig_mode(AmbigMode::Query);
    let matches = SearchEngine::default()
        .search_record(&SequenceRecord::new("seq", "ACGT"), &query)
        .unwrap();
    assert_eq!(
        values(&matches),
        vec![(1, 1, "A"), (2, 2, "C"), (3, 3, "G"), (4, 4, "T")]
    );
}

#[test]
fn query_mode_does_not_treat_target_n_as_wildcard() {
    let query = SearchQuery::new("ATGN").with_ambig_mode(AmbigMode::Query);
    let matches = SearchEngine::default()
        .search_record(&SequenceRecord::new("seq", "ATGNATGA"), &query)
        .unwrap();
    assert_eq!(values(&matches), vec![(5, 8, "ATGA")]);
}

#[test]
fn both_mode_matches_ambiguous_target() {
    let query = SearchQuery::new("AAAAA").with_ambig_mode(AmbigMode::Both);
    let matches = SearchEngine::default()
        .search_record(&SequenceRecord::new("uncertain", "NNANN"), &query)
        .unwrap();
    assert_eq!(values(&matches), vec![(1, 5, "NNANN")]);
}

#[test]
fn exact_mode_matches_n_literally() {
    assert_eq!(
        values(&search("NNNN", "NN")),
        vec![(1, 2, "NN"), (2, 3, "NN"), (3, 4, "NN")]
    );
}

#[test]
fn protein_mode_is_exact() {
    let query = SearchQuery::new("MTEYK").with_sequence_type(SequenceType::AminoAcid);
    let matches = SearchEngine::default()
        .search_record(
            &SequenceRecord::new("ras", "XXMTEYKLVVVGAGGVGKSALXX"),
            &query,
        )
        .unwrap();
    assert_eq!(values(&matches), vec![(3, 7, "MTEYK")]);
}

#[test]
fn protein_rejects_nucleotide_only_options() {
    let query = SearchQuery::new("MTE")
        .with_sequence_type(SequenceType::AminoAcid)
        .with_ambig_mode(AmbigMode::Query);
    assert!(SearchEngine::default()
        .search_record(&SequenceRecord::new("protein", "MTE"), &query)
        .unwrap_err()
        .to_string()
        .contains("nucleotide"));
}

#[test]
fn reverse_complement_match() {
    let query = SearchQuery::new("ACTTCAT").with_reverse_complement(true);
    let matches = SearchEngine::default()
        .search_record(&SequenceRecord::new("seq", "CCCATGAAGTCCC"), &query)
        .unwrap();
    assert!(matches.iter().any(|item| {
        item.strand == Strand::Reverse
            && item.start == 4
            && item.end == 10
            && item.matched == "ATGAAGT"
    }));
}

#[test]
fn circular_search_wraps_coordinates() {
    let query = SearchQuery::new("TTTCCC").with_circular(true);
    let matches = SearchEngine::default()
        .search_record(&SequenceRecord::new("plasmid", "CCCAAATTT"), &query)
        .unwrap();
    assert_eq!(matches.len(), 1);
    assert_eq!(matches[0].record, "plasmid");
    assert_eq!(matches[0].strand, Strand::Forward);
    assert_eq!(matches[0].start, 7);
    assert_eq!(matches[0].end, 3);
    assert_eq!(matches[0].matched, "TTTCCC");
    assert!(matches[0].circular);
}

#[test]
fn circular_pattern_can_be_longer_than_record() {
    let query = SearchQuery::new("ATGAT").with_circular(true);
    let matches = SearchEngine::default()
        .search_record(&SequenceRecord::new("tiny", "ATG"), &query)
        .unwrap();
    assert_eq!(values(&matches), vec![(1, 2, "ATGAT")]);
}

#[test]
fn parallel_matches_serial_for_all_modes() {
    let cases = [
        (
            SequenceRecord::new("nucleotide", "AAAAATGCAAAAA"),
            SearchQuery::new("TGCA"),
        ),
        (
            SequenceRecord::new("mixed", "A-CATGNA.C"),
            SearchQuery::new("A.C").with_ambig_mode(AmbigMode::Query),
        ),
        (
            SequenceRecord::new("uncertain", "NNANN"),
            SearchQuery::new("AAAAA").with_ambig_mode(AmbigMode::Both),
        ),
        (
            SequenceRecord::new("protein", "AAAAAMTEYKAAAAA"),
            SearchQuery::new("MTEYK").with_sequence_type(SequenceType::AminoAcid),
        ),
    ];

    let serial = SearchEngine::default();
    let parallel = SearchEngine::new(ExecutionConfig::new(3, 2).unwrap());

    for (record, query) in cases {
        assert_eq!(
            parallel.search_record(&record, &query).unwrap(),
            serial.search_record(&record, &query).unwrap()
        );
    }
}

#[test]
fn models_normalize_at_boundary() {
    assert_eq!(SequenceRecord::new("seq", "a c\ng").sequence(), "ACG");
    assert_eq!(SearchQuery::new("a t\ng").pattern(), "ATG");
}

#[test]
fn streaming_visitor_preserves_order() {
    let record = SequenceRecord::new("seq", "AAAA");
    let query = SearchQuery::new("AA");
    let engine = SearchEngine::new(ExecutionConfig::new(3, 1).unwrap());
    let mut starts = Vec::new();

    engine
        .visit_matches(&record, &query, |matched| {
            starts.push(matched.start);
            Ok(())
        })
        .unwrap();

    assert_eq!(starts, vec![1, 2, 3]);
}

#[test]
fn streaming_visitor_can_stop_parallel_search_with_an_error() {
    let record = SequenceRecord::new("seq", "A".repeat(10_000));
    let query = SearchQuery::new("A");
    let engine = SearchEngine::new(ExecutionConfig::new(4, 64).unwrap());
    let mut seen = 0;

    let error = engine
        .visit_matches(&record, &query, |_| {
            seen += 1;
            if seen == 3 {
                return Err(seqgrep::Error::InvalidOption("stop".to_owned()));
            }
            Ok(())
        })
        .unwrap_err();

    assert_eq!(seen, 3);
    assert_eq!(error.to_string(), "stop");
}

#[test]
fn symmetric_ambiguity_matches_n_to_n() {
    let query = SearchQuery::new("NNNNN").with_ambig_mode(AmbigMode::Both);
    let matches = SearchEngine::default()
        .search_record(&SequenceRecord::new("uncertain", "NNNNN"), &query)
        .unwrap();
    assert_eq!(values(&matches), vec![(1, 5, "NNNNN")]);
}

#[test]
fn protein_extended_symbols_remain_literal() {
    let query = SearchQuery::new("MAX").with_sequence_type(SequenceType::AminoAcid);
    let matches = SearchEngine::default()
        .search_record(&SequenceRecord::new("protein", "MAXMAB"), &query)
        .unwrap();
    assert_eq!(values(&matches), vec![(1, 3, "MAX")]);
}

#[test]
fn circular_protein_search_wraps() {
    let query = SearchQuery::new("MAC")
        .with_sequence_type(SequenceType::AminoAcid)
        .with_circular(true);
    let matches = SearchEngine::default()
        .search_record(&SequenceRecord::new("ring", "CDEMA"), &query)
        .unwrap();
    assert_eq!(values(&matches), vec![(4, 1, "MAC")]);
}

#[test]
fn parallel_reverse_complement_circular_search_matches_serial() {
    let record = SequenceRecord::new("plasmid", "CCCAAATTT");
    let query = SearchQuery::new("GGGAAA")
        .with_ambig_mode(AmbigMode::Query)
        .with_reverse_complement(true)
        .with_circular(true);
    let serial = SearchEngine::default()
        .search_record(&record, &query)
        .unwrap();
    let parallel = SearchEngine::new(ExecutionConfig::new(3, 2).unwrap())
        .search_record(&record, &query)
        .unwrap();
    assert_eq!(parallel, serial);
}

#[test]
fn prepared_query_is_reusable_across_records() {
    let engine = SearchEngine::default();
    let prepared = engine.prepare_query(SearchQuery::new("ATG")).unwrap();

    assert_eq!(
        values(
            &engine
                .search_prepared(&SequenceRecord::new("one", "ATGCCC"), &prepared)
                .unwrap()
        ),
        vec![(1, 3, "ATG")]
    );
    assert_eq!(
        values(
            &engine
                .search_prepared(&SequenceRecord::new("two", "CCCATG"), &prepared)
                .unwrap()
        ),
        vec![(4, 6, "ATG")]
    );
}

#[test]
fn empty_patterns_and_invalid_symbols_are_rejected() {
    let engine = SearchEngine::default();
    let record = SequenceRecord::new("seq", "ATGC");

    assert!(engine
        .search_record(&record, &SearchQuery::new(""))
        .unwrap_err()
        .to_string()
        .contains("must not be empty"));
    assert!(engine
        .search_record(&record, &SearchQuery::new("ATX"))
        .unwrap_err()
        .to_string()
        .contains("Unsupported nucleotide symbol"));
    assert!(engine
        .search_record(&SequenceRecord::new("bad", "ATX"), &SearchQuery::new("ATG"),)
        .unwrap_err()
        .to_string()
        .contains("Unsupported nucleotide symbol"));
}

#[test]
fn empty_records_have_no_matches() {
    let matches = SearchEngine::default()
        .search_record(&SequenceRecord::new("empty", ""), &SearchQuery::new("A"))
        .unwrap();
    assert!(matches.is_empty());
}
