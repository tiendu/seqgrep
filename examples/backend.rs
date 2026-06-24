use seqgrep::{AmbigMode, ExecutionConfig, SearchEngine, SearchQuery, SequenceRecord};

fn main() -> seqgrep::Result<()> {
    let engine = SearchEngine::new(ExecutionConfig::new(4, 1_000_000)?);
    let prepared =
        engine.prepare_query(SearchQuery::new("ATGN").with_ambig_mode(AmbigMode::Query))?;

    for record in [
        SequenceRecord::new("one", "ATGNATGA"),
        SequenceRecord::new("two", "CCCATGTCCC"),
    ] {
        engine.visit_prepared_matches(&record, &prepared, |matched| {
            println!(
                "{}\t{}\t{}\t{}\t{}",
                matched.record, matched.strand, matched.start, matched.end, matched.matched
            );
            Ok(())
        })?;
    }

    Ok(())
}
