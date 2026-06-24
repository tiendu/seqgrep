use std::fmt;
use std::str::FromStr;
use std::sync::mpsc::sync_channel;
use std::thread;

use crate::alphabet::normalize_symbols;
use crate::codec::{Codec, EncodedTarget};
use crate::error::{Error, Result};

/// Biological alphabet used for matching.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum SequenceType {
    /// Nucleotide symbols, including IUPAC codes and gaps.
    #[default]
    Nucleotide,
    /// Amino-acid symbols.
    AminoAcid,
}

impl fmt::Display for SequenceType {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Nucleotide => "nucleotide",
            Self::AminoAcid => "amino-acid",
        })
    }
}

impl FromStr for SequenceType {
    type Err = Error;

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "nucleotide" => Ok(Self::Nucleotide),
            "amino-acid" => Ok(Self::AminoAcid),
            _ => Err(Error::InvalidOption(format!(
                "invalid sequence type {value:?}; expected nucleotide or amino-acid"
            ))),
        }
    }
}

/// Interpretation of IUPAC nucleotide ambiguity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum AmbigMode {
    /// Match symbols exactly. `N` matches only `N`.
    #[default]
    None,
    /// Expand ambiguity in the query only. Ambiguous target symbols match nothing.
    Query,
    /// Expand ambiguity in both query and target.
    Both,
}

impl fmt::Display for AmbigMode {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::None => "none",
            Self::Query => "query",
            Self::Both => "both",
        })
    }
}

impl FromStr for AmbigMode {
    type Err = Error;

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "none" => Ok(Self::None),
            "query" => Ok(Self::Query),
            "both" => Ok(Self::Both),
            _ => Err(Error::InvalidOption(format!(
                "invalid ambiguity mode {value:?}; expected none, query, or both"
            ))),
        }
    }
}

/// Search strand reported for a match.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Strand {
    /// Forward query orientation.
    Forward,
    /// Reverse-complement query orientation.
    Reverse,
}

impl fmt::Display for Strand {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Forward => "+",
            Self::Reverse => "-",
        })
    }
}

/// One normalized biological sequence record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceRecord {
    name: String,
    sequence: String,
}

impl SequenceRecord {
    /// Create a normalized sequence record.
    pub fn new(name: impl Into<String>, sequence: impl AsRef<str>) -> Self {
        Self {
            name: name.into(),
            sequence: normalize_symbols(sequence.as_ref()),
        }
    }

    /// Return the record identifier without FASTA/FASTQ description text.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Return the normalized uppercase sequence.
    pub fn sequence(&self) -> &str {
        &self.sequence
    }

    /// Consume the record and return its identifier and normalized sequence.
    pub fn into_parts(self) -> (String, String) {
        (self.name, self.sequence)
    }
}

/// Search request for one pattern.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SearchQuery {
    pattern: String,
    reverse_complement: bool,
    circular: bool,
    sequence_type: SequenceType,
    ambig_mode: AmbigMode,
}

impl SearchQuery {
    /// Create an exact nucleotide query.
    pub fn new(pattern: impl AsRef<str>) -> Self {
        Self {
            pattern: normalize_symbols(pattern.as_ref()),
            reverse_complement: false,
            circular: false,
            sequence_type: SequenceType::Nucleotide,
            ambig_mode: AmbigMode::None,
        }
    }

    /// Set the biological alphabet.
    #[must_use]
    pub fn with_sequence_type(mut self, sequence_type: SequenceType) -> Self {
        self.sequence_type = sequence_type;
        self
    }

    /// Set nucleotide ambiguity semantics.
    #[must_use]
    pub fn with_ambig_mode(mut self, mode: AmbigMode) -> Self {
        self.ambig_mode = mode;
        self
    }

    /// Enable or disable reverse-complement search.
    #[must_use]
    pub fn with_reverse_complement(mut self, enabled: bool) -> Self {
        self.reverse_complement = enabled;
        self
    }

    /// Enable or disable circular search.
    #[must_use]
    pub fn with_circular(mut self, enabled: bool) -> Self {
        self.circular = enabled;
        self
    }

    /// Return the normalized search pattern.
    pub fn pattern(&self) -> &str {
        &self.pattern
    }

    /// Return whether reverse-complement search is enabled.
    pub fn reverse_complement_enabled(&self) -> bool {
        self.reverse_complement
    }

    /// Return whether circular matching is enabled.
    pub fn circular(&self) -> bool {
        self.circular
    }

    /// Return the selected biological alphabet.
    pub fn sequence_type(&self) -> SequenceType {
        self.sequence_type
    }

    /// Return the nucleotide ambiguity mode.
    pub fn ambig_mode(&self) -> AmbigMode {
        self.ambig_mode
    }

    /// Validate options and biological symbols without searching a record.
    pub fn validate(&self) -> Result<()> {
        PreparedQuery::new(self.clone()).map(|_| ())
    }

    fn validate_options(&self) -> Result<()> {
        if self.pattern.is_empty() {
            return Err(Error::InvalidOption("pattern must not be empty".to_owned()));
        }

        if self.sequence_type == SequenceType::AminoAcid {
            if self.ambig_mode != AmbigMode::None {
                return Err(Error::InvalidOption(
                    "--ambig-mode is only valid for nucleotide sequences".to_owned(),
                ));
            }
            if self.reverse_complement {
                return Err(Error::InvalidOption(
                    "--revcomp is only valid for nucleotide sequences".to_owned(),
                ));
            }
        }

        Ok(())
    }
}

/// Validated and encoded query reusable across many sequence records.
#[derive(Debug, Clone)]
pub struct PreparedQuery {
    query: SearchQuery,
    codec: Codec,
    exact: bool,
    patterns: Vec<PreparedPattern>,
}

#[derive(Debug, Clone)]
struct PreparedPattern {
    exact_symbols: Option<Vec<u8>>,
    encoded_symbols: Vec<u8>,
    strand: Strand,
}

impl PreparedQuery {
    /// Validate and prepare a query for repeated searches.
    pub fn new(query: SearchQuery) -> Result<Self> {
        query.validate_options()?;
        let codec = select_codec(&query);
        let exact = query.ambig_mode == AmbigMode::None;
        let mut raw_patterns = vec![(query.pattern.clone(), Strand::Forward)];
        if query.reverse_complement {
            raw_patterns.push((codec.reverse_complement(&query.pattern)?, Strand::Reverse));
        }

        let mut patterns = Vec::with_capacity(raw_patterns.len());
        for (pattern, strand) in raw_patterns {
            let encoded_symbols = codec.encode_query(&pattern)?;
            let exact_symbols = if exact {
                Some(codec.comparison_bytes(&pattern)?)
            } else {
                None
            };
            patterns.push(PreparedPattern {
                exact_symbols,
                encoded_symbols,
                strand,
            });
        }

        Ok(Self {
            query,
            codec,
            exact,
            patterns,
        })
    }

    /// Return the original normalized query configuration.
    pub fn query(&self) -> &SearchQuery {
        &self.query
    }
}

/// One 1-based inclusive search result.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct Match {
    /// Source record identifier.
    pub record: String,
    /// Query orientation.
    pub strand: Strand,
    /// One-based inclusive start.
    pub start: usize,
    /// One-based inclusive end; wraps for circular matches.
    pub end: usize,
    /// Sequence observed in the target.
    pub matched: String,
    /// Whether the match crosses the target boundary.
    pub circular: bool,
}

/// Parallel execution settings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExecutionConfig {
    jobs: usize,
    chunk_size: usize,
}

impl Default for ExecutionConfig {
    fn default() -> Self {
        Self {
            jobs: 1,
            chunk_size: 1_000_000,
        }
    }
}

impl ExecutionConfig {
    /// Validate and create execution settings.
    pub fn new(jobs: usize, chunk_size: usize) -> Result<Self> {
        if jobs == 0 {
            return Err(Error::InvalidOption("jobs must be at least 1".to_owned()));
        }
        if chunk_size == 0 {
            return Err(Error::InvalidOption(
                "chunk_size must be at least 1".to_owned(),
            ));
        }
        Ok(Self { jobs, chunk_size })
    }

    /// Return the worker-thread count used per record.
    pub fn jobs(self) -> usize {
        self.jobs
    }

    /// Return the number of candidate starts in one work chunk.
    pub fn chunk_size(self) -> usize {
        self.chunk_size
    }
}

/// Reusable search backend.
#[derive(Debug, Clone, Default)]
pub struct SearchEngine {
    config: ExecutionConfig,
}

impl SearchEngine {
    /// Create a search engine with validated execution settings.
    pub const fn new(config: ExecutionConfig) -> Self {
        Self { config }
    }

    /// Return the engine's execution settings.
    pub const fn config(&self) -> ExecutionConfig {
        self.config
    }

    /// Validate and encode a query for efficient reuse across records.
    pub fn prepare_query(&self, query: SearchQuery) -> Result<PreparedQuery> {
        PreparedQuery::new(query)
    }

    /// Search one record and return deterministic, coordinate-sorted matches.
    ///
    /// This convenience method prepares the query and collects all matches.
    /// Use [`Self::search_prepared`] when the same query is applied to many
    /// records, and [`Self::visit_prepared_matches`] for incremental output.
    pub fn search_record(
        &self,
        record: &SequenceRecord,
        query: &SearchQuery,
    ) -> Result<Vec<Match>> {
        let prepared = PreparedQuery::new(query.clone())?;
        self.search_prepared(record, &prepared)
    }

    /// Search one record with a previously prepared query and collect matches.
    pub fn search_prepared(
        &self,
        record: &SequenceRecord,
        prepared: &PreparedQuery,
    ) -> Result<Vec<Match>> {
        let mut matches = Vec::new();
        self.visit_prepared_matches(record, prepared, |matched| {
            matches.push(matched);
            Ok(())
        })?;
        Ok(matches)
    }

    /// Visit matches without collecting them.
    ///
    /// This convenience method prepares the query for this call. Prefer
    /// [`Self::visit_prepared_matches`] when searching multiple records.
    pub fn visit_matches<F>(
        &self,
        record: &SequenceRecord,
        query: &SearchQuery,
        visitor: F,
    ) -> Result<()>
    where
        F: FnMut(Match) -> Result<()>,
    {
        let prepared = PreparedQuery::new(query.clone())?;
        self.visit_prepared_matches(record, &prepared, visitor)
    }

    /// Visit matches from a prepared query in deterministic coordinate order.
    ///
    /// Forward-strand matches are emitted first, followed by reverse-strand
    /// matches. Parallel work is buffered by chunk so output remains ordered
    /// and memory use stays bounded.
    pub fn visit_prepared_matches<F>(
        &self,
        record: &SequenceRecord,
        prepared: &PreparedQuery,
        mut visitor: F,
    ) -> Result<()>
    where
        F: FnMut(Match) -> Result<()>,
    {
        if record.sequence.is_empty() {
            return Ok(());
        }

        if self.config.jobs == 1 && prepared.exact {
            return self.visit_exact(record, prepared, &mut visitor);
        }

        self.visit_encoded(record, prepared, &mut visitor)
    }

    fn visit_exact<F>(
        &self,
        record: &SequenceRecord,
        prepared: &PreparedQuery,
        visitor: &mut F,
    ) -> Result<()>
    where
        F: FnMut(Match) -> Result<()>,
    {
        let sequence_key = prepared.codec.comparison_bytes(&record.sequence)?;

        for pattern in &prepared.patterns {
            let symbols = pattern
                .exact_symbols
                .as_deref()
                .ok_or_else(|| Error::InvalidOption("exact query was not prepared".to_owned()))?;
            visit_kmp_starts(&sequence_key, symbols, prepared.query.circular, |start| {
                visitor(build_match(
                    record,
                    pattern.strand,
                    start,
                    symbols.len(),
                    prepared.query.circular,
                )?)
            })?;
        }
        Ok(())
    }

    fn visit_encoded<F>(
        &self,
        record: &SequenceRecord,
        prepared: &PreparedQuery,
        visitor: &mut F,
    ) -> Result<()>
    where
        F: FnMut(Match) -> Result<()>,
    {
        let target = prepared.codec.encode_target(&record.sequence)?;
        if target.len() != record.sequence.len() {
            return Err(Error::InvalidOption(
                "encoded target length does not match sequence length".to_owned(),
            ));
        }

        for pattern in &prepared.patterns {
            visit_scan_starts(
                &target,
                &pattern.encoded_symbols,
                prepared.codec,
                prepared.query.circular,
                self.config,
                |start| {
                    visitor(build_match(
                        record,
                        pattern.strand,
                        start,
                        pattern.encoded_symbols.len(),
                        prepared.query.circular,
                    )?)
                },
            )?;
        }
        Ok(())
    }
}

fn select_codec(query: &SearchQuery) -> Codec {
    match (query.sequence_type, query.ambig_mode) {
        (SequenceType::AminoAcid, _) => Codec::ProteinExact,
        (SequenceType::Nucleotide, AmbigMode::None) => Codec::NucleotideExact,
        (SequenceType::Nucleotide, AmbigMode::Query) => Codec::NucleotideIupac {
            target_ambiguity: false,
        },
        (SequenceType::Nucleotide, AmbigMode::Both) => Codec::NucleotideIupac {
            target_ambiguity: true,
        },
    }
}

fn build_match(
    record: &SequenceRecord,
    strand: Strand,
    zero_start: usize,
    pattern_len: usize,
    circular: bool,
) -> Result<Match> {
    let sequence_len = record.sequence.len();
    let zero_end = zero_start
        .checked_add(pattern_len - 1)
        .ok_or_else(|| Error::InvalidOption("match coordinate overflow".to_owned()))?;
    let wraps = circular && zero_end >= sequence_len;
    let end = if circular {
        zero_end % sequence_len + 1
    } else {
        zero_end + 1
    };

    Ok(Match {
        record: record.name.clone(),
        strand,
        start: zero_start + 1,
        end,
        matched: matched_sequence(&record.sequence, zero_start, pattern_len, circular)?,
        circular: wraps,
    })
}

fn matched_sequence(sequence: &str, start: usize, length: usize, circular: bool) -> Result<String> {
    if !circular {
        return sequence
            .get(start..start + length)
            .map(ToOwned::to_owned)
            .ok_or_else(|| Error::InvalidOption("invalid match bounds".to_owned()));
    }

    let bytes = sequence.as_bytes();
    let mut matched = Vec::with_capacity(length);
    let mut index = start;
    for _ in 0..length {
        matched.push(bytes[index]);
        index += 1;
        if index == bytes.len() {
            index = 0;
        }
    }
    String::from_utf8(matched)
        .map_err(|_| Error::InvalidOption("validated sequence was not ASCII".to_owned()))
}

fn visit_kmp_starts<F>(
    haystack: &[u8],
    pattern: &[u8],
    circular: bool,
    mut visitor: F,
) -> Result<()>
where
    F: FnMut(usize) -> Result<()>,
{
    if pattern.is_empty() {
        return Err(Error::InvalidOption("pattern must not be empty".to_owned()));
    }
    if haystack.is_empty() || (!circular && pattern.len() > haystack.len()) {
        return Ok(());
    }

    let total = if circular {
        haystack
            .len()
            .checked_add(pattern.len() - 1)
            .ok_or_else(|| Error::InvalidOption("circular search length overflow".to_owned()))?
    } else {
        haystack.len()
    };
    let maximum_start = if circular {
        haystack.len()
    } else {
        haystack.len() - pattern.len() + 1
    };
    let prefix = kmp_prefix(pattern);
    let mut matched = 0;

    for index in 0..total {
        let symbol = if circular {
            haystack[index % haystack.len()]
        } else {
            haystack[index]
        };

        while matched > 0 && symbol != pattern[matched] {
            matched = prefix[matched - 1];
        }
        if symbol == pattern[matched] {
            matched += 1;
        }
        if matched == pattern.len() {
            let start = index + 1 - pattern.len();
            if start < maximum_start {
                visitor(start)?;
            }
            matched = prefix[matched - 1];
        }
    }

    Ok(())
}

fn kmp_prefix(pattern: &[u8]) -> Vec<usize> {
    let mut prefix = vec![0; pattern.len()];
    let mut matched = 0;

    for index in 1..pattern.len() {
        while matched > 0 && pattern[index] != pattern[matched] {
            matched = prefix[matched - 1];
        }
        if pattern[index] == pattern[matched] {
            matched += 1;
            prefix[index] = matched;
        }
    }

    prefix
}

fn visit_scan_starts<F>(
    target: &EncodedTarget,
    query: &[u8],
    codec: Codec,
    circular: bool,
    config: ExecutionConfig,
    mut visitor: F,
) -> Result<()>
where
    F: FnMut(usize) -> Result<()>,
{
    if query.is_empty() {
        return Err(Error::InvalidOption("pattern must not be empty".to_owned()));
    }
    if !circular && query.len() > target.len() {
        return Ok(());
    }

    let total_starts = if circular {
        target.len()
    } else {
        target.len() - query.len() + 1
    };
    if total_starts == 0 {
        return Ok(());
    }

    let ranges = chunk_ranges(total_starts, config.chunk_size);
    if config.jobs == 1 || ranges.len() == 1 {
        for (begin, end) in ranges {
            for start in scan_range(target, query, codec, circular, begin, end) {
                visitor(start)?;
            }
        }
        return Ok(());
    }

    let worker_count = config.jobs.min(ranges.len());
    let partitions = contiguous_partitions(ranges.len(), worker_count);

    thread::scope(|scope| {
        let mut receivers = Vec::with_capacity(worker_count);
        let mut handles = Vec::with_capacity(worker_count);

        for (first, last) in partitions {
            let (sender, receiver) = sync_channel::<Vec<usize>>(1);
            receivers.push(receiver);
            let ranges_ref = &ranges;
            handles.push(scope.spawn(move || {
                for &(begin, end) in &ranges_ref[first..last] {
                    let starts = scan_range(target, query, codec, circular, begin, end);
                    if sender.send(starts).is_err() {
                        return;
                    }
                }
            }));
        }

        let mut visitor_error = None;
        'workers: for receiver in &receivers {
            while let Ok(starts) = receiver.recv() {
                for start in starts {
                    if let Err(error) = visitor(start) {
                        visitor_error = Some(error);
                        break 'workers;
                    }
                }
            }
        }

        drop(receivers);
        let mut worker_panicked = false;
        for handle in handles {
            if handle.join().is_err() {
                worker_panicked = true;
            }
        }

        if worker_panicked {
            return Err(Error::WorkerPanicked);
        }
        if let Some(error) = visitor_error {
            return Err(error);
        }
        Ok(())
    })
}

fn contiguous_partitions(total: usize, parts: usize) -> Vec<(usize, usize)> {
    debug_assert!(parts > 0 && parts <= total);
    let base = total / parts;
    let remainder = total % parts;
    let mut partitions = Vec::with_capacity(parts);
    let mut begin = 0;

    for part in 0..parts {
        let length = base + usize::from(part < remainder);
        let end = begin + length;
        partitions.push((begin, end));
        begin = end;
    }

    partitions
}

fn chunk_ranges(total_starts: usize, chunk_size: usize) -> Vec<(usize, usize)> {
    let mut ranges = Vec::with_capacity(total_starts.div_ceil(chunk_size));
    let mut begin = 0;
    while begin < total_starts {
        let end = begin.saturating_add(chunk_size).min(total_starts);
        ranges.push((begin, end));
        begin = end;
    }
    ranges
}

fn scan_range(
    target: &EncodedTarget,
    query: &[u8],
    codec: Codec,
    circular: bool,
    begin: usize,
    end: usize,
) -> Vec<usize> {
    let mut starts = Vec::new();
    for start in begin..end {
        if window_matches(target, query, codec, start, circular) {
            starts.push(start);
        }
    }
    starts
}

fn window_matches(
    target: &EncodedTarget,
    query: &[u8],
    codec: Codec,
    start: usize,
    circular: bool,
) -> bool {
    let mut index = start;
    for query_symbol in query.iter().copied() {
        if !codec.compatible(query_symbol, target.symbol_at(index)) {
            return false;
        }
        index += 1;
        if circular && index == target.len() {
            index = 0;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit_values(matches: &[Match]) -> Vec<(usize, usize, &str)> {
        matches
            .iter()
            .map(|item| (item.start, item.end, item.matched.as_str()))
            .collect()
    }

    #[test]
    fn exact_search_finds_overlaps() {
        let engine = SearchEngine::default();
        let matches = engine
            .search_record(&SequenceRecord::new("seq", "AAAA"), &SearchQuery::new("AA"))
            .unwrap();
        assert_eq!(
            hit_values(&matches),
            vec![(1, 2, "AA"), (2, 3, "AA"), (3, 4, "AA")]
        );
    }

    #[test]
    fn circular_pattern_can_exceed_sequence_length() {
        let engine = SearchEngine::default();
        let query = SearchQuery::new("ATGAT").with_circular(true);
        let matches = engine
            .search_record(&SequenceRecord::new("tiny", "ATG"), &query)
            .unwrap();
        assert_eq!(hit_values(&matches), vec![(1, 2, "ATGAT")]);
    }

    #[test]
    fn query_ambiguity_ignores_unknown_target_regions() {
        let engine = SearchEngine::default();
        let query = SearchQuery::new("NNNN").with_ambig_mode(AmbigMode::Query);
        let matches = engine
            .search_record(&SequenceRecord::new("gap", "N".repeat(64)), &query)
            .unwrap();
        assert!(matches.is_empty());
    }

    #[test]
    fn both_ambiguity_matches_unknown_target_regions() {
        let engine = SearchEngine::default();
        let query = SearchQuery::new("AAAAA").with_ambig_mode(AmbigMode::Both);
        let matches = engine
            .search_record(&SequenceRecord::new("uncertain", "NNANN"), &query)
            .unwrap();
        assert_eq!(hit_values(&matches), vec![(1, 5, "NNANN")]);
    }

    #[test]
    fn serial_and_parallel_results_are_identical() {
        let record = SequenceRecord::new("seq", "AAAAATGCAAAAA");
        let query = SearchQuery::new("TGCA");
        let serial = SearchEngine::default()
            .search_record(&record, &query)
            .unwrap();
        let parallel = SearchEngine::new(ExecutionConfig::new(3, 2).unwrap())
            .search_record(&record, &query)
            .unwrap();
        assert_eq!(parallel, serial);
    }
}
