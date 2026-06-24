use std::path::PathBuf;

use seqgrep::{AmbigMode, InputFormat, SearchQuery, SequenceType};

pub(crate) enum Command {
    Help,
    Version,
    Run(Arguments),
}

pub(crate) struct Arguments {
    pub(crate) input: PathBuf,
    pub(crate) query: SearchQuery,
    pub(crate) format: InputFormat,
    pub(crate) jobs: usize,
    pub(crate) chunk_size: usize,
    pub(crate) with_header: bool,
}

pub(crate) fn parse<I>(arguments: I) -> Result<Command, String>
where
    I: IntoIterator<Item = String>,
{
    let mut arguments = arguments.into_iter().peekable();
    let mut positionals = Vec::new();
    let mut sequence_type = SequenceType::Nucleotide;
    let mut ambig_mode = AmbigMode::None;
    let mut ambig_alias = false;
    let mut reverse_complement = false;
    let mut circular = false;
    let mut with_header = false;
    let mut format = InputFormat::Auto;
    let mut jobs = 1_usize;
    let mut chunk_size = 1_000_000_usize;
    let mut options_enabled = true;

    while let Some(argument) = arguments.next() {
        if options_enabled && argument == "--" {
            options_enabled = false;
            continue;
        }

        if options_enabled {
            if argument == "-h" || argument == "--help" {
                return Ok(Command::Help);
            }
            if argument == "--version" {
                return Ok(Command::Version);
            }
            if argument == "--ambig" {
                ambig_alias = true;
                continue;
            }
            if argument == "--revcomp" {
                reverse_complement = true;
                continue;
            }
            if argument == "--circular" {
                circular = true;
                continue;
            }
            if argument == "--with-header" {
                with_header = true;
                continue;
            }

            let (name, inline_value) = split_long_option(&argument);
            let value = match name {
                "-t" | "--sequence-type" | "--ambig-mode" | "--format" | "-j" | "--jobs"
                | "--chunk-size" => Some(match inline_value {
                    Some(value) => value.to_owned(),
                    None => arguments
                        .next()
                        .ok_or_else(|| format!("missing value for {name}"))?,
                }),
                _ => None,
            };

            if let Some(value) = value {
                match name {
                    "-t" | "--sequence-type" => {
                        sequence_type = value
                            .parse()
                            .map_err(|error: seqgrep::Error| error.to_string())?;
                    }
                    "--ambig-mode" => {
                        ambig_mode = value
                            .parse()
                            .map_err(|error: seqgrep::Error| error.to_string())?;
                    }
                    "--format" => {
                        format = value
                            .parse()
                            .map_err(|error: seqgrep::Error| error.to_string())?;
                    }
                    "-j" | "--jobs" => jobs = parse_positive_usize(name, &value)?,
                    "--chunk-size" => chunk_size = parse_positive_usize(name, &value)?,
                    _ => unreachable!(),
                }
                continue;
            }

            if argument.starts_with('-') {
                return Err(format!("unknown option: {argument}"));
            }
        }

        positionals.push(argument);
    }

    if positionals.len() != 2 {
        return Err("expected PATTERN and INPUT positional arguments".to_owned());
    }

    if ambig_alias {
        if ambig_mode == AmbigMode::Both {
            return Err("--ambig cannot be combined with --ambig-mode both".to_owned());
        }
        ambig_mode = AmbigMode::Query;
    }

    if sequence_type == SequenceType::AminoAcid {
        if ambig_mode != AmbigMode::None {
            return Err("--ambig-mode is only valid for nucleotide sequences".to_owned());
        }
        if reverse_complement {
            return Err("--revcomp is only valid for nucleotide sequences".to_owned());
        }
    }

    let query = SearchQuery::new(&positionals[0])
        .with_sequence_type(sequence_type)
        .with_ambig_mode(ambig_mode)
        .with_reverse_complement(reverse_complement)
        .with_circular(circular);

    Ok(Command::Run(Arguments {
        input: PathBuf::from(&positionals[1]),
        query,
        format,
        jobs,
        chunk_size,
        with_header,
    }))
}

fn split_long_option(argument: &str) -> (&str, Option<&str>) {
    if let Some((name, value)) = argument.split_once('=') {
        (name, Some(value))
    } else {
        (argument, None)
    }
}

fn parse_positive_usize(name: &str, value: &str) -> Result<usize, String> {
    let parsed = value
        .parse::<usize>()
        .map_err(|_| format!("{name} must be a positive integer"))?;
    if parsed == 0 {
        return Err(format!("{name} must be at least 1"));
    }
    Ok(parsed)
}

pub(crate) const HELP: &str = r#"Search nucleotide or amino-acid FASTA/FASTQ records.

Usage:
  seqgrep [OPTIONS] PATTERN INPUT

Arguments:
  PATTERN                     Sequence pattern, e.g. ATGNNRY or MTEYK
  INPUT                       FASTA/FASTQ input, optionally gzip-compressed

Options:
  -t, --sequence-type TYPE    nucleotide or amino-acid [default: nucleotide]
      --ambig-mode MODE       none, query, or both [default: none]
      --ambig                 Alias for --ambig-mode query
      --revcomp               Also search the nucleotide reverse complement
      --circular              Allow matches to cross the sequence boundary
      --with-header           Print a TSV header
      --format FORMAT         auto, fasta, or fastq [default: auto]
  -j, --jobs JOBS             Worker threads per record [default: 1]
      --chunk-size N          Candidate starts per work chunk [default: 1000000]
  -h, --help                  Print help
      --version               Print version
"#;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_to_exact_nucleotide() {
        let Command::Run(args) = parse(["ATG".to_owned(), "seq.fa".to_owned()]).unwrap() else {
            panic!("expected run command");
        };
        assert_eq!(args.query.sequence_type(), SequenceType::Nucleotide);
        assert_eq!(args.query.ambig_mode(), AmbigMode::None);
        assert_eq!(args.jobs, 1);
    }

    #[test]
    fn ambig_alias_selects_query_mode() {
        let Command::Run(args) =
            parse(["ATN".to_owned(), "seq.fa".to_owned(), "--ambig".to_owned()]).unwrap()
        else {
            panic!("expected run command");
        };
        assert_eq!(args.query.ambig_mode(), AmbigMode::Query);
    }

    #[test]
    fn rejects_protein_reverse_complement() {
        let error = parse([
            "MTE".to_owned(),
            "protein.fa".to_owned(),
            "-t".to_owned(),
            "amino-acid".to_owned(),
            "--revcomp".to_owned(),
        ])
        .err()
        .unwrap();
        assert!(error.contains("nucleotide"));
    }
}
