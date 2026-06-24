//! Streaming FASTA and FASTQ readers.

use std::ffi::OsStr;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::str::FromStr;

use flate2::read::MultiGzDecoder;

use crate::error::{Error, Result};
use crate::search::SequenceRecord;

/// Input sequence format.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum InputFormat {
    /// Infer format from the filename extension.
    #[default]
    Auto,
    /// FASTA input.
    Fasta,
    /// FASTQ input.
    Fastq,
}

impl fmt::Display for InputFormat {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Auto => "auto",
            Self::Fasta => "fasta",
            Self::Fastq => "fastq",
        })
    }
}

impl FromStr for InputFormat {
    type Err = Error;

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "auto" => Ok(Self::Auto),
            "fasta" => Ok(Self::Fasta),
            "fastq" => Ok(Self::Fastq),
            _ => Err(Error::InvalidOption(format!(
                "invalid input format {value:?}; expected auto, fasta, or fastq"
            ))),
        }
    }
}

/// Streaming FASTA/FASTQ reader.
///
/// FASTQ quality data is validated for length and then discarded because the
/// search backend only needs record identifiers and sequences.
pub struct FastxReader {
    inner: ReaderKind,
}

enum ReaderKind {
    Fasta(FastaReader<Box<dyn BufRead + Send>>),
    Fastq(FastqReader<Box<dyn BufRead + Send>>),
}

impl FastxReader {
    /// Open a plain-text or gzip-compressed FASTA/FASTQ path.
    pub fn from_path(path: impl AsRef<Path>, format: InputFormat) -> Result<Self> {
        let path = path.as_ref();
        let selected = infer_format(path, format)?;
        let file = File::open(path)?;
        let reader: Box<dyn BufRead + Send> = if has_gzip_suffix(path) {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        Self::from_boxed_reader(reader, selected)
    }

    /// Create a reader from any byte stream and an explicit format.
    pub fn from_reader<R>(reader: R, format: InputFormat) -> Result<Self>
    where
        R: Read + Send + 'static,
    {
        if format == InputFormat::Auto {
            return Err(Error::InvalidOption(
                "InputFormat::Auto requires a filesystem path".to_owned(),
            ));
        }
        Self::from_boxed_reader(Box::new(BufReader::new(reader)), format)
    }

    fn from_boxed_reader(reader: Box<dyn BufRead + Send>, format: InputFormat) -> Result<Self> {
        let inner = match format {
            InputFormat::Fasta => ReaderKind::Fasta(FastaReader::new(reader)),
            InputFormat::Fastq => ReaderKind::Fastq(FastqReader::new(reader)),
            InputFormat::Auto => {
                return Err(Error::InvalidOption(
                    "input format must be resolved before constructing a reader".to_owned(),
                ));
            }
        };
        Ok(Self { inner })
    }
}

impl Iterator for FastxReader {
    type Item = Result<SequenceRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match &mut self.inner {
            ReaderKind::Fasta(reader) => reader.next(),
            ReaderKind::Fastq(reader) => reader.next(),
        }
    }
}

/// Infer FASTA or FASTQ from a path, ignoring a final `.gz` suffix.
pub fn infer_format(path: impl AsRef<Path>, format: InputFormat) -> Result<InputFormat> {
    if format != InputFormat::Auto {
        return Ok(format);
    }

    let path = path.as_ref();
    let mut effective = PathBuf::from(path);
    if has_gzip_suffix(path) {
        let file_stem = path.file_stem().ok_or_else(|| {
            Error::UnknownFormat(
                "Could not infer input format. Use --format fasta or --format fastq.".to_owned(),
            )
        })?;
        effective = PathBuf::from(file_stem);
    }

    let suffix = effective
        .extension()
        .and_then(OsStr::to_str)
        .map(str::to_ascii_lowercase)
        .ok_or_else(|| {
            Error::UnknownFormat(
                "Could not infer input format from file extension. Use --format fasta or --format fastq."
                    .to_owned(),
            )
        })?;

    match suffix.as_str() {
        "fa" | "fasta" | "fna" | "ffn" | "faa" | "frn" => Ok(InputFormat::Fasta),
        "fq" | "fastq" => Ok(InputFormat::Fastq),
        _ => Err(Error::UnknownFormat(format!(
            "Could not infer input format from extension {suffix:?}. Use --format fasta or --format fastq."
        ))),
    }
}

fn has_gzip_suffix(path: &Path) -> bool {
    path.extension()
        .and_then(OsStr::to_str)
        .is_some_and(|suffix| suffix.eq_ignore_ascii_case("gz"))
}

struct FastaReader<R> {
    reader: R,
    line: String,
    line_number: usize,
    pending_header: Option<String>,
    finished: bool,
}

impl<R: BufRead> FastaReader<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            line: String::new(),
            line_number: 0,
            pending_header: None,
            finished: false,
        }
    }

    fn read_line(&mut self) -> Result<Option<String>> {
        self.line.clear();
        let count = self.reader.read_line(&mut self.line)?;
        if count == 0 {
            return Ok(None);
        }
        self.line_number += 1;
        Ok(Some(self.line.trim().to_owned()))
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = Result<SequenceRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        let header = if let Some(header) = self.pending_header.take() {
            header
        } else {
            loop {
                match self.read_line() {
                    Ok(Some(line)) if line.is_empty() => continue,
                    Ok(Some(line)) if line.starts_with('>') => break line,
                    Ok(Some(_)) => {
                        self.finished = true;
                        return Some(Err(Error::Parse {
                            format: "FASTA",
                            line: self.line_number,
                            message: "sequence appeared before the first header".to_owned(),
                        }));
                    }
                    Ok(None) => {
                        self.finished = true;
                        return None;
                    }
                    Err(error) => {
                        self.finished = true;
                        return Some(Err(error));
                    }
                }
            }
        };

        let name = header[1..]
            .split_whitespace()
            .next()
            .unwrap_or_default()
            .to_owned();
        if name.is_empty() {
            self.finished = true;
            return Some(Err(Error::Parse {
                format: "FASTA",
                line: self.line_number,
                message: "header is missing a record identifier".to_owned(),
            }));
        }
        let mut sequence = String::new();

        loop {
            match self.read_line() {
                Ok(Some(line)) if line.is_empty() => continue,
                Ok(Some(line)) if line.starts_with('>') => {
                    self.pending_header = Some(line);
                    return Some(Ok(SequenceRecord::new(name, sequence)));
                }
                Ok(Some(line)) => sequence.push_str(&line),
                Ok(None) => {
                    self.finished = true;
                    return Some(Ok(SequenceRecord::new(name, sequence)));
                }
                Err(error) => {
                    self.finished = true;
                    return Some(Err(error));
                }
            }
        }
    }
}

struct FastqReader<R> {
    reader: R,
    line: String,
    line_number: usize,
    finished: bool,
}

impl<R: BufRead> FastqReader<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            line: String::new(),
            line_number: 0,
            finished: false,
        }
    }

    fn read_raw_line(&mut self) -> Result<Option<String>> {
        self.line.clear();
        let count = self.reader.read_line(&mut self.line)?;
        if count == 0 {
            return Ok(None);
        }
        self.line_number += 1;
        Ok(Some(self.line.trim_end_matches(['\n', '\r']).to_owned()))
    }

    fn next_nonempty_line(&mut self) -> Result<Option<String>> {
        loop {
            match self.read_raw_line()? {
                Some(line) if line.is_empty() => continue,
                other => return Ok(other),
            }
        }
    }
}

impl<R: BufRead> Iterator for FastqReader<R> {
    type Item = Result<SequenceRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        let header = match self.next_nonempty_line() {
            Ok(Some(header)) => header,
            Ok(None) => {
                self.finished = true;
                return None;
            }
            Err(error) => {
                self.finished = true;
                return Some(Err(error));
            }
        };

        if !header.starts_with('@') {
            self.finished = true;
            return Some(Err(Error::Parse {
                format: "FASTQ",
                line: self.line_number,
                message: "record must start with '@'".to_owned(),
            }));
        }

        let name = header[1..]
            .split_whitespace()
            .next()
            .unwrap_or_default()
            .to_owned();
        if name.is_empty() {
            self.finished = true;
            return Some(Err(Error::Parse {
                format: "FASTQ",
                line: self.line_number,
                message: "header is missing a record identifier".to_owned(),
            }));
        }
        let mut raw_sequence = String::new();

        loop {
            match self.read_raw_line() {
                Ok(Some(line)) if line.starts_with('+') => break,
                Ok(Some(line)) => raw_sequence.push_str(line.trim()),
                Ok(None) => {
                    self.finished = true;
                    return Some(Err(Error::Parse {
                        format: "FASTQ",
                        line: self.line_number,
                        message: format!("record {name:?} is missing '+' quality header"),
                    }));
                }
                Err(error) => {
                    self.finished = true;
                    return Some(Err(error));
                }
            }
        }

        if raw_sequence.is_empty() {
            self.finished = true;
            return Some(Err(Error::Parse {
                format: "FASTQ",
                line: self.line_number,
                message: format!("record {name:?} has an empty sequence"),
            }));
        }

        let expected_quality = raw_sequence.len();
        let mut quality_length = 0;
        while quality_length < expected_quality {
            match self.read_raw_line() {
                Ok(Some(line)) => quality_length += line.len(),
                Ok(None) => {
                    self.finished = true;
                    return Some(Err(Error::Parse {
                        format: "FASTQ",
                        line: self.line_number,
                        message: format!("record {name:?} ended before quality block was complete"),
                    }));
                }
                Err(error) => {
                    self.finished = true;
                    return Some(Err(error));
                }
            }
        }

        if quality_length != expected_quality {
            self.finished = true;
            return Some(Err(Error::Parse {
                format: "FASTQ",
                line: self.line_number,
                message: format!(
                    "record {name:?} has quality length {quality_length}, expected {expected_quality}"
                ),
            }));
        }

        Some(Ok(SequenceRecord::new(name, raw_sequence)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn reads_multiline_fasta() {
        let input = Cursor::new(b">seq description\nATG\nCCC\n");
        let records = FastxReader::from_reader(input, InputFormat::Fasta)
            .unwrap()
            .collect::<Result<Vec<_>>>()
            .unwrap();
        assert_eq!(records, vec![SequenceRecord::new("seq", "ATGCCC")]);
    }

    #[test]
    fn reads_multiline_fastq() {
        let input = Cursor::new(b"@read description\nAT\nGC\n+\n!!\n!!\n");
        let records = FastxReader::from_reader(input, InputFormat::Fastq)
            .unwrap()
            .collect::<Result<Vec<_>>>()
            .unwrap();
        assert_eq!(records, vec![SequenceRecord::new("read", "ATGC")]);
    }

    #[test]
    fn rejects_short_fastq_quality() {
        let input = Cursor::new(b"@read\nATGC\n+\n!!!\n");
        let error = FastxReader::from_reader(input, InputFormat::Fastq)
            .unwrap()
            .next()
            .unwrap()
            .unwrap_err();
        assert!(error
            .to_string()
            .contains("ended before quality block was complete"));
    }
}
