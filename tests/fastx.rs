use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::write::GzEncoder;
use flate2::Compression;
use seqgrep::{FastxReader, InputFormat, SequenceRecord};

fn temp_dir() -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!("seqgrep-test-{}-{nonce}", std::process::id()));
    fs::create_dir_all(&path).unwrap();
    path
}

#[test]
fn reads_plain_and_gzip_fasta() {
    let directory = temp_dir();
    let plain = directory.join("example.fa");
    let gzip = directory.join("example.FA.GZ");
    let content = b">seq1 description\nATG\nCCC\n";
    fs::write(&plain, content).unwrap();

    let mut encoder = GzEncoder::new(File::create(&gzip).unwrap(), Compression::default());
    encoder.write_all(content).unwrap();
    encoder.finish().unwrap();

    let expected = vec![SequenceRecord::new("seq1", "ATGCCC")];
    assert_eq!(
        FastxReader::from_path(&plain, InputFormat::Auto)
            .unwrap()
            .collect::<seqgrep::Result<Vec<_>>>()
            .unwrap(),
        expected
    );
    assert_eq!(
        FastxReader::from_path(&gzip, InputFormat::Auto)
            .unwrap()
            .collect::<seqgrep::Result<Vec<_>>>()
            .unwrap(),
        expected
    );
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn reads_multiline_fastq() {
    let directory = temp_dir();
    let path = directory.join("reads.fastq");
    fs::write(&path, b"@read1 description\nAT\nGC\n+\n!!\n!!\n").unwrap();
    assert_eq!(
        FastxReader::from_path(&path, InputFormat::Auto)
            .unwrap()
            .collect::<seqgrep::Result<Vec<_>>>()
            .unwrap(),
        vec![SequenceRecord::new("read1", "ATGC")]
    );
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn rejects_quality_length_mismatch() {
    let directory = temp_dir();
    let path = directory.join("bad.fastq");
    fs::write(&path, b"@read1\nATGC\n+\n!!!\n").unwrap();
    let error = FastxReader::from_path(&path, InputFormat::Auto)
        .unwrap()
        .next()
        .unwrap()
        .unwrap_err();
    assert!(error.to_string().contains("quality block was complete"));
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn from_reader_requires_an_explicit_format() {
    let error = FastxReader::from_reader(std::io::Cursor::new(b">seq\nATG\n"), InputFormat::Auto)
        .err()
        .expect("auto format should fail without a path");
    assert!(error.to_string().contains("requires a filesystem path"));
}

#[test]
fn rejects_sequence_before_fasta_header() {
    let error = FastxReader::from_reader(
        std::io::Cursor::new(b"ATG\n>seq\nATG\n"),
        InputFormat::Fasta,
    )
    .unwrap()
    .next()
    .unwrap()
    .unwrap_err();
    assert!(error.to_string().contains("before the first header"));
}

#[test]
fn rejects_fastq_without_quality_header() {
    let error =
        FastxReader::from_reader(std::io::Cursor::new(b"@read\nATGC\n"), InputFormat::Fastq)
            .unwrap()
            .next()
            .unwrap()
            .unwrap_err();
    assert!(error.to_string().contains("missing '+'"));
}

#[test]
fn recognizes_amino_acid_fasta_suffix() {
    let directory = temp_dir();
    let path = directory.join("proteins.faa");
    fs::write(&path, b">protein\nMTEYK\n").unwrap();

    let records = FastxReader::from_path(&path, InputFormat::Auto)
        .unwrap()
        .collect::<seqgrep::Result<Vec<_>>>()
        .unwrap();
    assert_eq!(records, vec![SequenceRecord::new("protein", "MTEYK")]);
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn rejects_empty_record_identifiers() {
    for (input, format) in [
        (&b">\nATG\n"[..], InputFormat::Fasta),
        (&b"@\nATG\n+\n!!!\n"[..], InputFormat::Fastq),
    ] {
        let error = FastxReader::from_reader(std::io::Cursor::new(input), format)
            .unwrap()
            .next()
            .unwrap()
            .unwrap_err();
        assert!(error.to_string().contains("record identifier"));
    }
}
