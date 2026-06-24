use std::fs;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn unique_path(name: &str) -> std::path::PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!("seqgrep-cli-{}-{nonce}-{name}", std::process::id()))
}

#[test]
fn binary_prints_expected_tsv() {
    let input = unique_path("input.fa");
    fs::write(&input, b">seq\nATGNATGA\n").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_seqgrep"))
        .args(["ATGN", input.to_str().unwrap(), "--ambig-mode", "query"])
        .output()
        .unwrap();

    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        "seq\t+\t5\t8\tATGA\tfalse\n"
    );
    fs::remove_file(input).unwrap();
}

#[test]
fn binary_rejects_invalid_protein_options() {
    let input = unique_path("protein.fa");
    fs::write(&input, b">protein\nMTEYK\n").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_seqgrep"))
        .args([
            "MTE",
            input.to_str().unwrap(),
            "--sequence-type",
            "amino-acid",
            "--revcomp",
        ])
        .output()
        .unwrap();

    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8(output.stderr)
        .unwrap()
        .contains("nucleotide"));
    fs::remove_file(input).unwrap();
}

#[test]
fn binary_prints_header_when_requested() {
    let input = unique_path("header.fa");
    fs::write(&input, b">seq\nATG\n").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_seqgrep"))
        .args(["ATG", input.to_str().unwrap(), "--with-header"])
        .output()
        .unwrap();

    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        "record\tstrand\tstart\tend\tmatched\tcircular\nseq\t+\t1\t3\tATG\tfalse\n"
    );
    fs::remove_file(input).unwrap();
}

#[test]
fn binary_rejects_conflicting_ambiguity_options() {
    let input = unique_path("conflict.fa");
    fs::write(&input, b">seq\nATG\n").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_seqgrep"))
        .args([
            "ATG",
            input.to_str().unwrap(),
            "--ambig",
            "--ambig-mode",
            "both",
        ])
        .output()
        .unwrap();

    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8(output.stderr)
        .unwrap()
        .contains("cannot be combined"));
    fs::remove_file(input).unwrap();
}

#[test]
fn binary_rejects_zero_worker_values() {
    let input = unique_path("workers.fa");
    fs::write(&input, b">seq\nATG\n").unwrap();

    for arguments in [
        vec!["ATG", input.to_str().unwrap(), "--jobs", "0"],
        vec!["ATG", input.to_str().unwrap(), "--chunk-size", "0"],
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_seqgrep"))
            .args(arguments)
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(2));
    }

    fs::remove_file(input).unwrap();
}

#[test]
fn binary_reports_version() {
    let output = Command::new(env!("CARGO_BIN_EXE_seqgrep"))
        .arg("--version")
        .output()
        .unwrap();

    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        format!("seqgrep {}\n", seqgrep::VERSION)
    );
}
