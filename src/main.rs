mod cli;

use std::env;
use std::io::{self, BufWriter, Write};
use std::process::ExitCode;

use seqgrep::{ExecutionConfig, FastxReader, Match, SearchEngine};

fn main() -> ExitCode {
    match cli::parse(env::args().skip(1)) {
        Ok(cli::Command::Help) => {
            print!("{}", cli::HELP);
            ExitCode::SUCCESS
        }
        Ok(cli::Command::Version) => {
            println!("seqgrep {}", seqgrep::VERSION);
            ExitCode::SUCCESS
        }
        Ok(cli::Command::Run(arguments)) => match run(arguments) {
            Ok(()) => ExitCode::SUCCESS,
            Err(seqgrep::Error::Io(error)) if error.kind() == io::ErrorKind::BrokenPipe => {
                ExitCode::SUCCESS
            }
            Err(error) => {
                eprintln!("seqgrep: {error}");
                ExitCode::from(1)
            }
        },
        Err(message) => {
            eprintln!("seqgrep: {message}\n\n{}", cli::HELP);
            ExitCode::from(2)
        }
    }
}

fn run(arguments: cli::Arguments) -> seqgrep::Result<()> {
    let config = ExecutionConfig::new(arguments.jobs, arguments.chunk_size)?;
    let engine = SearchEngine::new(config);
    let prepared = engine.prepare_query(arguments.query)?;
    let reader = FastxReader::from_path(&arguments.input, arguments.format)?;
    let stdout = io::stdout();
    let mut output = BufWriter::new(stdout.lock());

    if arguments.with_header {
        writeln!(output, "record\tstrand\tstart\tend\tmatched\tcircular")?;
    }

    for record in reader {
        let record = record?;
        engine.visit_prepared_matches(&record, &prepared, |matched| {
            write_match(&mut output, &matched)?;
            Ok(())
        })?;
    }

    output.flush()?;
    Ok(())
}

fn write_match(output: &mut impl Write, matched: &Match) -> io::Result<()> {
    writeln!(
        output,
        "{}\t{}\t{}\t{}\t{}\t{}",
        matched.record,
        matched.strand,
        matched.start,
        matched.end,
        matched.matched,
        matched.circular
    )
}
