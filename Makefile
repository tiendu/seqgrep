SHELL := /bin/bash

CARGO ?= cargo
PREFIX ?= $(HOME)/.local
BINDIR ?= $(PREFIX)/bin
SEQGREP := target/release/seqgrep
ARGS ?=

.PHONY: help build release run install uninstall test test-chr21 lint \
	format format-check doc package check clean

help:
	@printf '%s\n' \
		'make build         Build the debug binary and library' \
		'make release       Build the optimized release binary' \
		'make run           Run seqgrep; pass arguments with ARGS="..."' \
		'make install       Install the release binary under ~/.local/bin' \
		'make uninstall     Remove the installed binary' \
		'make test          Run unit and integration tests' \
		'make test-chr21    Run the GRCh38 chromosome 21 integration test' \
		'make lint          Run Clippy with warnings denied' \
		'make format        Format all Rust source' \
		'make format-check  Verify formatting without changing files' \
		'make doc           Build library documentation with warnings denied' \
		'make package       Validate the publishable Cargo package' \
		'make check         Run formatting, linting, tests, and docs' \
		'make clean         Remove Cargo build output'

build:
	$(CARGO) build --locked --all-targets

release:
	$(CARGO) build --locked --release

run:
	$(CARGO) run --locked -- $(ARGS)

install: release
	install -d "$(BINDIR)"
	install -m 0755 "$(SEQGREP)" "$(BINDIR)/seqgrep"

uninstall:
	rm -f "$(BINDIR)/seqgrep"

test:
	$(CARGO) test --locked --all-targets
	$(CARGO) test --locked --doc

test-chr21: release
	SEQGREP_BIN="$(CURDIR)/$(SEQGREP)" tests/test_chr21.sh

lint:
	$(CARGO) clippy --locked --all-targets --all-features -- -D warnings

format:
	$(CARGO) fmt --all

format-check:
	$(CARGO) fmt --all -- --check

doc:
	RUSTDOCFLAGS="-D warnings" $(CARGO) doc --locked --no-deps

package:
	$(CARGO) package --locked --allow-dirty

check: format-check lint test doc

clean:
	$(CARGO) clean
