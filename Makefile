SHELL := /bin/bash

VENV := .venv
PYTHON := $(VENV)/bin/python
PIP := $(PYTHON) -m pip
SEQGREP := $(VENV)/bin/seqgrep
INSTALL_STAMP := $(VENV)/.seqgrep-installed

BOOTSTRAP_PYTHON ?= python3
ARGS ?=

.PHONY: help venv install run test test-chr21 lint format format-check \
	typecheck check clean

help:
	@printf '%s\n' \
		'make install        Create or repair .venv and install development dependencies' \
		'make run            Run seqgrep; pass arguments with ARGS="..."' \
		'make test           Run the unit test suite' \
		'make test-chr21     Run the chromosome 21 integration test' \
		'make lint           Run Ruff linting' \
		'make format         Format source and tests' \
		'make format-check   Check formatting without changing files' \
		'make typecheck      Run strict mypy checks' \
		'make check          Run tests, linting, formatting, and type checking' \
		'make clean          Remove the virtual environment and generated artifacts'

venv:
	@if [ ! -x "$(PYTHON)" ] || ! "$(PYTHON)" -m pip --version >/dev/null 2>&1; then \
		printf '%s\n' "Creating or repairing $(VENV)"; \
		rm -rf "$(VENV)"; \
		"$(BOOTSTRAP_PYTHON)" -m venv "$(VENV)"; \
		if ! "$(PYTHON)" -m pip --version >/dev/null 2>&1; then \
			"$(PYTHON)" -m ensurepip --upgrade; \
		fi; \
		"$(PYTHON)" -m pip install --upgrade pip; \
	fi

install: venv
	@if [ ! -f "$(INSTALL_STAMP)" ] \
		|| [ ! -x "$(SEQGREP)" ] \
		|| [ pyproject.toml -nt "$(INSTALL_STAMP)" ]; then \
		printf '%s\n' "Installing seqgrep development environment"; \
		$(PIP) install -e ".[dev]"; \
		touch "$(INSTALL_STAMP)"; \
	fi

run: install
	$(SEQGREP) $(ARGS)

test: install
	$(PYTHON) -m pytest

test-chr21: install
	SEQGREP_BIN="$(SEQGREP)" tests/test_chr21.sh

lint: install
	$(PYTHON) -m ruff check .

format: install
	$(PYTHON) -m ruff format .
	$(PYTHON) -m ruff check --fix .

format-check: install
	$(PYTHON) -m ruff format --check .

typecheck: install
	$(PYTHON) -m mypy src

check: test lint format-check typecheck

clean:
	rm -rf "$(VENV)"
	rm -rf .pytest_cache .ruff_cache .mypy_cache
	rm -rf build dist *.egg-info src/*.egg-info
	find . -type d -name "__pycache__" -prune -exec rm -rf {} +
	find . -type f \( -name "*.pyc" -o -name "*.pyo" \) -delete
