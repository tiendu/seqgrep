VENV := .venv
PYTHON := $(VENV)/bin/python
PIP := $(VENV)/bin/pip
SEQGREP := $(VENV)/bin/seqgrep

.PHONY: install test lint format format-check typecheck check clean run help

$(PYTHON):
	python3 -m venv $(VENV)
	$(PYTHON) -m pip install --upgrade pip

install: $(PYTHON)
	$(PIP) install -e ".[dev]"

run:
	$(SEQGREP)

test:
	$(PYTHON) -m pytest

lint:
	$(PYTHON) -m ruff check .

format:
	$(PYTHON) -m ruff format .

format-check:
	$(PYTHON) -m ruff format --check .

typecheck:
	$(PYTHON) -m mypy src

check: test lint format-check typecheck

clean:
	rm -rf $(VENV)
	rm -rf .pytest_cache .ruff_cache .mypy_cache
	rm -rf build dist *.egg-info src/*.egg-info
	find . -type d -name "__pycache__" -exec rm -rf {} +

help:
	@echo "Available targets:"
	@echo "  make install       Create venv and install development dependencies"
	@echo "  make run           Run the seqgrep CLI"
	@echo "  make test          Run pytest"
	@echo "  make lint          Run Ruff linting"
	@echo "  make format        Format source and tests"
	@echo "  make format-check  Check formatting"
	@echo "  make typecheck     Run strict mypy"
	@echo "  make check         Run all checks"
	@echo "  make clean         Remove environments, builds, and caches"
