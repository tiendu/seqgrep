VENV := .venv
PYTHON := $(VENV)/bin/python
PIP := $(VENV)/bin/pip
SEQGREP := $(VENV)/bin/seqgrep

.PHONY: install test lint typecheck check clean run help

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

typecheck:
	$(PYTHON) -m mypy src

check: test lint typecheck

clean:
	rm -rf $(VENV)
	rm -rf .pytest_cache .ruff_cache .mypy_cache
	rm -rf build dist *.egg-info src/*.egg-info
	find . -type d -name "__pycache__" -exec rm -rf {} +

help:
	@echo "Available targets:"
	@echo "  make install     Create venv and install package with dev dependencies"
	@echo "  make run         Run seqgrep CLI"
	@echo "  make test        Run pytest"
	@echo "  make lint        Run ruff"
	@echo "  make typecheck   Run mypy"
	@echo "  make check       Run test, lint, and typecheck"
	@echo "  make clean       Remove venv and build/cache files"
