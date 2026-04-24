# Contributing to Quantarhei

Thank you for your interest in contributing to Quantarhei!

## Quick setup

```bash
git clone https://github.com/tmancal74/quantarhei
cd quantarhei
pip install -e .
pip install -r requirements_devel.txt
pre-commit install
```

`pre-commit install` wires up the code quality hooks so they run automatically on every commit. You only need to do this once per clone.

## Running tests

```bash
# Unit tests — fast feedback, run before every push
pytest tests/unit

# Doc tests
pytest --doctest-modules \
  quantarhei/core quantarhei/builders quantarhei/qm/corfunctions \
  quantarhei/spectroscopy quantarhei/qm/liouvillespace quantarhei/functions \
  quantarhei/qm/hilbertspace quantarhei/qm/propagators

# Behave acceptance tests
coverage run -m behave tests/behave/features
```

CI runs all three suites on Python 3.10, 3.11, and 3.12. Please verify unit tests pass locally before submitting.

## Pre-commit hooks

The repository uses [pre-commit](https://pre-commit.com) to catch common issues before they reach review:

```bash
# Run hooks manually across all files (same as CI)
pre-commit run --all-files
```

Hooks check for: valid YAML/TOML, valid Python AST, merge conflict markers, accidentally large files, and leftover debug statements.

## Submitting changes

1. Fork the repository and create a feature branch
2. Make your changes with tests where applicable
3. Ensure all tests pass: `pytest tests/unit`
4. Open a pull request against `master` with a clear description

## Reporting bugs

Please use the [bug report template](https://github.com/tmancal74/quantarhei/issues/new?template=bug_report.md) and include a minimal reproducible example.

## Full contributing guide

For detailed guidelines on the development workflow, writing tests, and the project structure, see the [full contributing documentation](https://quantarhei.readthedocs.io/en/latest/contributing.html).

## Code of Conduct

This project follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this standard.
