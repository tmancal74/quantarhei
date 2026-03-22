# Contributing to Quantarhei

Thank you for your interest in contributing to Quantarhei!

## Quick setup

```bash
git clone https://github.com/tmancal74/quantarhei
cd quantarhei
pip install -e .
pip install -r requirements_devel.txt
```

## Running tests

```bash
# Unit tests (recommended for quick feedback)
pytest tests/unit

# Full test suite (as run in CI)
paver test
```

Tests are run against Python 3.10, 3.11, and 3.12 in CI. Please verify locally on at least one of these versions before submitting.

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
