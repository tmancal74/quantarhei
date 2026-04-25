# Contributing to Quantarhei

Thank you for your interest in contributing to Quantarhei!

## Quick setup

```bash
git clone https://github.com/tmancal74/quantarhei
cd quantarhei
pip install -e .
pip install -r requirements_devel.txt
```

Optionally, activate the commit-msg hook that enforces issue numbers at the start of every commit message (e.g. `#123 fix: description`). The hook script lives in `.githooks/` in the repo, but git does not activate it automatically — each developer runs this once per clone to opt in:

```bash
git config core.hooksPath .githooks
```

Before pushing, make sure the code is formatted and lint-clean:

```bash
ruff format quantarhei/
ruff check quantarhei/ --fix
```

CI enforces both — a push with unformatted or unlinted code will fail the lint job.

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

## Developer tools

The repository uses three tools to enforce code quality. All three are installed by `pip install -r requirements_devel.txt`.

### pre-commit (optional)

[pre-commit](https://pre-commit.com) can run all checks automatically on every commit. It is **optional** — CI is the authoritative gate, not pre-commit.

**Why you might skip it:** the formatter and linter hooks auto-fix files but then block the commit (because the files changed), so you have to `git add` the fixes and commit again. Many developers find this annoying, especially on quick fixup commits. The recommended alternative is to configure format-on-save in your editor (the [Ruff VS Code extension](https://marketplace.visualstudio.com/items?itemName=charliermarsh.ruff) does this) and run `ruff format` + `ruff check --fix` manually before pushing.

If you do want hooks, install them once per clone:

```bash
pre-commit install
```

Or run the full suite manually at any time:

```bash
pre-commit run --all-files
```

### ruff

[ruff](https://docs.astral.sh/ruff/) is the linter and formatter. It runs automatically via pre-commit, but you can also run it directly:

```bash
# Check for violations
ruff check quantarhei/

# Auto-fix violations (safe fixes only)
ruff check quantarhei/ --fix
```

Configuration lives in `[tool.ruff]` in `pyproject.toml`.

### mypy

[mypy](https://mypy.readthedocs.io) checks static types. It runs automatically via pre-commit, but you can also run it directly:

```bash
mypy quantarhei/
```

Configuration lives in `[tool.mypy]` in `pyproject.toml`. The entire package is annotated and `mypy quantarhei/` should exit with `Success: no issues found`.

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
