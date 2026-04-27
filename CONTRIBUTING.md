# Contributing to Quantarhei

Thank you for your interest in contributing to Quantarhei!

## Quick setup

The recommended way to set up a development environment is with [uv](https://docs.astral.sh/uv/), which manages Python versions, virtual environments, and dependencies in one tool.

```bash
git clone https://github.com/tmancal74/quantarhei
cd quantarhei
uv sync --extra dev
```

That's it — `uv sync` creates a virtual environment, installs all runtime and development dependencies from `uv.lock`, and installs quantarhei in editable mode. uv picks the highest Python version available that satisfies `requires-python = ">=3.10"`.

To run any command in the environment:

```bash
uv run pytest tests/unit
uv run mypy quantarhei/
```

### Pinning a specific Python version

If you want the virtualenv to use a specific Python version (e.g. 3.12 while your system default is 3.13):

```bash
uv sync --extra dev --python 3.12
```

uv will download and manage that Python version automatically if it is not already installed. The `.venv` created in the project directory will use it. Subsequent `uv run` calls use whatever Python is in `.venv`, so no further flags are needed.

### Alternative: plain pip

If you prefer not to use uv:

```bash
pip install -e ".[dev]"
```

Optionally, activate the commit-msg hook that enforces issue numbers at the start of every commit message (e.g. `#123 fix: description`). The hook script lives in `.githooks/` in the repo, but git does not activate it automatically — each developer runs this once per clone to opt in:

```bash
git config core.hooksPath .githooks
```

Before pushing, make sure the code is formatted and lint-clean:

```bash
uv run ruff format .
uv run ruff check quantarhei/ --fix
```

CI enforces both — a push with unformatted or unlinted code will fail the lint job.

## Running tests

```bash
# Unit tests — fast feedback, run before every push
uv run pytest tests/unit

# Doc tests
uv run pytest --doctest-modules \
  quantarhei/core quantarhei/builders quantarhei/qm/corfunctions \
  quantarhei/spectroscopy quantarhei/qm/liouvillespace quantarhei/functions \
  quantarhei/qm/hilbertspace quantarhei/qm/propagators

# Behave acceptance tests
uv run coverage run -m behave tests/behave/features
```

CI runs all three suites on Python 3.10, 3.11, and 3.12. Please verify unit tests pass locally before submitting.

## Developer tools

All developer tools are installed as part of the `dev` extra (`uv sync --extra dev`).

### pre-commit (optional)

[pre-commit](https://pre-commit.com) can run all checks automatically on every commit. It is **optional** — CI is the authoritative gate, not pre-commit.

**Why you might skip it:** the formatter and linter hooks auto-fix files but then block the commit (because the files changed), so you have to `git add` the fixes and commit again. Many developers find this annoying, especially on quick fixup commits. The recommended alternative is to configure format-on-save in your editor (the [Ruff VS Code extension](https://marketplace.visualstudio.com/items?itemName=charliermarsh.ruff) does this) and run `ruff format` + `ruff check --fix` manually before pushing.

If you do want hooks, install them once per clone:

```bash
uv run pre-commit install
```

Or run the full suite manually at any time:

```bash
uv run pre-commit run --all-files
```

### ruff

[ruff](https://docs.astral.sh/ruff/) is the linter and formatter. It runs automatically via pre-commit, but you can also run it directly:

```bash
# Check for violations
uv run ruff check quantarhei/

# Auto-fix violations (safe fixes only)
uv run ruff check quantarhei/ --fix
```

Configuration lives in `[tool.ruff]` in `pyproject.toml`.

### mypy

[mypy](https://mypy.readthedocs.io) checks static types. It runs automatically via pre-commit, but you can also run it directly:

```bash
uv run mypy quantarhei/
```

Configuration lives in `[tool.mypy]` in `pyproject.toml`. The entire package is annotated and `mypy quantarhei/` should exit with `Success: no issues found`.

## Submitting changes

1. Fork the repository and create a feature branch
2. Make your changes with tests where applicable
3. Ensure all tests pass: `uv run pytest tests/unit`
4. Open a pull request against `master` with a clear description

## Reporting bugs

Please use the [bug report template](https://github.com/tmancal74/quantarhei/issues/new?template=bug_report.md) and include a minimal reproducible example.

## Full contributing guide

For detailed guidelines on the development workflow, writing tests, and the project structure, see the [full contributing documentation](https://quantarhei.readthedocs.io/en/latest/contributing.html).

## Code of Conduct

This project follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this standard.
