# Developing Quantarhei

Quantarhei is a scientific Python package for molecular open quantum systems
and spectroscopy. Read CONTRIBUTING.md and pyproject.toml for the development
workflow; .github/workflows/python-package.yml defines CI checks.

## Environment and commands

- Use the repository's uv-managed .venv and editable package installation.
  Setup: `uv sync --extra dev`; add `--extra qutip` for CI-equivalent unit
  dependencies. Python 3.12 matches PR CI; retain Python 3.10 compatibility.
- For an already prepared environment, use `uv run --no-sync` to avoid
  unexpectedly changing dependencies or an in-progress uv.lock.
- Focused spectroscopy tests:
  `MPLBACKEND=Agg uv run --no-sync pytest tests/unit/spectroscopy/twod_test.py`.
- Unit suite: `MPLBACKEND=Agg uv run --no-sync pytest tests/unit`.
- Check changed Python files with `uv run --no-sync ruff check <paths>` and
  `uv run --no-sync ruff format --check <paths>`.
- Type checks: `uv run --no-sync mypy quantarhei/`.
- Doctest and acceptance commands are in CONTRIBUTING.md and Makefile.
  `make lint` invokes pre-commit hooks that can modify files.

## Code and verification

- Core axes, units, and basis management are in quantarhei/core; molecular
  construction in builders; quantum dynamics in qm; spectra in spectroscopy.
- Tests live in tests/unit (both test_*.py and *_test.py are used), acceptance
  tests in tests/behave/features, and executable examples in examples.
- Follow existing NumPy-style docstrings and the Ruff configuration.
  Format changed files; avoid unrelated repository-wide cleanup.
- Preserve physical units, basis conventions, array shapes, complex dtypes,
  normalization, and public API compatibility. Explain scientific assumptions.
- For numerical changes, add focused regression checks with justified
  tolerances. Do not regenerate reference data merely to make tests pass.
- Start with affected tests, then broaden to related core/builders/qm tests
  when changes cross those boundaries. Report exactly what was verified.
- Preserve unrelated work and lockfile changes. Do not commit generated
  spectra, coverage output, or machine-specific settings.
- PR titles must start with `#<issue-number> `, as required by CI.
