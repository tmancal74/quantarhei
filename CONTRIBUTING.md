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

### VS Code and Codex

Open the repository root in VS Code and install the workspace's recommended
extensions. Run **Python: Select Interpreter** and select `.venv/bin/python`
(Windows: `.venv/Scripts/python.exe`). An existing interpreter selection takes
precedence over the workspace default. The current local environment can be
kept; for a new environment matching PR CI, use
`uv sync --extra dev --extra qutip --python 3.12`.

Workspace settings enable Ruff format-on-save and pytest in the Testing sidebar.
Ruff uses the environment's version when available; rules remain in
`pyproject.toml`. Mypy remains the project's type checker; Pylance provides
navigation and completion without a second set of type-checking diagnostics.
Automatic test discovery on save is disabled to avoid repeated imports of the
scientific stack; use **Test: Refresh Tests** after adding tests.

Use **Tasks: Run Task** for the existing Makefile tasks. For quick feedback
without environment synchronization or coverage reports:

```bash
MPLBACKEND=Agg uv run --no-sync pytest tests/unit/spectroscopy/twod_test.py -x
```

The **Python: TwoD calculator tests** launch profile runs focused calculator
tests under the debugger; set breakpoints in `twodcalculator.py` and press F5.
The Testing sidebar also supports debugging individual tests. **Python: Current
file** runs scripts from the repository root; examples requiring data files may
need a different working directory. Test debug profiles use the noninteractive
Matplotlib backend; current-file debugging keeps interactive plotting available.

Codex project instructions live in `AGENTS.md`. Start Codex in the repository
root and give each task a concrete behavior, relevant files, and acceptance
criteria. For example: "Investigate this TwoDResponseCalculator behavior, add a
focused numerical regression test, and run the affected spectroscopy tests."
Use planning for scientific or architectural changes and review numerical
assumptions and tolerances in the resulting diff. Keep model and personal
permission preferences in your Codex user settings.

Official references: [Codex IDE extension](https://learn.chatgpt.com/docs/codex/ide),
[AGENTS.md](https://learn.chatgpt.com/docs/agent-configuration/agents-md), and
[Ruff editor setup](https://docs.astral.sh/ruff/editors/setup/).

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

## Docstring style

All new and updated docstrings should follow **NumPy style**. This is the convention used by NumPy, SciPy, and the rest of the scientific Python ecosystem. It renders correctly with the `numpydoc` Sphinx extension already enabled in `docs/sphinx/conf.py`.

```python
def some_method(self, val, units="1/cm"):
    """Convert a value to internal units.

    Parameters
    ----------
    val : float or numpy.ndarray
        Value to convert.
    units : str, optional
        Unit string. Default is ``"1/cm"``.

    Returns
    -------
    float or numpy.ndarray
        Value in internal units.

    Raises
    ------
    UnitsError
        If ``units`` is not a recognized unit string.
    """
```

Key conventions:
- One-line summary on the opening `"""` line
- If a description follows, separate it from the summary with a blank line
- Section names (`Parameters`, `Returns`, `Raises`, `Notes`, `Examples`) are underlined with dashes of the same length
- Parameter types go after a colon; mark optional parameters with `, optional`
- Use double backticks for inline code inside docstrings

Full reference: [NumPy docstring guide](https://numpydoc.readthedocs.io/en/latest/format.html)

## Submitting changes

1. Fork the repository and create a feature branch
2. Make your changes with tests where applicable
3. Ensure all tests pass: `uv run pytest tests/unit`
4. Open a pull request against `master` with a clear description

Label pull requests with `bug`, `enhancement`, or `tests` where it applies.
These labels decide which section of the release notes a PR appears in (see
`.github/release.yml`). Unlabeled PRs go under "Other Changes".

## Releasing

Versions follow an odd/even scheme: `master` carries an odd development
version (for example `0.0.71`), and releases use the next even version
(`0.0.72`). To release `0.0.72`:

1. In a pull request, set `version = "0.0.72"` in `pyproject.toml`, run
   `uv lock`, and add a `## [0.0.72]` section to `CHANGELOG.md`. Merge it into
   `master`.
2. Check that `master` has the release version:
   `git show origin/master:pyproject.toml | grep '^version'`.
3. Publish a GitHub Release with tag `v0.0.72` on `master`, using generated
   notes:

   ```bash
   gh release create v0.0.72 --target master --title v0.0.72 --generate-notes
   ```

   In the web UI, use **Releases > Draft a new release**, create the tag
   `v0.0.72` on `master`, click **Generate release notes**, and review. You can
   save it as a draft first; nothing runs until you click **Publish release**.
4. Publishing starts `.github/workflows/publish-to-pypi.yml`. It checks the tag
   against `pyproject.toml`, builds the sdist and wheel, runs `twine check`,
   installs the wheel and compares `quantarhei.__version__` with the tag, then
   uploads to PyPI with trusted publishing. Confirm the new version on
   <https://pypi.org/project/quantarhei/>.
5. Bump `master` to the next development version (`0.0.73`) and run
   `uv lock`.

`CHANGELOG.md` is the authoritative, hand-written summary of each version. The
generated release notes list the merged PRs and serve as the GitHub Release
body.

Pushing a tag on its own does not publish anything. Pre-releases (the
**Set as a pre-release** box) are never uploaded to PyPI; to trial a build, run
the "Publish Python package to TestPyPI" workflow by hand from the **Actions**
tab.

The GitHub Release and its tag exist before the workflow runs. If a check
fails, nothing is uploaded to PyPI; the failed run prints the recovery steps.
Delete the release and tag (`gh release delete v0.0.72 --cleanup-tag --yes`),
fix `master`, and publish again. If only the upload failed for a transient
reason, use **Re-run failed jobs** on the workflow run.

## Reporting bugs

Please use the [bug report template](https://github.com/tmancal74/quantarhei/issues/new?template=bug_report.md) and include a minimal reproducible example.

## Full contributing guide

For detailed guidelines on the development workflow, writing tests, and the project structure, see the [full contributing documentation](https://quantarhei.readthedocs.io/en/latest/contributing.html).

## Code of Conduct

This project follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this standard.
