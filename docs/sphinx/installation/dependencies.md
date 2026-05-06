(dependencies-label)=

# Quantarhei's Dependencies

All runtime dependencies are declared in [pyproject.toml] and installed
automatically by `pip install quantarhei`.

## Runtime Dependencies

- [numpy] — numerical arrays and linear algebra
- [scipy] — scientific algorithms
- [matplotlib] — plotting
- [terminaltables] — formatted console output
- [dill] — extended pickling for serialization
- [pyyaml] — YAML configuration file support
- [packaging] — version parsing utilities
- [gherkin-official] — Gherkin parser for acceptance tests

## Optional Dependencies

- [sympy] — symbolic mathematics (install with `pip install quantarhei[symbolic]`)

## Development Dependencies

Install with `pip install -e ".[dev]"`:

- [pytest] — unit test runner
- [pytest-cov] — test coverage
- [behave] — acceptance (BDD) tests
- [pre-commit] — Git hook management
- [ruff] — linter and formatter
- [mypy] — static type checker

[behave]: https://behave.readthedocs.io
[dill]: https://dill.readthedocs.io
[gherkin-official]: https://pypi.org/project/gherkin-official/
[matplotlib]: https://matplotlib.org
[mypy]: https://mypy-lang.org
[numpy]: https://numpy.org
[packaging]: https://packaging.pypa.io
[pre-commit]: https://pre-commit.com
[pyproject.toml]: https://github.com/tmancal74/quantarhei/blob/master/pyproject.toml
[pytest]: https://pytest.org
[pytest-cov]: https://pytest-cov.readthedocs.io
[pyyaml]: https://pyyaml.org
[ruff]: https://docs.astral.sh/ruff/
[scipy]: https://scipy.org
[sympy]: https://www.sympy.org
[terminaltables]: https://robpol86.github.io/terminaltables/
