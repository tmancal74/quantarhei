.. _dependencies-label:

Quantarhei's Dependencies
=========================

All runtime dependencies are declared in `pyproject.toml`_ and installed
automatically by `pip install quantarhei`.

Runtime Dependencies
--------------------

- `numpy`_ — numerical arrays and linear algebra
- `scipy`_ — scientific algorithms
- `matplotlib`_ — plotting
- `terminaltables`_ — formatted console output
- `dill`_ — extended pickling for serialization
- `pyyaml`_ — YAML configuration file support
- `packaging`_ — version parsing utilities
- `gherkin-official`_ — Gherkin parser for acceptance tests

Optional Dependencies
---------------------

- `sympy`_ — symbolic mathematics (install with ``pip install quantarhei[symbolic]``)

Development Dependencies
------------------------

Install with ``pip install -e ".[dev]"``:

- `pytest`_ — unit test runner
- `pytest-cov`_ — test coverage
- `behave`_ — acceptance (BDD) tests
- `pre-commit`_ — Git hook management
- `ruff`_ — linter and formatter
- `mypy`_ — static type checker


.. _`pyproject.toml`: https://github.com/tmancal74/quantarhei/blob/master/pyproject.toml
.. _`numpy`: https://numpy.org
.. _`scipy`: https://scipy.org
.. _`matplotlib`: https://matplotlib.org
.. _`terminaltables`: https://robpol86.github.io/terminaltables/
.. _`dill`: https://dill.readthedocs.io
.. _`pyyaml`: https://pyyaml.org
.. _`packaging`: https://packaging.pypa.io
.. _`gherkin-official`: https://pypi.org/project/gherkin-official/
.. _`sympy`: https://www.sympy.org
.. _`pytest`: https://pytest.org
.. _`pytest-cov`: https://pytest-cov.readthedocs.io
.. _`behave`: https://behave.readthedocs.io
.. _`pre-commit`: https://pre-commit.com
.. _`ruff`: https://docs.astral.sh/ruff/
.. _`mypy`: https://mypy-lang.org
