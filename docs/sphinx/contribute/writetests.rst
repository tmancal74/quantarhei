.. _write-tests:

How to Write Tests for Quantarhei
=================================

Every contribution must include appropriate tests. Quantarhei uses three
levels of testing: unit tests, doctests, and acceptance (BDD) tests.

Running Tests
-------------

Run the unit test suite:

.. code:: bash

    $ pytest tests/unit

Run with coverage:

.. code:: bash

    $ pytest tests/unit --cov=quantarhei --cov-report=term-missing

Run all checks (lint + format + type check + tests) via pre-commit:

.. code:: bash

    $ pre-commit run --all-files

Writing Unit Tests
------------------

Unit tests live in ``tests/unit/``. Use `pytest`_ conventions — test
files are named ``test_*.py`` and test functions start with ``test_``:

.. code:: python

    def test_molecule_energy_levels():
        import quantarhei as qr
        m = qr.Molecule([0.0, 1.0])
        assert m.Nel == 2

Place your test file in the subdirectory that mirrors the module being
tested, e.g. ``tests/unit/builders/test_molecules.py`` for code in
``quantarhei/builders/molecules.py``.

Writing Doctests
----------------

Doctests are examples embedded directly in docstrings. They are run as
part of the Sphinx build and the pytest suite. Follow NumPy docstring
conventions:

.. code:: python

    def get_energy(self, state):
        """Return the energy of the given state.

        Parameters
        ----------
        state : int
            Index of the electronic state.

        Returns
        -------
        float
            Energy in internal units.

        Examples
        --------
        >>> import quantarhei as qr
        >>> m = qr.Molecule([0.0, 1.0])
        >>> m.get_energy(1)
        1.0

        """

Writing Acceptance Tests
------------------------

Acceptance tests use `Behave`_ (BDD) and live in
``quantarhei/testing/resources/behave/``. Features are written in
Gherkin syntax (``*.feature`` files) with corresponding step
definitions in Python.

.. code:: gherkin

    Feature: Molecule creation
      Scenario: Create a two-level molecule
        Given I create a molecule with energies [0.0, 1.0]
        Then the molecule has 2 electronic levels

Run acceptance tests:

.. code:: bash

    $ behave quantarhei/testing/resources/behave/


.. _`pytest`: https://pytest.org
.. _`Behave`: https://behave.readthedocs.io
