.. _how-to-contribute:

Process of Contributing Code to Quantarhei
==========================================

Thank you for considering a contribution to Quantarhei. The process is
straightforward:

1. **Open an issue** on the `issue tracker`_ describing the bug or feature.
   For non-trivial changes, discuss the approach with the maintainer before
   writing code.

2. **Fork and branch** — see :ref:`Using the GitHub Repository <using-the-github-repository>`
   for the mechanics. Work on a dedicated branch, never on ``master``.

3. **Write tests first** — every change must come with tests. See
   :ref:`write-tests` for guidance on unit tests, doctests, and acceptance
   tests.

4. **Keep changes focused** — one logical change per pull request. Large
   changes are hard to review and slow to merge.

5. **Pass CI** — your PR must pass the full CI pipeline (lint, format,
   type check, unit tests on Python 3.10/3.11/3.12) before it can be
   merged.

6. **Commit message format** — every commit must start with the GitHub
   issue number:

   .. code:: bash

       #123 fix: correct energy calculation in aggregate

7. **Open a draft PR** early — you can open a draft pull request before
   the work is finished to get early feedback.

Code Style
----------

Quantarhei uses `ruff`_ for linting and formatting. Run before committing:

.. code:: bash

    $ ruff check quantarhei/ --fix
    $ ruff format quantarhei/

Type annotations are encouraged for new code. `mypy`_ is run in CI:

.. code:: bash

    $ mypy quantarhei/


.. _`issue tracker`: https://github.com/tmancal74/quantarhei/issues
.. _`ruff`: https://docs.astral.sh/ruff/
.. _`mypy`: https://mypy-lang.org
