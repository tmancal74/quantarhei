Using the GitHub Repository
===========================

Quantarhei is developed on `GitHub`_. This page covers the mechanics of
forking, cloning, and submitting pull requests.

Forking and Cloning
-------------------

1. **Fork** the repository on GitHub — click the Fork button on the
   `project page`_. This creates your own copy at
   ``https://github.com/YOUR-USERNAME/quantarhei``.

2. **Clone** your fork locally:

   .. code:: bash

       $ git clone https://github.com/YOUR-USERNAME/quantarhei.git
       $ cd quantarhei

3. **Add the upstream remote** so you can pull in future changes:

   .. code:: bash

       $ git remote add upstream https://github.com/tmancal74/quantarhei.git

Setting Up for Development
--------------------------

Install in editable mode with development dependencies:

.. code:: bash

    $ pip install -e ".[dev]"
    $ pre-commit install

The `pre-commit install` step sets up Git hooks that automatically run
the linter, formatter, and type checker before each commit.

Configuring Git
---------------

Make sure Git knows who you are:

.. code:: bash

    $ git config --global user.name "Your Name"
    $ git config --global user.email you@example.com

Keeping Your Fork Up to Date
-----------------------------

Fetch and merge upstream changes into your local ``master``:

.. code:: bash

    $ git fetch upstream
    $ git checkout master
    $ git merge upstream/master

If you are working on a feature branch, rebase it on top of the updated master:

.. code:: bash

    $ git checkout my-feature
    $ git rebase master

Making Changes
--------------

Work on a dedicated branch — never commit directly to ``master``:

.. code:: bash

    $ git checkout -b my-feature

After making changes, run the test suite and linter before committing:

.. code:: bash

    $ pytest tests/unit
    $ ruff check quantarhei/ --fix
    $ ruff format quantarhei/

Commit your changes:

.. code:: bash

    $ git add <changed files>
    $ git commit -m "#<issue-number> short description"

.. note::

    Commit messages must start with the GitHub issue number, e.g.
    ``#123 fix: correct energy calculation in aggregate``.

Push your branch to your fork:

.. code:: bash

    $ git push -u origin my-feature

Pull Requests
-------------

Open a pull request from your branch to ``tmancal74/quantarhei master``
on the `GitHub website`_. Keep pull requests small and focused — one
logical change per PR makes review faster.

Your PR will be reviewed and merged by the maintainer once it passes CI
and review. See :ref:`how-to-contribute` for contribution guidelines and
:ref:`write-tests` for the testing requirements.


.. _`GitHub`: https://github.com/tmancal74/quantarhei
.. _`project page`: https://github.com/tmancal74/quantarhei
.. _`GitHub website`: https://github.com
