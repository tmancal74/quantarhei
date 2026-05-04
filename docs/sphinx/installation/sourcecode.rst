Installation from Source Code
=============================

To install from source, clone the repository and install in editable mode.

1. Python
---------

Make sure you have **Python 3.10 or later** installed:

.. code:: bash

    $ python --version
    Python 3.12.3

2. Clone the Repository
------------------------

.. code:: bash

    $ git clone https://github.com/tmancal74/quantarhei.git
    $ cd quantarhei

Alternatively, download a zip from the `releases page`_ and unzip it.

3. Install in Editable Mode
----------------------------

Using `pip`:

.. code:: bash

    $ pip install -e .

Using `uv`_ (faster):

.. code:: bash

    $ uv pip install -e .

This installs all required runtime dependencies automatically.

4. Install Development Dependencies
-------------------------------------

To run tests and use the linter/formatter:

.. code:: bash

    $ pip install -e ".[dev]"
    $ pre-commit install

5. Testing the Installation
----------------------------

Run the unit test suite:

.. code:: bash

    $ pytest tests/unit

Verify the package:

.. code:: bash

    $ python -c "import quantarhei as qr; print(qr.Manager().version)"
    0.0.69


.. _`releases page`: https://github.com/tmancal74/quantarhei/releases
.. _`uv`: https://github.com/astral-sh/uv
