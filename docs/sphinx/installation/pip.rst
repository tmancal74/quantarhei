Installation from PyPI
======================

The simplest way to install Quantarhei is via `pip` from `PyPI`_.

Requirements
------------

Quantarhei requires **Python 3.10 or later**. Check your version:

.. code:: bash

    $ python --version
    Python 3.12.3

Installation
------------

Install Quantarhei and all its dependencies:

.. code:: bash

    $ pip install quantarhei

If you use `uv`_ (recommended for faster, reproducible installs):

.. code:: bash

    $ uv pip install quantarhei

Testing the Installation
------------------------

Verify the installation:

.. code:: bash

    $ python -c "import quantarhei as qr; print(qr.Manager().version)"
    0.0.69


.. _`PyPI`: https://pypi.org/project/quantarhei/
.. _`uv`: https://github.com/astral-sh/uv
