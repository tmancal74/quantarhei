# Installation from PyPI

The simplest way to install Quantarhei is via `pip` from [PyPI].

## Requirements

Quantarhei requires **Python 3.10 or later**. Check your version:

```bash
$ python --version
Python 3.12.3
```

## Installation

Install Quantarhei and all its dependencies:

```bash
$ pip install quantarhei
```

If you use [uv] (recommended for faster, reproducible installs):

```bash
$ uv pip install quantarhei
```

## Testing the Installation

Verify the installation:

```bash
$ python -c "import quantarhei as qr; print(qr.Manager().version)"
0.0.69
```

[pypi]: https://pypi.org/project/quantarhei/
[uv]: https://github.com/astral-sh/uv
