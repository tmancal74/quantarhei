# Installation from Source Code

To install from source, clone the repository and install in editable mode.

## 1. Python

Make sure you have **Python 3.10 or later** installed:

```bash
$ python --version
Python 3.12.3
```

## 2. Clone the Repository

```bash
$ git clone https://github.com/tmancal74/quantarhei.git
$ cd quantarhei
```

Alternatively, download a zip from the [releases page] and unzip it.

## 3. Install in Editable Mode

Using `pip`:

```bash
$ pip install -e .
```

Using [uv] (faster):

```bash
$ uv pip install -e .
```

This installs all required runtime dependencies automatically.

## 4. Install Development Dependencies

To run tests and use the linter/formatter:

```bash
$ pip install -e ".[dev]"
$ pre-commit install
```

## 5. Testing the Installation

Run the unit test suite:

```bash
$ pytest tests/unit
```

Verify the package:

```bash
$ python -c "import quantarhei as qr; print(qr.Manager().version)"
0.0.69
```

[releases page]: https://github.com/tmancal74/quantarhei/releases
[uv]: https://github.com/astral-sh/uv
