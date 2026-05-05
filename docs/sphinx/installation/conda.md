# Installation with Anaconda

Quantarhei can be installed into an [Anaconda] or [Miniconda] environment.
This is a good option if you already use conda for managing your scientific
Python stack.

## 1. Create or activate an environment

Create a dedicated environment with Python 3.10 or later:

```bash
$ conda create -n qrhei python=3.12
$ conda activate qrhei
```

## 2. Install Quantarhei

Install from PyPI using `pip` inside the active conda environment:

```bash
$ pip install quantarhei
```

:::{note}
Quantarhei is not currently distributed via a conda channel.
Use `pip` inside your conda environment as shown above.
:::

## 3. Testing the Installation

```bash
$ python -c "import quantarhei as qr; print(qr.Manager().version)"
0.0.69
```

[anaconda]: https://www.anaconda.com
[miniconda]: https://docs.anaconda.com/miniconda/
