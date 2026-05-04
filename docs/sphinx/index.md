---
substitutions:
  Qrhei: '**Quantarhei**'
---

% Quantarhei documentation master file, created by
% sphinx-quickstart on Mon Oct  3 10:48:47 2016.
% You can adapt this file completely to your liking, but it should at least
% contain the root `toctree` directive.

# Welcome to Quantarhei's Documentation!

[Qrhei] is a Molecular Open Quantum Systems Simulator written predominantly
in Python. Its name is derived from the famous aphorism "Panta rhei" of the
Greek philosopher Heraclitus of Ephesus. "Panta rhei" means "Everything flows"
or "Everything is in flux", which is quite fitting when you change Panta into
Quanta.

In "Quantarhei" the last four letter ("rhei") should be written in Greek,
i.e. (using LateX convention) "\\rho \\epsilon \\iota".

This page is meant to be the complete source of documentation
for the [Qrhei] package. As the documentation project progresses, we will
describe [Qrhei]'s main features and
the philosophy behind them, together with a complete description of its
functionality and content.

Large part of this documentation is based directly
on [Qrhei]'s source code. Inspecting the source code should be an
important part of learning to use [Qrhei]. Not only that you can learn
to use [Qrhei]
better by inspecting the source code, but, in a better case, you learn something
usefull from how open quantum systems' problems are solved in [Qrhei].
In a worse case, you will be motivated to fix [Qrhei]'s deficiencies.

There are two types of deficiencies that one can expect in [Qrhei] - a mild one:
things work well but programming style is terrible, or things are not
implemented in a general enough manner. In this case you are most welcome
to fix the code. Make sure that your improvement is equipped with tests and that
the coded passes all existing automatic tests.

A more serious defficiency is when you find that something is really
implemented wrongly. The best approach then is to write an alternative
test code which demonstrates the errors. Submit this to the maintainers
and when they agree that the error is real, go ahead to fix it (or get
it fixed by the maintainers).

# Current status of {{ Qrhei }}

[Qrhei]'s source code is available on [Github]. Binary packages are
published on [PyPI] and installable with `pip install quantarhei`.
The test suite runs on Python 3.10, 3.11, and 3.12 via GitHub Actions,
with coverage reported on [Codecov].
Documentation is hosted on [Readthedocs] and built directly from the
source code.
[Qrhei] is Open Source software published under the MIT license.

# Detailed Documentation

Dive into [Qrhei] documentation below:

```{toctree}
:maxdepth: 2

Get started <getstarted>
Installation <installation>
Examples <examples>
API Reference <classes>
Advanced API <advclasses>
Internals <internals>
Contributing <contribute>
```

# Indices and tables

- {ref}`genindex`
- {ref}`modindex`
- {ref}`search`

[codecov]: https://codecov.io/gh/tmancal74/quantarhei
[github]: https://github.com/tmancal74/quantarhei
[pypi]: https://pypi.org/project/quantarhei/
[qrhei]: http://github.com/tmancal74/quantarhei
[readthedocs]: https://quantarhei.readthedocs.io/en/latest/
