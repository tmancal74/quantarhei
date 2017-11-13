QUANTArhei: Open Quantum System Theory for Molecular Systems
============================================================

Quantarhei is a Molecular Open Quantum Systems Simulator written predominantly
in Python. Its name is derived from the famous aphorism "Panta rhei" of the
Greek philosopher Heraclitus of Ephesus. "Panta rhei" means "Everything flows"
or "Everything is in flux" which is quite fitting when you change Panta into
Quanta.

In "Quantarhei" the last four letter ("rhei") should be written in Greek,
i.e. (using LateX convention) "\\rho \\epsilon \\iota". 

----

Quantarhei is in flux, but it already provides helper classes to define
molecules, their aggregates and their interaction with external environment.
It can calculate absorption spectra of individual molecules and their
aggregates and excitation energy transfer dynamics using various types
of Redfield and Foerster theories.

Quantarhei provides Python code (optimized with Numpy) for all its implemented
methods and theories, and allows extensions and replacements of the reference
Python code with optimised routines written in C, Fortran or other lower level
languages.

In the first development stage, we concentrate on bringing to you tools
to quickly build essential components of a quantum mechanical simulation,
such as Hamiltonian, relaxation tensors, various initial
conditions for density matrix etc.

Quantarhei is at its experimental stage. 
Current version is 0.0.31

Quantarhei is available in source form on GitHub and from PyPI for installation
with the pip command.


Acknowledgements
================

The work on Quantarhei is supported by

|NFN|_

.. |NFN| replace:: **Neuron Fund for Support of Science**
.. _NFN: http://www.nfneuron.cz

through the Impuls grant in physics 2014 (2015-2017)

and

|GACR|_

.. |GACR| replace:: **Czech Science Foundation (GACR)**
.. _GACR: http://www.gacr.cz
                                               

through grants: 14-25752S (2014-2016) and 17-22160S (2017- )



New in version 0.0.31
=====================

For users:

- Arbitrary time independent Lindblad form 
- quantarhei.wizard module which contains IPython magic commands and some helpful Python console commands
- Simulation templates which can be fetched into IPython notebooks or console by %template  magic command (IPython) or fetch_template (console and IPython)
- Part of the test suit available for installed Quantarhei package
- Some small improvements and bug fixes

For developers:

- Makefile is back in the package root directory
- examples directory depleted in favor of quantarhei/wizard/examples directory
- New tests under quantarhei/tests directory (mostly unit tests which contain plots)
- pytest required to run newtests with matplotlib plots
 
