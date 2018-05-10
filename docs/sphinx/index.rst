.. Quantarhei documentation master file, created by
   sphinx-quickstart on Mon Oct  3 10:48:47 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

|Build Status| |DocBuild Status| |Coverage|

Welcome to Quantarhei's Documentation!
======================================

|Qrhei|_ is a Molecular Open Quantum Systems Simulator written predominantly
in Python. Its name is derived from the famous aphorism "Panta rhei" of the
Greek philosopher Heraclitus of Ephesus. "Panta rhei" means "Everything flows"
or "Everything is in flux" which is quite fitting when you change Panta into
Quanta.

In "Quantarhei" the last four letter ("rhei") should be written in Greek,
i.e. (using LateX convention) "\\rho \\epsilon \\iota".

This page should one day be the complete source of documentation
for |Qrhei|_ package. We describe |Qrhei|_'s main features and 
the philosophy behind them. 

Large part of this documentation is based directly
on the source code. Inspecting the source code should be an important part
of learning to use |Qrhei|_. Not only that one can learn to use |Qrhei|_
better by inspecting the source code. In a better case you learn something
usefull from how open quantum systems' problems are solved in |Qrhei|_,
in a worse case, you will be motivated to fix |Qrhei|_'s deficiencies. 

There two types of deficiencies that one can expect in |Qrhei|_ - a mild one:
things work well but programming style is terrible, or things are not
implemented in a general enough manner. In this case you are most welcome
to fix the code. Make sure that your improvement passes all automatic tests.
Here we provide a detailed description of the classes provided by |Qrhei|_.

A more serious defficiency is when you find that something is really
implemented wrongly. The best approach then is to write an alternative
test code which demonstrates the errors. Submit this to the maintainers
and when they agree that the error is real, go ahead to fix it (or get
it fixed by the maintainers).


Dive into |Qrhei|_ documentation below:

.. toctree::
   :maxdepth: 2

   getstarted
   installation
   examples
   classes
   advclasses
   mngrclasses
   internals
   contribute


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |DocBuild Status| image:: https://readthedocs.org/projects/quantarhei/badge/?version=latest
   :target: http://quantarhei.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
   
.. |Build Status| image:: https://travis-ci.com/tmancal74/quantarhei.svg?branch=master
   :target: https://travis-ci.com/tmancal74/quantarhei
   :alt: Build Status
 
.. |Coverage| image:: https://img.shields.io/codecov/c/github/tmancal74/quantarhei.svg
   :target: https://codecov.io/gh/tmancal74/quantarhei
   
.. |Qrhei| replace:: **Quantarhei**
.. _Qrhei: http://github.com/tmancal74/quantarhei
