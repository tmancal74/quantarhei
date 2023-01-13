.. Quantarhei documentation master file, created by
   sphinx-quickstart on Mon Oct  3 10:48:47 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



Welcome to Quantarhei's Documentation!
======================================

|Qrhei|_ is a Molecular Open Quantum Systems Simulator written predominantly
in Python. Its name is derived from the famous aphorism "Panta rhei" of the
Greek philosopher Heraclitus of Ephesus. "Panta rhei" means "Everything flows"
or "Everything is in flux", which is quite fitting when you change Panta into
Quanta.

In "Quantarhei" the last four letter ("rhei") should be written in Greek,
i.e. (using LateX convention) "\\rho \\epsilon \\iota".

This page is meant be the complete source of documentation
for |Qrhei|_ package. We describe |Qrhei|_'s main features and 
the philosophy behind them, together with a complete description of its
functionality and content.

Large part of this documentation is based directly
on |Qrhei|_'s source code. Inspecting the source code should be an 
important part of learning to use |Qrhei|_. Not only that you can learn 
to use |Qrhei|_
better by inspecting the source code, but, in a better case, you learn something
usefull from how open quantum systems' problems are solved in |Qrhei|_.
In a worse case, you will be motivated to fix |Qrhei|_'s deficiencies. 

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


Current status of |Qrhei| 
=========================

|Qrhei|_'s source code is available from `Github`_, source and binary bundles
for easy installation can be found on `Pypi`_ and `Anaconda Cloud`_. Latest
builds are tested on `Travis CI`_ and the st coverage is measured 
by `coverage` package and displayed on `Codecov`_
Documentation for |Qrhei|_ can be found on `Readthedocs`_. Also the documentaion
is build from the source code. Besides that we have some
interactive examples installed on `Crosscompute`_. 
|Qrhei|_ is an Open Source software published under the MIT license, which is
short and very non-restrictive.
At present we are in the alpha stage of the development. Our interim goal is 
version 0.1.0, for which we are developing a definition (see our `Github Wiki
page`_)
We expect our package to run on more than the most current version of Python,
but at the moment we do not explicitely test it.
When installing from Anaconda Clound (recommended) the package supports
all major platforms.


Detailed Documentation
======================

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


Parts of Documentation Considered Complete
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1
   
   classes/fulldoc


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

 
.. |Qrhei| replace:: **Quantarhei**
.. _Qrhei: http://github.com/tmancal74/quantarhei

.. _`Github`: http://github.com/tmancal74/quantarhei
.. _`Travis CI`: https://travis-ci.com/tmancal74/quantarhei
.. _`Readthedocs`: http://quantarhei.readthedocs.io/en/latest/
.. _`Pypi`: https://pypi.org/project/quantarhei/
.. _`Anaconda Cloud`: https://anaconda.org/tmancal74/quantarhei
.. _`Crosscompute`:
.. _`Codecov`: https://codecov.io/gh/tmancal74/quantarhei
.. _`Github Wiki page`: https://github.com/tmancal74/quantarhei/wiki/Quantarhei-0.1