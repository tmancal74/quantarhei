Installation from Source Code
=============================

To install from the source code, please follow the instructions below. 

1. Python
---------

Make sure you have |Python3|_ installed (versions 3.4, 3.5 and 3.6 are supported). Installing Python is easy on all standard platforms. We 
recommend that you install scientific python distribution from Anaconda_.

2. Clone or Download the Package
--------------------------------
 
Clone the repository on |github|_ or download a source bundle from |github_rel|_. 

Cloning from Github
~~~~~~~~~~~~~~~~~~~

To clone the Quantarhei repository, you need Git to be installed on your computer. If you have it, just type

.. code:: bash

    $ git clone https://github.com/tmancal74/quantarhei.git

Downloading the Package
~~~~~~~~~~~~~~~~~~~~~~~

Alternatively, if you do not want to install Git, you can download the package
in zip format. Got to |github_rel|_, download the zip file and unzip it in
a suitable directory.  

  
3. Dependencies
---------------
  
You need to first have the dependencies installed. Check the package dependencies in section |dependencies-label| and install them.

This is usually as easy as typing

.. code:: bash

    $ pip install numpy

to install the numpy package, which is one of the most important dependencies
of Quantarhei. Do this similarly for all required packages.


4. Installation
---------------

With all dependencies installed there are two ways how you can install the
Quantarhei package. 

Installation with `pip`
~~~~~~~~~~~~~~~~~~~~~~~

First you can create a source code distribution by runnning
Quantarhei setup script. Enter the directory containing the package (it will
contain a Python script called `setup.py`) and type

.. code:: bash

    $ python setup.py sdist
    
This will create a directory called `dist` (if it did not already exist), and
it will put in a distribution file, e.g. `quantarhei-0.0.36.tar.gz`, into it.
This file can be then installed by the `pip` command

.. code:: bash

    $ pip install dist/quantarhei-0.0.36.tar.gz
    
Now Quantarhei is installed on your system. You can check the installation
by `pip`. Just type

.. code:: bash

    $ pip show quantarhei
    
and you get the following message

.. code:: bash

    Name: quantarhei
    Version: 0.0.36
    Summary: Quantarhei: Open Quantum System Theory for Molecular Systems
    Home-page: https://github.com/tmancal74/quantarhei
    Author: Tomas Mancal
    Author-email: mancal@karlov.mff.cuni.cz
    License: MIT
    Location: /Users/tomas/anaconda3/lib/python3.6/site-packages
    Requires: h5py, numpy, terminaltables, scipy, matplotlib
    
You can always remove the installation by typing

.. code:: bash

    $ pip uninstall quantarhei
    
    
Installation by `setup.py`
~~~~~~~~~~~~~~~~~~~~~~~~~~

Alternatively, you can do without `pip` in the installation of Quantarhei
(except the you have to install its dependencies which would be quite tedious
to install manually). In the Quantarhei directory type
    
.. code:: bash

    $ python setup.py install
    
and all is done. The message from `pip` will now read

.. code:: bash

    $ pip show quantarhei
    Name: quantarhei
    Version: 0.0.36
    Summary: Quantarhei: Open Quantum System Theory for Molecular Systems
    Home-page: https://github.com/tmancal74/quantarhei
    Author: Tomas Mancal
    Author-email: mancal@karlov.mff.cuni.cz
    License: MIT
    Location: /Users/tomas/anaconda3/lib/python3.6/site-packages/quantarhei-0.0.36-py3.6.egg
    Requires: numpy, scipy, matplotlib, h5py, terminaltables

and the only difference seems to be that the package is install in form of 
an `egg`. 

Removal of the installation works the same way as above.


.. |github| replace:: github.org

.. _github: http://github.com/tmancal74/quantarhei

.. |github_rel| replace:: Quantarhei release page on github.org

.. _github_rel: https://github.com/tmancal74/quantarhei/releases


.. |dependencies-label| replace:: :ref:`dependencies-label`

.. _Python3: http://www.python.org

.. |Python3| replace:: Python 3

.. _Anaconda: http://www.anaconda.com

