Installation with Anaconda
==========================

The recommended way of using Quantarhei is with the scientific Python 
distribution from Anaconda_. With Anaconda_ you get a well maintained
multiplatform Python distribution, where the core Quantarhei dependencies,
namely `scipy` and `numpy` packages are already installed. The `conda`
command which comes with Ananconda distribution helps you to install
Quantarhei smoothly on the main computer platforms.


1. Anaconda Scientific Python
-----------------------------

As a first step, install the Anaconda distribution of Python. Download and 
installation instructions for all major platforms can be found on the
Anaconda_ website.


2. Installation using `conda` command
-------------------------------------

The second (and the last) step is to invoke the `conda` command and to install
Quantarhei package for your platform

.. code:: bash

    $ conda install quantarhei
    
The `conda` command will install all the dependencies of Quantarhei.


3. Testing the Installation
---------------------------

To test that Quantarhei is installed on your system, you can run the Python
interpretter

.. code:: bash

    $ python
    
This starts the Python console

.. code:: bash

    Python 3.6.2 |Anaconda custom (64-bit)| (default, Sep 21 2017, 18:29:43) 
    [GCC 4.2.1 Compatible Clang 4.0.1 (tags/RELEASE_401/final)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>>
    
Type the following line of code into Python console

.. code:: python

    >>> import quantarhei as qr
    >>> qr.Manager().version
    
This will result in the console printing your version of Quantarhei, e.g.

.. code:: python

    '0.0.36'
    
Leave the Python console by typing

.. code:: python

    >>> quit()
    
    
You have just installed Quantarhei package successfully.
     


 .. _Anaconda: http://www.anaconda.com