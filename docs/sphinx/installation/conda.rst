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

Make sure that you really run Python installed with Anaconda. On some systems,
notably on some Linuxes, Python will already be installed when you first
log in. It might however by Python 2. Try the following command on command
line

.. code:: bash

    $ python --version
    
You should see something like this

.. code:: bash

    Python 3.6.4 :: Anaconda, Inc.
    
    
**OpenSuSE note**: When installing, Anaconda asks whether it can prepend its
$HOME/anaconda3/bin path to the $PATH variable in your `.bashrc` file. 
Only in this case you propertly get Anaconda Python called when typing `python`
on the command line. On (at least one version of) OpenSuSE, this leads to
X-Windows not loading after you re-log in. In this case, it is recommended not
to prepend the path to Anaconda Python in .bashrc, but only top do it by 
hand when you first open the terminal:

.. code:: bash

    $ export PATH=$HOME/anaconda3/bin:$PATH 


2. Installation using `conda` command
-------------------------------------

The second (and the last) step is to invoke the `conda` command and to install
Quantarhei package for your platform

.. code:: bash

    $ conda install -c tmancal74 quantarhei 
    
The `-c` option specifies the chanell from which Quantarhei is installed. The
value `tmancal74` is here to stay for a foreseable future.
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