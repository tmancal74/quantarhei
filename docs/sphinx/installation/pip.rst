Installation from Pypi 
======================

If you have standard Python 3 installed, you most likely have the `pip`
command installed with it. The `pip` command allows you to manage local
installations of packages and connects you to the Python Package Index (PIPY),
where you can find large collection of packages including Quantarhei and 
all its dependencies.

Note: Installing Python is easy on all standard platforms. We 
recommend that you install scientific python distribution from Anaconda_. With
this Python distribution you can follow the installation description below,
but it may be more convenient to install with the `conda` command (next section).

1. `pip` command
----------------

Check for the presence of the `pip` command type on the command line:

.. code:: bash

    $ pip 
    
If you get an output like this

.. code:: bash

    Usage:   
      pip <command> [options]
    
    Commands:
      install                     Install packages.
      download                    Download packages.
      ...
      
      
etc., you have the `pip` command installed.

2. Installation
---------------

Now you have everything need to install Quantarhei. Just type

.. code:: bash

    $ pip install quantarhei
    
on the command line, and `pip` will install all Quantarhei dependencies (which
may take quite a while if you only have a bare Python 3 installed.)

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
 