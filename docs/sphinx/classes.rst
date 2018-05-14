Quantarhei User Level Classes
=============================

This is a list of classes that are available through the import of
the `Quantarhei` package. You can import them all as

.. code:: python

    from quantarhei import *
    
although the recommended way of importing and using Quantarhei is to rename
it, e.g. as

.. code:: python

    import quantarhei as qr
    
after which you can use all of its classes via a prefex `qr`. This way you will
not inferere with other packages which might accidentally use the same class
names. Use the `TimeAxis` class of Quantarhei e.g. as

.. code:: python

    timea = qr.TimeAxis(0.0, 1000, 1.0)
    
We call these classes **user level classes** because they are intended for usage
in user scripts and programs. There are also **advanced level classes**, which
should be used by experienced, advanced users to enable them e.g. to changed details
of implemented simulation methods. Advanced level classes are one level below
the user level classes and their usage is for instance as follows

.. code:: python

    import quantarhei as qr
    
    sbi = qr.qm.SystemBathInteraction


The other classes in Quantarhei are considered
**expert level classes** and they should be tinkered with on when you develop
Quantarhei (away, this is just a recommendation).


Classes Representing General Concepts
-------------------------------------

.. toctree::
   :maxdepth: 2

   classes/value
   classes/time
   classes/frequency
   classes/dfunction


Molecular Systems
-----------------

.. toctree::
   :maxdepth: 2
   
   classes/molecule
   classes/modes
   classes/aggregates   
   
Molecular Environment
---------------------

.. toctree::
   :maxdepth: 2
   
   classes/correlationfunction
      
   
Quantum Mechanics
-----------------

.. toctree::
   :maxdepth: 2   

   classes/statevector


   
 