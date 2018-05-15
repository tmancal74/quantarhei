Absorption Spectrum
===================

Here we describe a group of classes which enable calculation, handling and
analysis of absorption spectra. All spectroscopy methods are represented
by *calculator*, *spectrum* and *spectrum container* classes. Often the main
class, in this classes the class `AbsSpectrum` inherits from a chain of classes
which implement some of its more basic functionality. In case of absorption
spectra, we have the following classes that represent the spectrum

.. toctree::
   :maxdepth: 1
   
   abss/absbase
   abss/abs
   
   
To calculate absorption spectrum

.. toctree::
   :maxdepth: 1
   
   abss/abscalculator
   
For storage, plotting and further analysis of groups of spectra, Quantarhei
offers a container class

.. toctree::
   :maxdepth: 1
   
   abss/abscontainer
   
   