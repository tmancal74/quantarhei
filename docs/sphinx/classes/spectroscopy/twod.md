# Coherent 2D Spectrum

Quantarhei provides a set of classes for calculating and representing
coherent two-dimensional (2D) electronic and vibrational spectra. The main
workflow begins with `TwoDResponseCalculator`, which computes the nonlinear
optical response of a molecular system and stores the results in a
`TwoDResponseContainer`. The container can then be converted into a
`TwoDSpectrumContainer` (where the pulse-field convolution is applied), from
which individual `TwoDSpectrum` objects are retrieved for plotting and
analysis.

```{eval-rst}
.. automodule:: quantarhei.spectroscopy.twodspect
    :members:
```

```{eval-rst}
.. automodule:: quantarhei.spectroscopy.twodcontainer
    :members:
```

```{eval-rst}
.. automodule:: quantarhei.spectroscopy.twodcalculator
    :members:
```
