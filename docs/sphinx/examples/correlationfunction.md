# How to Use `CorrelationFunction`

`CorrelationFunction` describes the bath-induced fluctuations of a molecular
transition frequency. It requires a `TimeAxis` and a parameter dictionary that
specifies the functional form and its physical parameters.

## Creating a correlation function

```python
from quantarhei import CorrelationFunction, TimeAxis, energy_units

# define the time axis: start, number of steps, step size (in fs)
ta = TimeAxis(0.0, 1000, 1.0)

temperature = 300.0  # Kelvin
params = {
    "ftype": "OverdampedBrownian",
    "reorg": 20.0,    # reorganization energy in 1/cm
    "cortime": 100.0, # correlation time in fs
    "T": temperature,
    "matsubara": 20,
}

# supply energy parameters in wavenumbers
with energy_units("1/cm"):
    cf = CorrelationFunction(ta, params)

# inspect the reorganization energy
print("Reorganization energy:", cf.reorganization_energy)
```

## Accessing the spectral density

```python
from quantarhei import SpectralDensity, TimeAxis, energy_units

ta = TimeAxis(0.0, 1000, 1.0)
params = {
    "ftype": "OverdampedBrownian",
    "reorg": 20.0,
    "cortime": 100.0,
    "T": 300.0,
    "matsubara": 20,
}

with energy_units("1/cm"):
    sd = SpectralDensity(ta, params)

print(sd)
```

## Defining a correlation function through M(t)

Legacy line-shape models sometimes specify a normalized relaxation function
`M(t)` rather than a spectral density.  Such a model can be supplied through
the standard parameter dictionary with `ftype="M-defined"`.  The `M` array
has to be real, sampled on the supplied time axis, and normalized to
`M[0] == 1`.

```python
import numpy
from quantarhei import CorrelationFunction, TimeAxis, energy_units

ta = TimeAxis(0.0, 5000, 1.0)
tau = 130.0
m_values = numpy.exp(-(ta.data/tau)**2)

params = {
    "ftype": "M-defined",
    "M": m_values,
    "reorg": 140.0,       # 1/cm
    "T": 300.0,           # K
    "cutoff-time": 650.0, # fs, optional
}

with energy_units("1/cm"):
    cf = CorrelationFunction(ta, params)
```

Internally Quantarhei constructs the spectral density according to

```text
J(omega) = 2*reorg*omega*integral(M(t)*cos(omega*t), t=0..infinity)
```

and converts it to a numerical correlation function through the existing
value-defined mechanism.  The time axis should extend far enough that `M(t)`
has decayed at its upper boundary.  Insufficient time range causes truncation
and Fourier-transform artifacts.
