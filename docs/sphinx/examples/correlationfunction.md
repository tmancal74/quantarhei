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
