# Laboratory setup

`LabSetup` owns the experimental pulse definitions; `LabField` is the
per-pulse view used when a calculation needs a time-domain field.  Pulse
parameters and carrier frequencies are stored in Quantarhei internal units.
An `energy_units` context therefore controls only how numeric input and output
arguments are interpreted, not the physical field stored by the object.

## Pulse and field evaluation

For a finite pulse, use the explicit evaluators on `LabField`:

```python
field = lab.get_labfield(0)
envelope = field.envelope_at(times)
positive = field.field_p_at(times, rwa_frequency=omega_rwa)
negative = field.field_m_at(times, rwa_frequency=omega_rwa)
real_field = field.real_field_at(times)
```

The positive-frequency analytic field is

\[
E^{(+)}(t; \Omega) = A(t)\exp[-i(\omega-\Omega)(t-t_c)+i\phi],
\]

where the configured phase is defined at the pulse centre `t_c`.  The
negative-frequency component is its complex conjugate and `real_field_at()`
returns their half-sum.  All evaluators accept either a scalar or an array;
numeric pulses evaluate to zero outside their sampled support.

Use `derivative_at()` when a calculation requires a time derivative:

```python
d_envelope = field.derivative_at(times, component="envelope")
d_positive = field.derivative_at(times)  # component="positive"
d_negative = field.derivative_at(times, component="negative")
d_real = field.derivative_at(times, component="real")
```

Time-defined Gaussian envelopes are differentiated analytically.  Numeric
pulses, and pulses originally specified in frequency space, are differentiated
on their derived time grid and interpolated to the requested times.  The
positive-frequency derivative includes both the envelope derivative and the
carrier term.  `LabSetup.get_field_derivative()` returns the corresponding
sum over its pulses.  The historical `LabField.get_field_derivative()` remains
as a deprecated compatibility wrapper.

Delta pulses are area-defined objects for the impulsive workflow.  They do
not have a pointwise envelope, field, or derivative; request their area for an
impulsive calculation instead.

```{eval-rst}
.. automodule:: quantarhei.spectroscopy.labsetup
    :members:
```
