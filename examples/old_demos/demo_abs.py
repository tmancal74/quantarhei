#
###############################################################################
#
# ABSORPTION SPECTRUM CALCULATION
#
# Two ways to calculate absorption spectrum in Quantarhei
#
# 1) Secular calculation with realistic lineshapes
# 2) Secular calculation with an effective Gaussian or Lorenzian lineshape
#
###############################################################################
#

import quantarhei as qr

#
# Define spectroscopic lab settings 
#
from quantarhei.utils.vectors import X, Y, Z

lab = qr.LabSetup()
lab.set_polarizations(pulse_polarizations=[X,X,X], detection_polarization=X)

###############################################################################
#
# Calculation with a realistic lineshape using bath correlation function
#
###############################################################################

#
# Set up molecules composing an aggregate
#
with qr.energy_units("1/cm"):
    # molecule 1
    # molecular states 
    m1 = qr.Molecule([0.0, 10000.0])
    # transition dipole moment
    m1.set_dipole(0,1,[1.0, 0.0, 0.0])
    
    # molecule 2
    m2 = qr.Molecule([0.0, 11000.0])
    m2.set_dipole(0,1,[1.0, 0.0, 0.0])
    
#
# Create aggregate of the two molecules
#
agg = qr.Aggregate(molecules=[m1, m2])

# 
# Define time span of the calculation
#
time = qr.TimeAxis(0.0, 10000, 1.0)

#
# Define bath correlation function 
#
cpar = dict(ftype="OverdampedBrownian", cortime=30, reorg=200, T=300)
with qr.energy_units("1/cm"):
    cf = qr.CorrelationFunction(time, cpar)

#
# Set the correlation function to the transitions on the molecules
#
m1.set_transition_environment((0,1), cf)
m2.set_transition_environment((0,1), cf)

#
# Build the aggregate
#
agg.build()

#
# Calculate absorption spectrum
#
asc = qr.AbsSpectrumCalculator(time, system=agg)
asc.bootstrap(rwa=qr.convert(10050,"1/cm", "int"), lab=lab)
abs1 = asc.calculate(raw=False)
abs1.normalize2()

#
# Plot the absorption spectrum
#
with qr.energy_units("1/cm"):
    abs1.plot(axis=[9000,12000,0.0,1.01], color="b", show=False)
    

###############################################################################
#
# Calculation using effective Gaussian lineshape
#
###############################################################################

#
# Set up molecules composing an aggregate
#
with qr.energy_units("1/cm"):
    m3 = qr.Molecule([0.0, 9900.0])
    m3.set_dipole(0,1,[1.0, 0.0, 0.0])
    m4 = qr.Molecule([0.0, 10900.0])
    m4.set_dipole(0,1,[1.0, 0.0, 0.0])
   
#
# The molecules get a transition width which effectively describes optical
# dephasing.
#
m3.set_transition_width((0,1), qr.convert(280,"1/cm", "int"))
m4.set_transition_width((0,1), qr.convert(280,"1/cm", "int"))
m3.set_transition_dephasing((0,1), 1.0/200.0)
m4.set_transition_dephasing((0,1), 1.0/200.0)

#
# Create the aggregate of the two molecules
# 
agg2 = qr.Aggregate(molecules=[m3, m4])

#
# Build the aggregate
#
agg2.build()


#
#  Absorption Spectrum by pathway method
#
mac = qr.MockAbsSpectrumCalculator(time, system=agg2)
mac.bootstrap(rwa=qr.convert(10000.0,"1/cm","int"), 
              shape="Gaussian", lab=lab)

abs1 = mac.calculate(raw=False)
abs1.normalize2()

#
# Plot the absorption spectrum
#
with qr.energy_units("1/cm"):
    abs1.plot(axis=[9000,12000,0.0,1.01], color="r", show=False)


mac.bootstrap(rwa=qr.convert(10000.0,"1/cm","int"), 
              shape="Voigt", lab=lab)
abs2 = mac.calculate()
abs2.normalize2()    
with qr.energy_units("1/cm"):
    abs2.plot(axis=[9000,12000,0.0,1.01], color="g")

#
# EOF
#