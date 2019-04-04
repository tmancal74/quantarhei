# -*- coding: utf-8 -*-

import quantarhei as qr
from quantarhei import LabSetup
from quantarhei.utils.vectors import X 
import quantarhei.functions as func

use_vib = True
vib_loc = "up"

make_movie = False

HR = 0.05
omega = 500.0
dE = 500.0
JJ = 30.0


def run(omega, HR, dE, JJ, vib_loc="up", use_vib=True, make_movie=False):

    E0 = 10000.0
    dip1 = [1.5, 0.0, 0.0]
    dip2 = [-1.0, -1.0, 0.0]
    width = 200.0
    rate = 1.0/200.0
    
    # dimer of molecules
    with qr.energy_units("1/cm"):
        mol1 = qr.Molecule([0.0, E0])
        mol2 = qr.Molecule([0.0, E0+dE])
        
        mod1 = qr.Mode(omega)
        mod2 = qr.Mode(omega)
    
    mol1.set_transition_width((0,1), qr.convert(width, "1/cm", "int"))
    mol1.set_dipole(0,1, dip1)
    
    mol2.set_transition_width((0,1), qr.convert(width, "1/cm", "int"))
    mol2.set_dipole(0,1, dip2)
    
    agg = qr.Aggregate([mol1, mol2])
    
    with qr.energy_units("1/cm"):
        agg.set_resonance_coupling(0,1,JJ)
    
    agg.save("agg.qrp")
    
    if use_vib:
    
        if vib_loc == "down":
            set_vib = [True, False]
        elif vib_loc == "up":
            set_vib = [False, True]
        elif vib_loc == "both":
            set_vib = [True, True]
        else:
            raise Exception("Unknown location of the vibrations")
            
        if set_vib[0]:
            mol1.add_Mode(mod1)
            mod1.set_nmax(0, 3)
            mod1.set_nmax(1, 3)
            mod1.set_HR(1, HR)
        
        if set_vib[1]:
            mol2.add_Mode(mod2)
            mod2.set_nmax(0, 3)
            mod2.set_nmax(1, 3)
            mod2.set_HR(1, HR)
    
        
    agg_el = qr.load_parcel("agg.qrp")
    
    agg.save("agg_vib.qrp")
    agg3 = qr.load_parcel("agg_vib.qrp")
    
    agg.build(mult=1)
    agg_el.build(mult=1)
    
    HH = agg.get_Hamiltonian()
    He = agg_el.get_Hamiltonian()
    
    #
    # Laboratory setup
    #
    lab = LabSetup()
    lab.set_polarizations(pulse_polarizations=[X,X,X], detection_polarization=X)
    
    time2 = qr.TimeAxis(0.0, 100, 10.0)
    
    cont = qr.TwoDSpectrumContainer(t2axis=time2)
    #
    # spectra will be indexed by the times in the time axis `time2`
    #
    cont.use_indexing_type(time2)
    
    #
    # We define two-time axes, which will be FFTed and will define the omega_1 and
    # omega_3 axes of the 2D spectrum
    #
    t1_N_steps = 100
    t1_time_step = 10.0
    t3_N_steps = 100
    t3_time_step = 10.0
    t1axis = qr.TimeAxis(0.0, t1_N_steps, t1_time_step)
    t3axis = qr.TimeAxis(0.0, t3_N_steps, t3_time_step)
    
    #
    # This calculator calculated 2D spectra from the effective width defined above
    #
    msc = qr.MockTwoDSpectrumCalculator(t1axis, time2, t3axis)
    msc.bootstrap(rwa=qr.convert(E0+(dE/2.0),"1/cm","int"), 
                  all_positive=False, shape="Gaussian")
    
    operators = []
    rates = []
    with qr.eigenbasis_of(He):
        operators.append(qr.qm.ProjectionOperator(1, 2, dim=He.dim))
    rates.append(rate)
    
    #
    # System-bath interaction including vibrational states
    #
    sbi = qr.qm.SystemBathInteraction(sys_operators=operators,
                                      rates=rates)
    sbi.set_system(agg)
    
    LF = qr.qm.ElectronicLindbladForm(HH, sbi, as_operators=True)
    
    #
    # Pure dephasing
    #
    p_deph = qr.qm.ElectronicPureDephasing(agg, dtype="Gaussian")
    
    
    eUt = qr.qm.EvolutionSuperOperator(time2, HH, relt=LF, pdeph=p_deph,
                                       mode="all")
    eUt.set_dense_dt(100)
    eUt.calculate(show_progress=True)
    
    agg3.build(mult=2)
    agg3.diagonalize()
    
    
    print("Dims: ", eUt.dim, HH.dim)
    for t2 in time2.data:
        print("t2 =", t2)
        twod = msc.calculate_one_system(t2, agg3, HH, eUt, lab)
        cont.set_spectrum(twod)
    
    if make_movie:
        with qr.energy_units("1/cm"):
            cont.make_movie("mov.mp4")
        
    #
    # Window function for subsequenty FFT
    #
    window = func.Tukey(time2, r=0.3, sym=False)
    
    #
    # FFT with the window function
    #
    # Specify REPH, NONR or `total` to get different types of spectra
    #
    print("\nCalculating FFT of the 2D maps")
    fcont = cont.fft(window=window, dtype="total") #, dpart="real", offset=0.0)
    
    show_omega = omega
    
    #
    # Have a look which frequencies we actually have
    #
    Ndat = len(fcont.axis.data)
    print("\nNumber of frequency points:", Ndat)
    print("In 1/cm they are:")
    with qr.energy_units("1/cm"):
        for k_i in range(Ndat):
            print(k_i, fcont.axis.data[k_i])
    
    with qr.frequency_units("1/cm"):
        sp1, show_Npoint1 = fcont.get_nearest(show_omega)
        sp2, show_Npoint2 = fcont.get_nearest(-show_omega)
    
    units = "1/cm"
    with qr.energy_units(units):
        print("\nPlotting spectrum at frequency:", 
              fcont.axis.data[show_Npoint1], units)
        #sp.plot(Npos_contours=10, window=[9500,11500, 9500, 11500],
        #          stype="total", spart="abs")
        sp1.plot(Npos_contours=10, 
                stype="total", spart="abs")   
        fftfile = "twod_fft_map_1.png"
        sp1.savefig(fftfile)#os.path.join(pre_out, fftfile))
        print("... saved into: ", fftfile)
        print("\nPlotting spectrum at frequency:", 
              fcont.axis.data[show_Npoint2], units)
        #sp.plot(Npos_contours=10, window=[9500,11500, 9500, 11500],
        #          stype="total", spart="abs")
        sp2.plot(Npos_contours=10, 
                stype="total", spart="abs")   
        fftfile = "twod_fft_map_2.png"
        sp2.savefig(fftfile)#os.path.join(pre_out, fftfile))
        print("... saved into: ", fftfile) 

run(omega, HR, dE, JJ, vib_loc,use_vib, make_movie)