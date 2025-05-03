import contextlib
import time

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

import qutip
import quantarhei as qr


from qutip import (
    basis,
    brmesolve,
    destroy,
    expect,
    liouvillian,
    qeye,
    sigmax,
    sigmay,
    sigmaz,
    spost,
    spre,
    tensor,
    Qobj,
    coherent,
    coherent_dm,
)

from qutip.solver.heom import (
    BosonicBath,
    DrudeLorentzBath,
    DrudeLorentzPadeBath,
    HEOMSolver,
    HSolverDL,
)

from qutip.core.dimensions import (
    to_tensor_rep,
    from_tensor_rep,
)

@contextlib.contextmanager
def timer(label):
    """Simple utility for timing functions:

    with timer("name"):
        ... code to time ...
    """
    start = time.time()
    yield
    end = time.time()
    print(f"{label}: {end - start}")

from random import randint
import time



"""
fmo_ham = numpy.array([
            [410.0,  -87.7,    5.5,   -5.9,    6.7,  -13.7,   -9.9],
            [-87.7,  530.0,   30.8,    8.2,    0.7,   11.8,    4.3],
            [  5.5,   30.8,  210.0,  -53.5,   -2.2,   -9.6,    6.0],
            [ -5.9,    8.2,  -53.5,  320.0,  -70.7,  -17.0,  -63.3],
            [  6.7,    0.7,   -2.2,  -70.7,  480.0,   81.1,   -1.3],
            [-13.7,   11.8,   -9.6,  -17.0,   81.1,  630.0,   39.7],
            [ -9.9,    4.3,    6.0,  -63.3,   -1.3,   39.7,  440.0]
        ], dtype=qr.REAL)

"""


def run_simulation(Ham, sbi, heom_params, rhoi):

    timea = sbi.get_time_axis()

    id = int(time.time())

    # number of states 
    Ngr = 1
    Ndim = Ham.dim
    Nex = Ndim  - Ngr  # This has to be checked when we use aggregate

    with qr.energy_units("1/cm"):
        hdata = Ham.data

    #
    # Use Quantarhei data
    #
    Hsys = Qobj(hdata) # cm^-1
    T_qr = sbi.get_temperature() #300 # from cfce

    # eigen energies and eigen states
    eigenenergies, ekets = Hsys.eigenstates()


    NC = heom_params["hierarchy_level"] #4 # max hierarchy level
    Nk = heom_params["number_of_exponentials"] # 0 # bath exponential terms (0 = single correlation term) + we can add terminator

    Nt = timea.length 
    end_fs = timea.data[Nt-1]

    t_slices = Nt #1000 # at which time points to evaluate (how many)
    end_time_ps = end_fs / 1000 # 0.5

    xlims = [0, end_time_ps * 1000]   # in femto seconds
    ylims = [0, 1]
    tlist = np.linspace(0, end_time_ps * 1e-12 * 2*np.pi*3e10, t_slices) # 1e-12 s -> cm^-1

    # describing initial condition 
    """
    onesover7 = np.ones((7, 7))/7
    diagonly = False
    # ones in sites
    a = Qobj(onesover7)
    # into exc
    b = a.transform(ekets, False)
        # delete non dia
    c = None
    if (diagonly == True):
        c = Qobj(np.diag(b.diag()))
    else:
        # no delete of non dia - keep full
        c = b
    # back into site
    d = c.transform(ekets, True)

    rho0 = d
    """

    # initial excitation on site 1
    rho0 = Qobj(rhoi.data) #basis(Ndim, 1) * basis(Ndim, 1).dag()

    # describe what to put in graphs
    POPSONLY = True
    COHERENCESONLY = not POPSONLY

    """
    cfce = sbi.get_correlation_function((0,0))
    with qr.energy_units("1/cm"):
        lam_qr = cfce.get_reorganization_energy() # from cfce

    gamma_qr = 1.0/cfce.get_correlation_time() # from cfce
    print("Corelation function parameters: ", T_qr, lam_qr, gamma_qr)
    """

    Q_list = []
    baths = []

    Ltot = liouvillian(Hsys)

    for m in range(Nex):

        # projector
        Q = basis(Ndim, m+Ngr) * basis(Ndim, m+Ngr).dag()

        # FIXME: Get it from SBI here
        Qm = np.zeros((Ndim, Ndim), dtype=qr.REAL)
        Qm[m+Ngr, m+Ngr] = 1.0
        Qalt = Qobj(Qm)

        Q_list.append(Qalt)

        # bath correlation function 

        cfce = sbi.get_correlation_function((m,m))
        with qr.energy_units("1/cm"):
            lam_qr = cfce.get_reorganization_energy() # from cfce
        gamma_qr = 1.0/cfce.get_correlation_time() # from cfce

        T = T_qr * 0.6949 # K -> cm^-1
        gamma = gamma_qr * (1e15/ (2*np.pi*3e10))  # fs^-1 -> cm^-1
        lam = lam_qr # cm^-1

        baths.append(
            DrudeLorentzBath(
                Q, lam=lam, gamma=gamma, T=T, Nk=Nk,
                tag=str(m)
            )
        )
        _, terminator = baths[-1].terminator()
        Ltot += terminator


    # numerical options 
    options = {
        "nsteps": 5000,
        "store_states": True,
        "rtol": 1e-12,
        "atol": 1e-12,
        "min_step": 1e-18,
        "method": "vern9",
        "progress_bar": "enhanced",
    }


    with timer("RHS construction time"):
        HEOMMats = HEOMSolver(Ltot, baths, NC, options=options)

    with timer("ODE solver time"):
        output = HEOMMats.run(rho0, tlist)

    rho_t = qr.DensityMatrixEvolution(timeaxis=timea, rhoi=rhoi)
    for t in range(t_slices):
        rho_t.data[t,:,:] = output.states[t][:,:]

    return rho_t

    #
    # CONSTRUCTING EVOLUTION SUPEROPERATOR
    #
    evosuperops = []
    for t in range(t_slices):
        evosuperops.append(np.zeros((Ndim,Ndim,Ndim,Ndim), dtype=complex))

    for m in range(Ndim):
        for n in range(m+1):
            rho0_nm = basis(Ndim, n) * basis(Ndim, m).dag()
            print(rho0_nm)
            with timer("ODE solver time"):
                outputFMO_HEOM = HEOMMats.run(rho0_nm, tlist)

            for t in range(t_slices):
                evosuperops[t][n,m] += to_tensor_rep(outputFMO_HEOM.states[t])
                evosuperops[t][m,n] += to_tensor_rep(outputFMO_HEOM.states[t].dag())

    for n in range(Ndim):
        for t in range(t_slices):
            evosuperops[t][n,n] /= 2

    evosuperops_transformed = []
    for t in range(t_slices):
        evosuperops_transformed.append(np.zeros((Ndim,Ndim,Ndim,Ndim), dtype=complex))


    #
    #   TRANSFORMATIONS
    #

    # make transformation matrices
    vs = []
    for i in range(Ndim):
        vs.append(to_tensor_rep(ekets[i]))

    Smatrix = np.hstack(vs)
    SmatrixInv = np.linalg.inv(Smatrix)

    # example usage of transformation matrices and transform function
    """
    Hsysnp = to_tensor_rep(Hsys)
    print(Hsysnp)
    # second parameter is "Inverse", i.e. Inverse = False -> transform into said basis, Inverse = True -> transform from said basis
    HsysDiag = Hsys.transform(ekets, False).tidyup(1e-10)
    print(HsysDiag)
    ham_np = SmatrixInv @ Hsysnp @ Smatrix
    ham_np[ham_np < 1e-10] = 0
    print(ham_np)
    print(np.allclose(ham_np, to_tensor_rep(HsysDiag)))
    quit()
    """

    S = Smatrix
    S_inv = SmatrixInv

    for t in range(t_slices):
        print("t: ", t)
        dim = S.shape[0]
        """evosuperops_copy = evosuperops[t]
        evosuperops_reshaped = evosuperops_copy.reshape(dim**2, dim**2)
        S_kron = np.kron(S, S)
        S_inv_kron = np.kron(S_inv, S_inv)
        evosuperops_transformed[t] = S_kron @ evosuperops_reshaped @ S_inv_kron.T
        evosuperops_transformed[t] = evosuperops_transformed[t].reshape(dim, dim, dim, dim)
        """
        """
        for k_prime in range(dim):
            for l_prime in range(dim):
                for i_prime in range(dim):
                    for j_prime in range(dim):
                        for k in range(dim):
                            for l in range(dim):
                                for i in range(dim):
                                    for j in range(dim):
                                        evosuperops_transformed[t][k_prime, l_prime, i_prime, j_prime] += (
                                            np.conj(S[k, k_prime])
                                            * np.conj(S[l, l_prime])
                                            * evosuperops[t][k, l, i, j]
                                            * S[i, i_prime]
                                            * S[j, j_prime]
                                        )
        """
        evosuperops_transformed[t] = np.einsum(
            'km,ln,klij,ip,jq->mnpq',
            np.conj(S),  # Transformation for the first input index
            np.conj(S),  # Transformation for the second input index
            evosuperops[t], # Original evosuperops
            S,           # Transformation for the first output index
            S            # Transformation for the second output index
        )


    #
    #    PLOTTING
    #

    colors = ['r', 'g', 'b', 'y', 'c', 'm', 'k']

    if True:
        rho_final = []
        # Transformed
        fig, axes = plt.subplots(1, 1, figsize=(12, 8))
        fig2, axes2 = plt.subplots(1, 1, figsize=(12, 8))

        rho0numpy = SmatrixInv @ to_tensor_rep(rho0) @ Smatrix
        
        for t in range(t_slices):
            rho_f = np.zeros((Ndim,Ndim), dtype=complex)
            for i in range(Ndim):
                for j in range(Ndim):
                    for k in range(Ndim):
                        for l in range(Ndim):
                            rho_f[i,j] += evosuperops_transformed[t][k,l,i,j] * rho0numpy[k,l]

            rho_f_to_append = rho_f
            rho_final.append(rho_f_to_append)
        
        for n in range(Ndim):
            for m in range(Ndim):
                
                if POPSONLY == True and (n != m):
                    continue

                if COHERENCESONLY == True:
                    if (n == m):
                        continue

                    # only coherences with first site/exc)
                    if (n != 0):
                        continue
                    
                if True:
                    # sitebasis == True
                    # excbasis == False

                    """
                    # qutip calculation
                    Q = basis(7, n) * basis(7, m).dag()
                    vals = []
                    for t in range(t_slices):
                        # transform rho back into site basis
                        vals.append(expect(Qobj(rho_final[t]).transform(ekets, True), Q))
                    """

                    # numpy calculation
                    Q_np = to_tensor_rep(basis(Ndim, n) * basis(Ndim, m).dag())
                    vals = []
                    for t in range(t_slices):
                        # transform rho back into site basis
                        vals.append(np.trace(Smatrix @ rho_final[t] @ SmatrixInv @ Q_np))
                    
                    ls = "-"
                    axes.plot(
                        np.array(tlist) * 1e15 / 3e10 / 2 / np.pi,
                        np.real(vals),
                        label=f"{n + 1}{m+1}",
                        color=colors[m % len(colors)],
                        linestyle=ls,
                    )
                    axes.set_xlabel(r'$t$ (fs)', fontsize=30)
                    axes.set_ylabel(r"vals", fontsize=30)
                    axes.locator_params(axis='y', nbins=6)
                    axes.locator_params(axis='x', nbins=6)
                    axes.set_xlim(xlims)
                    axes.set_ylim(ylims)

                if True:
                    # sitebasis == False
                    # excbasis == True

                    """
                    # qutip calculation
                    Q = basis(7, n) * basis(7, m).dag()
                    vals = []
                    for t in range(t_slices):
                        vals.append(expect(Qobj(rho_final[t]), Q))
                    """
                    
                    # numpy calculation
                    Q_np = to_tensor_rep(basis(Ndim, n) * basis(Ndim, m).dag())
                    vals = []
                    for t in range(t_slices):
                        vals.append(np.trace(rho_final[t] @ Q_np))
                    
                    ls = "-"
                    axes2.plot(
                        np.array(tlist) * 1e15 / 3e10 / 2 / np.pi,
                        np.real(vals),
                        label=f"{n + 1}{m+1}",
                        color=colors[m % len(colors)],
                        linestyle=ls,
                    )
                    axes2.set_xlabel(r'$t$ (fs)', fontsize=30)
                    axes2.set_ylabel(r"vals", fontsize=30)
                    axes2.locator_params(axis='y', nbins=6)
                    axes2.locator_params(axis='x', nbins=6)
                    axes2.set_xlim(xlims)
                    axes2.set_ylim(ylims)

        axes.legend(loc=1)
        fig.savefig("fullevo_sitebasis_" + str(NC) + "_" + str(Nk)  + "_" + str(id) + ".png")
        axes2.legend(loc=1)
        fig2.savefig("fullevo_exctbasis_" + str(NC) + "_" + str(Nk)  + "_" + str(id) + ".png")

    evosuperops_transformed_secular = []
    for t in range(t_slices):
        evosuperops_transformed_secular.append(np.zeros((Ndim,Ndim,Ndim,Ndim), dtype=complex))

    for t in range(t_slices):
        for i in range(Ndim):
            for j in range(Ndim):
                for k in range(Ndim):
                    for l in range(Ndim):
                        if (k == l) and (i == j):
                            evosuperops_transformed_secular[t][k,l,i,j] += evosuperops_transformed[t][k,l,i,j]
                        if (k == i) and (l == j):
                            evosuperops_transformed_secular[t][k,l,i,j] += evosuperops_transformed[t][k,l,i,j]
                            if (k == l) and (i == j):
                                evosuperops_transformed_secular[t][k,l,i,j] -= evosuperops_transformed[t][k,l,i,j]

    if True:
        rho_final = []
        # Transformed SECULAR
        
        fig, axes = plt.subplots(1, 1, figsize=(12, 8))
        fig2, axes2 = plt.subplots(1, 1, figsize=(12, 8))

        rho0numpy = SmatrixInv @ to_tensor_rep(rho0) @ Smatrix

        for t in range(t_slices):
            rho_f = np.zeros((Ndim,Ndim), dtype=complex)
            for i in range(Ndim):
                for j in range(Ndim):
                    for k in range(Ndim):
                        for l in range(Ndim):
                            rho_f[i,j] += evosuperops_transformed_secular[t][k,l,i,j] * rho0numpy[k,l]
            rho_f_to_append = rho_f
            rho_final.append(rho_f_to_append)
        
        for n in range(Ndim):
            for m in range(Ndim):

                if POPSONLY == True and (n != m):
                    continue

                if COHERENCESONLY == True:
                    if (n == m):
                        continue

                    # only coherences with first site/exc)
                    if (n != 0):
                        continue
                
                if True:
                    # sitebasis == True
                    # excbasis == False

                    """
                    # qutip calculation
                    Q = basis(7, n) * basis(7, m).dag()
                    vals = []
                    for t in range(t_slices):
                        # transform rho back into site basis
                        vals.append(expect(Qobj(rho_final[t]).transform(ekets, True), Q))
                    """

                    # numpy calculation
                    Q_np = to_tensor_rep(basis(Ndim, n) * basis(Ndim, m).dag())
                    vals = []
                    for t in range(t_slices):
                        vals.append(np.trace(Smatrix @ rho_final[t] @ SmatrixInv @ Q_np))
                    
                    ls = "--"
                    axes.plot(
                        np.array(tlist) * 1e15 / 3e10 / 2 / np.pi,
                        np.real(vals),
                        label=f"{n + 1}{m+1}",
                        color=colors[m % len(colors)],
                        linestyle=ls,
                    )
                    axes.set_xlabel(r'$t$ (fs)', fontsize=30)
                    axes.set_ylabel(r"vals", fontsize=30)
                    axes.locator_params(axis='y', nbins=6)
                    axes.locator_params(axis='x', nbins=6)
                    axes.set_xlim(xlims)
                    axes.set_ylim(ylims)


                if True:
                    # sitebasis == False
                    # excbasis == True

                    """
                    # qutip calculation
                    Q = basis(7, n) * basis(7, m).dag()
                    vals = []
                    for t in range(t_slices):
                        vals.append(expect(Qobj(rho_final[t]), Q))
                    """
                    
                    # numpy calculation
                    Q_np = to_tensor_rep(basis(Ndim, n) * basis(Ndim, m).dag())
                    vals = []
                    for t in range(t_slices):
                        vals.append(np.trace(rho_final[t] @ Q_np))
                    
                    ls = "--"
                    axes2.plot(
                        np.array(tlist) * 1e15 / 3e10 / 2 / np.pi,
                        np.real(vals),
                        label=f"{n + 1}{m+1}",
                        color=colors[m % len(colors)],
                        linestyle=ls,
                    )
                    axes2.set_xlabel(r'$t$ (fs)', fontsize=30)
                    axes2.set_ylabel(r"vals", fontsize=30)
                    axes2.locator_params(axis='y', nbins=6)
                    axes2.locator_params(axis='x', nbins=6)
                    axes2.set_xlim(xlims)
                    axes2.set_ylim(ylims)

        
        axes.legend(loc=1)
        fig.savefig("secexc_sitebasis_" + str(NC) + "_" + str(Nk)  + "_" + str(id) + ".png")
        axes2.legend(loc=1)
        fig2.savefig("secexc_exctbasis_" + str(NC) + "_" + str(Nk)  + "_" + str(id) + ".png")
                    

    if True:
        # Normal HEOM
        fig, axes = plt.subplots(1, 1, figsize=(12, 8))
        fig2, axes2 = plt.subplots(1, 1, figsize=(12, 8))
        with timer("ODE solver time"):
            outputFMO_HEOM = HEOMMats.run(rho0, tlist)    
        
        for n in range(Ndim):
            for m in range(Ndim):
                
                if POPSONLY == True and (n != m):
                    continue

                if COHERENCESONLY == True:
                    if (n == m):
                        continue

                    # only coherences with first site/exc)
                    if (n != 0):
                        continue
                
                if True:
                    # sitebasis == True
                    # excbasis == False

                    # qutip calculation
                    Q = basis(Ndim, n) * basis(Ndim, m).dag()
                    vals = []
                    for t in range(t_slices):
                        vals.append(expect(outputFMO_HEOM.states[t], Q))
                    
                    axes.plot(
                        np.array(tlist) * 1e15 / 3e10 / 2 / np.pi,
                        np.real(vals),
                        label=f"{m + 1}{n+1}",
                        color=colors[n % len(colors)],
                        linestyle=ls,
                    )
                    axes.set_xlabel(r'$t$ (fs)', fontsize=30)
                    axes.set_ylabel(r"vals", fontsize=30)
                    axes.locator_params(axis='y', nbins=6)
                    axes.locator_params(axis='x', nbins=6)
                    axes.set_xlim(xlims)
                    axes.set_ylim(ylims)

                if True:
                    # sitebasis == False
                    # excbasis == True

                    # qutip calculation
                    Q = basis(Ndim, n) * basis(Ndim, m).dag()
                    vals = []
                    for t in range(t_slices):
                        # transform into exc basis
                        vals.append(expect(outputFMO_HEOM.states[t].transform(ekets, False), Q))
                    
                    axes2.plot(
                        np.array(tlist) * 1e15 / 3e10 / 2 / np.pi,
                        np.real(vals),
                        label=f"{m + 1}{n+1}",
                        color=colors[n % len(colors)],
                        linestyle=ls,
                    )
                    axes2.set_xlabel(r'$t$ (fs)', fontsize=30)
                    axes2.set_ylabel(r"vals", fontsize=30)
                    axes2.locator_params(axis='y', nbins=6)
                    axes2.locator_params(axis='x', nbins=6)
                    axes2.set_xlim(xlims)
                    axes2.set_ylim(ylims)
        
        axes.legend(loc=1)
        fig.savefig("heom_sitebasis_" + str(NC) + "_" + str(Nk)  + "_" + str(id) + ".png")
        axes2.legend(loc=1)
        fig2.savefig("heom_exctbasis_" + str(NC) + "_" + str(Nk)  + "_" + str(id) + ".png")