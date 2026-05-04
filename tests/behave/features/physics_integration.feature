Feature: Physics integration tests pinning numerical outcomes of core simulations
    These scenarios encode the expected physical behaviour of the core simulation
    workflows.  They are intentionally coarse-grained: they do not test internal
    implementation details but instead assert that each end-to-end workflow
    produces results that are consistent with known analytic or physical constraints.

    # -----------------------------------------------------------------------
    # A. Lindblad two-level system – coherence decay matches prescribed rate
    # -----------------------------------------------------------------------

    Scenario: Lindblad coherence decays exponentially at the prescribed rate
        Given a two-level molecule with transition energy 12000 invcm
          And a Lindblad relaxation rate of 1 per 200 fs
         When I propagate the density matrix for 1000 fs with time step 1 fs
         Then the off-diagonal element decays as exp(-t/200) within 1 percent

    # -----------------------------------------------------------------------
    # B. Heterodimer absorption – peak positions match exciton eigenvalues
    # -----------------------------------------------------------------------

    Scenario: Heterodimer absorption peaks coincide with exciton eigenvalues
        Given a heterodimer with site energies 12000 and 12300 invcm and coupling 100 invcm
          And an overdamped Brownian bath with reorganisation energy 20 invcm and correlation time 100 fs at 300 K
         When I calculate the absorption spectrum
         Then there are two peaks whose positions match the exciton eigenvalues within 50 invcm

    # -----------------------------------------------------------------------
    # C. Redfield homodimer – detailed balance at long times
    # -----------------------------------------------------------------------

    Scenario: Redfield homodimer populations satisfy detailed balance at long times
        Given a homodimer-2-env test aggregate with coupling 100 invcm
         When I propagate with standard Redfield theory for 5000 fs with time step 1 fs
         Then the ratio of forward to backward rates satisfies exp(-dE/kT) within 5 percent

    # -----------------------------------------------------------------------
    # D. Trace conservation under Redfield propagation
    # -----------------------------------------------------------------------

    Scenario: Density matrix trace is conserved during Redfield evolution
        Given a heterodimer with site energies 12000 and 12300 invcm and coupling 100 invcm
          And an overdamped Brownian bath with reorganisation energy 20 invcm and correlation time 100 fs at 300 K
         When I propagate with standard Redfield theory for 1000 fs with time step 1 fs
         Then the trace of the density matrix equals 1 at every time step within 1e-6

    # -----------------------------------------------------------------------
    # E. Förster rate matches analytic formula in the weak-coupling limit
    # -----------------------------------------------------------------------

    Scenario: Foerster transfer rate matches the analytic weak-coupling formula
        Given a heterodimer with site energies 12000 and 12500 invcm and coupling 10 invcm
          And an overdamped Brownian bath with reorganisation energy 20 invcm and correlation time 100 fs at 300 K
         When I calculate the Foerster rate matrix
         Then the forward rate matches 2*pi*J^2*spectral_overlap within 10 percent

    # -----------------------------------------------------------------------
    # F. 2DES dimer – cross-peaks are present and rephasing differs from
    #    non-rephasing
    # -----------------------------------------------------------------------

    Scenario: 2DES dimer has cross-peaks in the total spectrum
        Given a heterodimer with site energies 12000 and 12300 invcm and coupling 100 invcm
          And a Lindblad relaxation rate of 1 per 200 fs
         When I calculate the 2D electronic spectrum at population time 0 fs
         Then cross-peaks are present in the total 2D spectrum
