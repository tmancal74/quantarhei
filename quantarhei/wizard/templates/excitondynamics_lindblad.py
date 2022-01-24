"""

Simulation definition file for Quantarhei package. You can run the simulation
in the terminal using Quantarhei command line tool ``qrhei``. Alternatively you
can instantiate the class and run it in python/IPython console or in Jupyter 
notebook.


Usage in terminal
-----------------

Make sure you have ``quantarhei`` package installed, e.g. by typing

    > pip show quantarhei
    
Then type in the terminal

    > qrhei --simulation filename
  
Console/Notebook usage
----------------------

Load the definition file, instantiate it and run it as:

    >>> ms = mySimulation()
    >>> ms.run()

    
"""

from quantarhei.wizard.simulations import ExcitonDynamics

class mySimulation(ExcitonDynamics):
    """Example file for ExcitonDynamics simulation
    
    Edit this file to perform your own simulation within the predefined set
    of tasks. You can add molecules and interactions between them, specify 
    interaction with environment and the initial condition for the propagation.
    You can choose from several theories to handle system-bath interaction.
    
    
    Relaxation theories
    -------------------
    
    Specify the ``relaxation_theory`` attribute, e.g.
    
        self.relaxation_theory = "Lindblad_form"
    
    Implemented theories
    
    * Lindblad_form
    
    
    Simulation Tasks
    ----------------
    
    Implemented tasks
    
    * Population_dynamics
    
    """
    
    def setup(self):
        """Simulation definition
        
        
        In this method we define the simulation. We have to define molecular
        system, type of system-bath interaction and the theory to handel it, 
        and the data to be calculated. Below, the definition sections are 
        described in detail together with the definition options and required
        Quantarhei imports.
        
        
        
        """
        ###############################################################
        #
        #     General section
        #
        ###############################################################
        from quantarhei import TimeAxis
        
        # time interval of propagation
        self.timeaxis = TimeAxis(0.0, 1000, 1.0)
        
        ###############################################################
        #
        #     Molecular system definition section
        #
        ###############################################################
        from quantarhei import Molecule
        from quantarhei import energy_units
        
        # define as many molecules as you want
        # use energy units of your choice
        with energy_units("1/cm"):
            m1 = Molecule([0.0, 12000.0])
            m2 = Molecule([0.0, 12200.0])
            m3 = Molecule([0.0, 12000.0])
        
        # optionally you can specify positions (default units are Angstroms)
        m1.position = [0.0, 0.0, 0.0]
        m2.position = [10.0, 0.0, 0.0]
        m3.position = [0.0, 10.0, 0.0]
        
        # transition dipole moments for selected transitions
        m1.set_dipole(0, 1, [12.0, 0.0, 0.0])
        m2.set_dipole(0, 1, [8.0, 8.0, 0.0])
        m3.set_dipole(0, 1, [0.0, 8.0, 8.0])
        
        # couplings between molecules
        #
        #self.coupling_matrix = [[0.0, 100.0],
        #                        [100.0, 0.0]]
        #self.energy_units = "1/cm"
        #
        # or
        #
        self.coupling_matrix = "by_dipole_dipole"
        #
        
        # setting molecules is compulsory
        self.molecules = [m1, m2, m3]
        
        ##############################################################
        #
        #    Type of excitonic problem section
        #
        ##############################################################
        
        # here we define what kind of problem we calculate
        self.exciton_type = "electronic"
        #self.exciton_type = "vibronic"
        
        # exciton multiplicity specifies if we handle only single
        # exciton band (value 1) or also include two-exciton band (value 2)
        self.exciton_multiplicity = 2
        
        #############################################################
        #
        #    Relaxation theory
        #
        #############################################################
        
        #
        # Relaxation using Lindblad form
        #
        from quantarhei.qm import Operator
        
        o1 = Operator(dim=7, real=True)
        o1.data[1,2] = 1.0
        o1.data[2,3] = 1.0
        
        self.sys_operators = [o1]
        self.rates = (1.0/300,)
        self.relaxation_theory = "Lindblad_form"
        
        # Use it in site or excitonic basis
        self.use_sitebasis = True
        
        #
        # Redfield theory 
        #
        #self.relaxation_theory = "Standard_Redfield"
        #
        # Nothings else needs to be specified; it was defined on molecules
        #
        # Redfield is calculated strictly in excitonic basis
        # The  ``use_sitebasis`` attribute does not affect anything
        #
        
        
        #############################################################
        #
        #     Initial condition
        #
        #############################################################
        from quantarhei import ReducedDensityMatrix
        
        rho0 = ReducedDensityMatrix(dim=7)
        rho0.data[3,3] = 1.0
        
        self.rho0 = rho0

        #############################################################
        #
        #     Simulation tasks
        #
        #############################################################
        
        #task_1 = dict(task="density_matrix_dynamics", 
        #              object_file="rdm"+self._get_timestamp(filename=True)
        #              +".hdf5")

        task_1 = dict(task="density_matrix_dynamics") 
        
        self.tasks = [task_1]
