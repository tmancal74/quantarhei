# -*- coding: utf-8 -*-

from ...builders.aggregates import Aggregate
from .simulation import Simulation

class ExcitonDynamics(Simulation):
    """Excitonic Dynamics Simulation
    
    
    """        
        
    def _build(self):
        """Here we prepare objects for the simulation
        
        """
        self._printlog("Building aggregate", loglevel=0)
                
        try:
            self.aggregate = Aggregate(molecules=self.molecules)
        except:
            raise Exception("Aggregate construction failed")
            
        #try:
        #    if self.coupling_matrix
            
        try:
            self.aggregate.build(mult=self.exciton_multiplicity)
        except AttributeError:
            raise Exception("Aggregate build failed: exciton_multiplicity not defined")
            
        self._incr_indent_level()
        self._printlog("Number of monomers          :", self.aggregate.nmono, loglevel=0)
        self._printlog("Number of electronic states :", self.aggregate.Nel, loglevel=0)
        self._printlog("Total number of states      :", self.aggregate.Ntot, loglevel=0)
        self._decr_indent_level()
        
        #
        #  Relaxation theory section
        #
        #  Lindblad is done all here. Redfield is calculated in implementation method
        #
        if self.relaxation_theory == "Lindblad_form":
            
            #from ..qm import LindbladForm
            from ...qm import SystemBathInteraction
            
            ham = self.aggregate.get_Hamiltonian()
            try:
                sysops = self.sys_operators
            except AttributeError:
                raise Exception("``sys_operators`` have to be speficied for Lindblad form")
            try:
                rates = self.rates
            except AttributeError:
                raise Exception("``rates`` have to be speficied for Lindblad form")
                
            self._printlog("Relaxation theory in Lindblad form", loglevel=0)
            sbi = SystemBathInteraction(sys_operators=sysops, rates=rates)
            self.aggregate.set_SystemBathInteraction(sbi)
            # FIXME: use self.use_sitebasis to allow definition of Lindblad in excitonic basis
            self._relaxation_tensor, ham = self.aggregate.get_RelaxationTensor(self.timeaxis,
                                                                               relaxation_theory="Lindblad_form")
            self._ham = ham
            self._printlog("...tensor constructed", loglevel=0)

            
        
    def _implementation(self):
        """Here the main simulation tasks are performed
        
        """
        self._printlog("Simulating excitonic dynamics", loglevel=0)
        
        if self.relaxation_theory == "Standard_Redfield":
            # set up Redfield tensor
            pass
        
        elif self.relaxation_theory == "Standard_Foerster":
            # set up Foerster tensor
            pass
        
        elif self.relaxation_theory == "Lindblad_form":
            # This is already set up
            self._relaxation_tensor.convert_2_tensor() 
        
        tasks = self.tasks
        
        for task in tasks:
            
            if task["task"] == "density_matrix_dynamics":
                
                from ...qm import ReducedDensityMatrixPropagator
                
                self._printlog("Calculating density matrix dynamics")
                rdmprop = ReducedDensityMatrixPropagator(timeaxis=self.timeaxis,
                                                         Ham=self._ham, 
                                                         RTensor=self._relaxation_tensor)
                
                rhot = rdmprop.propagate(self.rho0)
                self._rhot = rhot
                
                try:
                    objfile = task["object_file"]
                    self._printlog("Saving RDM dynamics into file:", objfile, loglevel=0)
                    self._rhot.save(objfile)
                except:
                    pass
                
                self._printlog("...done")
                
            

