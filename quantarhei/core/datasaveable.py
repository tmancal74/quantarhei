# -*- coding: utf-8 -*-
"""
    Module datasaveable

    Defines the class DataSaveable. DataSaveable objects can save their data
    to and load them from a file in various formats. Data is exclusively
    represented as a `data` property of the class.
    

"""
import os
import numpy
import scipy.io as io

class DataSaveable:
    """This class defines saving and loading procedure for the data property
    
    """
    
    
    def save_data(self, name):
        """Saves the data into a format determined by the file name extension

        Parameters
        ----------

        name : str
            Name of the file to be saved into
            
        
        Notes
        -----
        
        This method knows the following file extensions
        
        .dat 
            text
            
        .txt
            text, same as .dat
            
        .npy
            binary numpy format, no compression
            
        .npz
            compressed numpy format
            
        .mat
            Matlab format

        """
        filename, extension = os.path.splitext(name)

        if extension not in [".dat",".txt",".npy",".npz",".mat"]:
            raise Exception("Unknown data format")

        if (extension == ".dat") or (extension == ".txt"):
            self._exportDataToText(name)

        elif extension == ".npy":
            self._saveBinaryData(name)

        elif extension == ".npz":
            self._saveBinaryData_compressed(name)
            
        elif extension == ".mat":
            self._saveMatlab(name)



    def load_data(self, name):
        """Loads the data in a format determined by the file name extension

        Parameters
        ----------

        name : str
            Name of the file to be loaded from

        """
        filename, extension = os.path.splitext(name)
        
        if extension not in [".dat",".txt",".npy",".npz", ".mat"]:
            raise Exception("Unknown data format")        

        if (extension == ".dat") or (extension == ".txt"):
            self._importDataFromText(name)

        elif extension == ".npy":
            self._loadBinaryData(name)

        elif extension == ".npz":
            self._loadBinaryData_compressed(name)

        elif extension == ".mat":
            self._loadMatlab(name)

           
    def set_data_writable(self):
        """Implement this method to lift existing protection of data property
        
        """
        pass


    def set_data_protected(self):
        """Implement this method to put protections on data property
        
        """        
        pass


    def _saveBinaryData(self, file):
        """Saves uncompressed binary data to an file

        """
        numpy.save(file, self.data)


    def _saveBinaryData_compressed(self, file):
        """Saves compressed binary data to an file

        """
        numpy.savez_compressed(file, data=self.data)


    def _loadBinaryData(self, filename):
        """Imports binary data from a file

        """
        
        self.set_data_writable()
        self.data = numpy.load(filename)
        self.set_data_protected()


    def _loadBinaryData_compressed(self, filename):
        """Imports binary data from a file

        """  
        self.set_data_writable()        
        self.data = numpy.load(filename)["data"]
        self.set_data_protected()

    
    def _exportDataToText(self, file):
        """Saves textual data to a file

        """
        numpy.savetxt(file, self.data)


    def _importDataFromText(self, filename):
        """Imports textual data to a file

        """
        self.set_data_writable()
        try:        
            self.data = numpy.loadtxt(filename)
        except ValueError:
            self.data = numpy.loadtxt(filename, dtype=complex)
        self.set_data_protected()            


    def _saveMatlab(self, file):
        """Saves data as a Matlab file
        
        """
        io.savemat(file, {"data":self.data})

    
    def _loadMatlab(self, file):
        """Loads a matrix called `data` from a matlab file
        
        """
        self.set_data_writable()
        self.data = io.loadmat(file)["data"]
        self.set_data_protected()

