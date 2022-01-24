# -*- coding: utf-8 -*-
"""
    Basic data type of Quantarhei


"""
import numpy
import os

class MatrixData:
    """MatrixData type


    Parameters
    ----------

    dims : list or tuple
        Dimensions of the data in form of a tuple or list

    name : str
        Name of the object

    data : array like
        Object data

    """

    def __init__(self, dims=(0), name=None, data=None):

        if name is None:
            self.name = ""
        else:
            self.name = name

        if data is None:
            self.data = numpy.zeros(dims)
        else:
            if isinstance(data, list) or isinstance(data, tuple):
                self.data = numpy.array(data)
            else:
                self.data = data


    def set_name(self, name):
        """Sets the object name

        Parameters
        ----------

        name : str
            Name of of the object

        """
        self.name = name

    def get_name(self):
        """Returns the object name

        """
        return self.name

    def get_dim(self, n):
        """Returns dimension of the nth index of the matrix

        Parameters
        ----------

        n : int
            matrix index

        """
        return self.data.shape[n]

    def get_shape(self):
        """Returns the data shape

        """
        return self.data.shape

    def get_rank(self):
        """Returns the matrix rank, i.e. the number of its indices

        """
        return self.data.ndim

    def set_data(self, data):
        """Sets the data of the object

        Parameters
        ----------

        data : array like
            Data to be set in the object

        """
        self.data = data

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

        """
        filename, extension = os.path.splitext(name)

        if extension not in [".dat",".txt",".npy",".npz"]:
            raise Exception("Unknown data format")

        if (extension == ".dat") or (extension == ".txt"):
            self._exportDataToText(name)

        elif extension == ".npy":
            self._saveBinaryData(name)

        elif extension == ".npz":
            self._saveBinaryData_compressed(name)



    def load_data(self, name):
        """Loads the data in a format determined by the file name extension

        Parameters
        ----------

        name : str
            Name of the file to be loaded from

        """
        filename, extension = os.path.splitext(name)
        
        if extension not in [".dat",".txt",".npy",".npz"]:
            raise Exception("Unknown data format")        

        if (extension == ".dat") or (extension == ".txt"):
            self._importDataFromText(name)

        elif extension == ".npy":
            self._loadBinaryData(name)

        elif extension == ".npz":
            self._loadBinaryData_compressed(name)


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
        self.data = numpy.load(filename)

    def _loadBinaryData_compressed(self, filename):
        """Imports binary data from a file

        """          
        self.data = numpy.load(filename)["data"]
        
    def _exportDataToText(self, file):
        """Saves textual data to a file

        """
        numpy.savetxt(file, self.data)

    def _importDataFromText(self, filename):
        """Imports textual data to a file

        """        
        self.data = numpy.loadtxt(filename)
