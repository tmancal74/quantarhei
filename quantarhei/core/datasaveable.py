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
    
    
    def save_data(self, name, with_axis=None):
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
            self._exportDataToText(name, with_axis)

        elif extension == ".npy":
            self._saveBinaryData(name, with_axis)

        elif extension == ".npz":
            self._saveBinaryData_compressed(name, with_axis)
            
        elif extension == ".mat":
            self._saveMatlab(name, with_axis)



    def load_data(self, name, with_axis=None):
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
            self._importDataFromText(name, with_axis)

        elif extension == ".npy":
            self._loadBinaryData(name, with_axis)

        elif extension == ".npz":
            self._loadBinaryData_compressed(name, with_axis)

        elif extension == ".mat":
            self._loadMatlab(name, with_axis)

           
    def set_data_writable(self):
        """Implement this method to lift existing protection of data property
        
        """
        pass


    def set_data_protected(self):
        """Implement this method to put protections on data property
        
        """        
        pass


    def _data_with_axis(self, axis):
        """Constructs data array which contains also data from the axis
        
        """
        shpl = list(self.data.shape)
        
        if len(shpl) == 2:
            shpl[1] += 1
            shp = tuple(shpl)
            data = numpy.zeros(shp,dtype=self.data.dtype)
            data[:,1:] = self.data
            data[:,0] = axis.data     
        elif len(shpl) == 1:
            shpl.append(2)
            shp = tuple(shpl)
            data = numpy.zeros(shp,dtype=self.data.dtype)
            data[:,1] = self.data
            data[:,0] = axis.data
        else:
            raise Exception("Other shapes than (N,) and (N,M) not implemented")
        return data


    def _extract_data_with_axis(self, data, axis):
        """Extracts data part and the axis data from the `data` array 
        
        """
        if axis is None:
            return data
        else:
            if len(data.shape) == 2:
                
                print("Data shape:", data.shape)
                if data.shape[1] == 2:
                    print("Extracting from two columns")
                    axis.data = data[:,0]
                    return data[:,1]
                elif data.shape[1] > 2:
                    axis.data = data[:,0]
                    return data[:,1:]
                else:
                    raise Exception()
                    
            else:
                raise Exception("Other shapes than (N,) and (N,M)"+
                                " not implemented")



    def _saveBinaryData(self, file, with_axis=None):
        """Saves uncompressed binary data to an file

        """
        if with_axis is not None:
            data = self._data_with_axis(with_axis)
            numpy.save(file, data)
        else:
            numpy.save(file, self.data)


    def _saveBinaryData_compressed(self, file, with_axis=None):
        """Saves compressed binary data to an file

        """
        if with_axis is not None:
            data = self._data_with_axis(with_axis)
            numpy.save_compressed(file, data=data)
        else:
            numpy.savez_compressed(file, data=self.data)


    def _loadBinaryData(self, filename, with_axis=None):
        """Imports binary data from a file

        """
        
        self.set_data_writable()
        _data = numpy.load(filename)
        self.data = self._extract_data_with_axis(_data, with_axis)
        self.set_data_protected()


    def _loadBinaryData_compressed(self, filename, with_axis=None):
        """Imports binary data from a file

        """  
        self.set_data_writable()        
        _data = numpy.load(filename)["data"]
        self.data = self._extract_data_with_axis(_data, with_axis)
        self.set_data_protected()

    
    def _exportDataToText(self, file, with_axis=None):
        """Saves textual data to a file

        """
        if with_axis is not None:
            data = self._data_with_axis(with_axis)
            numpy.savetxt(file, data)
        else:
            numpy.savetxt(file, self.data)


    def _importDataFromText(self, filename, with_axis=None):
        """Imports textual data to a file

        """
        self.set_data_writable()
        try:        
            _data = numpy.loadtxt(filename)
        except ValueError:
            _data = numpy.loadtxt(filename, dtype=complex)
        
        self.data = self._extract_data_with_axis(_data, with_axis)
        self.set_data_protected()            


    def _saveMatlab(self, file, with_axis=None):
        """Saves data as a Matlab file
        
        """
        if with_axis is not None:
            data = self._data_with_axis(with_axis)
            io.savemat(file, {"data":data})
        else:
            io.savemat(file, {"data":self.data})

    
    def _loadMatlab(self, file, with_axis=None):
        """Loads a matrix called `data` from a matlab file
        
        """
        self.set_data_writable()
        _data = io.loadmat(file)["data"]
        self.data = self._extract_data_with_axis(_data, with_axis)
        self.set_data_protected()

