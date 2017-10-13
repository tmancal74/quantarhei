# -*- coding: utf-8 -*-
import numpy
import h5py

from .managers import energy_units

class Saveable:
    """Class implementing object persistence through saving to hdf5 files
    
    
    
    
    """
    
    def save(self, file):
        """Saves the object to a file
        
        
        This method should be implemented in classes that inherit from Savable.
        Prepare the following dictionaries
        
        strings
        numeric
        boolean
        numdata
        objects
        
        Always use the name of the attribute you are saving as a key. If your
        class has an integer attribute N, like this one
        
            class Test():
            
                def __init__(self, id):
                    self.id = id
                
        you should have something this in your save method
        
            num["id"] = self.id
        
        This dictionary will be submitted to the self._do_save method together 
        with all other dictionaries.
        
        ```data``` should be numpy arrays and objects must inherit from the 
        Savable class. If all these requirements are met, calling self._do_save
        method with all the above arguments, e.g.
        
            self._do_save(filename=fname, strings=strs, numeric=num, ... )
        
        
        """
         
        strings = {}
        numeric = {}
        boolean = {}
        numdata = {}
        objects = {}
        lists = {}
        dictionaries = {}
        
        
        self._do_save(file=file,
                      strings=strings, 
                      numeric=numeric, 
                      boolean=boolean,
                      numdata=numdata,
                      objects=objects,
                      lists=lists,
                      dictionaries=dictionaries)
                              
            
    def _do_save(self, file="", strings={}, numeric={},
                 boolean={}, numdata={}, objects={},
                 lists={}, dictionaries={}):
        """Performs the save to hdf5 file
        
        
        """
        if file == "":
            raise Exception("No file specified")
            
        file_openned = False
    
        # check if we should open a file
        if isinstance(file, str):
            file_openned = True
            f = h5py.File(file,"w")
            root = self._create_root_group(f,self.__class__)
                
        else:
            root = self._create_root_group(file,self.__class__)

        # do the save of all dictionaries
        with energy_units("int"):            
            self._save_classid(root, self.__class__)
            self._save_strings(root, strings)
            self._save_numeric(root, numeric)
            self._save_boolean(root, boolean)
            self._save_numdata(root, numdata)
            self._save_objects(root, objects)
            self._save_lists(root, lists)
            self._save_dictionaries(root, dictionaries)
    
    
        # if file openned here, close it
        if file_openned:
            f.close()
        
    def load(self, file):
        """This method should be implemented in classes that inherit from here
        
        """        
        strings = {}
        numeric = {}
        boolean = {}
        numdata = {}
        objects = {}
        lists = {}
        dictionaries = {}        
        
        self._do_load(file=file,
                      strings=strings, 
                      numeric=numeric, 
                      boolean=boolean,
                      numdata=numdata,
                      objects=objects,
                      lists=lists,
                      dictionaries=dictionaries)

        
        
        
        
    def _do_load(self, file="", strings={}, numeric={},
                 boolean={}, numdata={}, objects={},
                 lists={}, dictionaries={}):
        """Loads data from hdf5 file
        
        
        """        
        if file == "":
            raise Exception("No file specified")
            
        file_openned = False
    
        # check if we should open a file
        if isinstance(file, str):
            file_openned = True
            f = h5py.File(file,"r")
            root = f[self.__class__]
                
        else:
            root = file[self.__class__]
            
        cid = self._load_classid(root)
        if cid != self.__class__:
            raise Exception("File contains an object of unexpected type")
            
        self._load_strings(root, strings)
        self._load_numeric(root, numeric)
        self._load_boolean(root, boolean)
        self._load_numdata(root, numdata)
        self._load_objects(root, objects)
        self._load_lists(root, lists)
        self._load_dictionaries(root, dictionaries)            
        
        # if file openned here, close it
        if file_openned:
            f.close()
        
    #
    # helper methods
    #
    
    def _create_root_group(self, start, name):
        """Creates the root "directory" of the tree to save
        
        """
        return start.create_group(name) 
    
    #
    # saving methods
    #
    
    def _save_classid(self, loc, classid):
        """Saves a string as "classid" attribute
        
        """
        loc.attrs.create("classid",classid)
    
    
    def _save_strings(self, loc, dictionary):
        """Saves a dictionary of strings under the group "strings"
        
        """
        strs = self._create_root_group(loc,"strings")
        for key in dictionary.keys(): 
            strs.attrs.create(key, numpy.string_(dictionary[key]))
    
    def _save_numeric(self, loc, dictionary):
        """Saves a dictionary of numeric values under the group "numeric"
        
        """
        strs = self._create_root_group(loc,"numeric")
        for key in dictionary.keys(): 
            strs.attrs.create(key, dictionary[key])
    
    def _save_boolean(self, loc, dictionary):
        """Saves a dictionary of boolean values under the group "boolean"
        
        """
        strs = self._create_root_group(loc,"boolean")
        for key in dictionary.keys(): 
            strs.attrs.create(key, dictionary[key])
            
            
    def _save_numdata(self, loc, dictionary):
        """Saves a dictionary of numpy arrays under the group "numdata"
        
        """
        strs = self._create_root_group(loc,"numdata")
        for key in dictionary.keys(): 
            strs.create_dataset(key, dictionary[key])
    
    def _save_objects(self, loc, dictionary):
        """Saves a dictionary of objects under the group "objects"
        
        """
        strs = self._create_root_group(loc,"objects")
        for key in dictionary.keys(): 
            obj = dictionary[key]
            obj.save(strs)
    
    def _save_lists(self, loc, lists):
        """Saves a dictionary of lists under the group "lists"
        
        """
        strs = self._create_root_group(loc,"lists")
        for key in lists.keys(): 
            clist = lists[key]
            self._save_a_list(strs, clist)
 
    def _save_dictionaries(self, loc, dictionaries):
        """Saves a dictionary of dictionaries under the group "dicts"
        
        """
        strs = self._create_root_group(loc,"dicts")
        for key in dictionaries.keys(): 
            cdict = dictionaries[key]
            self._save_a_dictionary(strs, cdict)           
            
    def _save_a_list(self, loc, clist):
        """Saves a list of values
        
        All simple values (strings, numeric, boolean), numpy arrays, Savable
        objects, and also lists and dictionaries should be acceptable
        
        """
        pass
    
    def _save_a_dictionary(self, loc, cdict):
        """Saves a dictionary of values
        
        All simple values (strings, numeric, boolean), numpy arrays, Savable
        objects, and also lists and dictionaries should be acceptable
        
        """
        pass
    
    
    
    #
    # loading methods
    #   
    
    def _load_classid(self, loc):
        """Loads a string as "classid" attribute
        
        """
        return loc.attrs["classid"]
    
    
    def _load_strings(self, loc, dictionary):
        """Load a dictionary of strings from the group "strings"
        
        """
        strs = loc["strings"]
        for key in strs.attrs.keys(): 
            dictionary[key] = strs.attrs[key].decode("utf-8")
    
    def _load_numeric(self, loc, dictionary):
        """Saves a dictionary of numeric values under the group "numeric"
        
        """
        strs = self._create_root_group(loc,"numeric")
        for key in dictionary.keys(): 
            strs.attrs.create(key, dictionary[key])
    
    def _load_boolean(self, loc, dictionary):
        """Saves a dictionary of boolean values under the group "boolean"
        
        """
        strs = self._create_root_group(loc,"boolean")
        for key in dictionary.keys(): 
            strs.attrs.create(key, dictionary[key])
            
            
    def _load_numdata(self, loc, dictionary):
        """Saves a dictionary of numpy arrays under the group "numdata"
        
        """
        strs = self._create_root_group(loc,"numdata")
        for key in dictionary.keys(): 
            strs.create_dataset(key, dictionary[key])
    
    def _load_objects(self, loc, dictionary):
        """Saves a dictionary of objects under the group "objects"
        
        """
        strs = self._create_root_group(loc,"objects")
        for key in dictionary.keys(): 
            obj = dictionary[key]
            obj.save(strs)
    
    def _load_lists(self, loc, lists):
        """Saves a dictionary of lists under the group "lists"
        
        """
        strs = self._create_root_group(loc,"lists")
        for key in lists.keys(): 
            clist = lists[key]
            self._save_a_list(strs, clist)
 
    def _load_dictionaries(self, loc, dictionaries):
        """Saves a dictionary of dictionaries under the group "dicts"
        
        """
        strs = self._create_root_group(loc,"dicts")
        for key in dictionaries.keys(): 
            cdict = dictionaries[key]
            self._save_a_dictionary(strs, cdict)           
            
    def _load_a_list(self, loc, clist):
        """Saves a list of values
        
        All simple values (strings, numeric, boolean), numpy arrays, Savable
        objects, and also lists and dictionaries should be acceptable
        
        """
        pass
    
    def _load_a_dictionary(self, loc, cdict):
        """Saves a dictionary of values
        
        All simple values (strings, numeric, boolean), numpy arrays, Savable
        objects, and also lists and dictionaries should be acceptable
        
        """
        pass    
     
        
    
class TestObject(Saveable):
    """
    
    

    
    """
    def __init__(self):
        
        self.name = "Ahoj"
        self.class_name = "TestObject"
        self.N = 23
        
        
    def save(self, filename):
        
        # strings
        strs = dict(name=self.name,
                    class_name=self.class_name) 
        # integers
        ints = dict(N=self.N)
        
        with h5py.File(filename,"w") as file:

            #
            # Saving all types
            #
            self.save_string_attributes(file, strs)
            self.save_numeric_attributes(file, ints)
        
        
        
    def load(self, file):
        
        string_attr_names = ["name", "class_name"]
        string_attrs = self.load_string_attributes(file, string_attr_names)
        self.name = string_attrs["name"]
        self.class_name = string_attrs["class_name"]
        
        numeric_attr_names = ["N"]
        numeric_attrs = self.load_numeric_attributes(file, numeric_attr_names)
        self.N = numeric_attrs["N"]
        
                