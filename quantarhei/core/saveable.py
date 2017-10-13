# -*- coding: utf-8 -*-
import numpy
import h5py
import inspect
import fnmatch
import re
from numbers import Number

       
from .managers import energy_units

def _isattr(obj):
    return not inspect.ismethod(obj)
        
class Saveable:
    """Class implementing object persistence through saving to hdf5 files
    
    
    
    
    """
    
    def save(self, file):
        """Saves the object to a file
        
        
        The save is done by inspection of the object. Attributes are 
        automatically detected. Simple attributes (numbers, strigs, etc.) are
        saved, lists, dictionaries and tuples are saved recursively. All
        objects which are subclasses of Savable are saved too. Everything
        else is ignored.
        
    
        >>> class Test(Saveable):
        ...     a = 0.0
        ...     def __init__(self):
        ...       self.b = "ahoj"
        ...       self.c = ["A","B"]
        ...       self.d = dict(neco=0.5, jine="ahoj")
        ...     def set_T(self,TT):
        ...       self.obj = TT

        >>> T1 = Test()
        >>> T2 = Test()
        >>> T1.set_T(T2)
 
       
        """
         
        strings = {}
        numeric = {}
        boolean = {}
        numdata = {}
        objects = {}
        lists = {}
        tuples = {}
        dictionaries = {}
        
        # get all attributes to save
        
        rg = fnmatch.translate("__*__")
        prog = re.compile(rg)
        
        info = inspect.getmembers(self, _isattr)
        attr = []
        for fo in info:
            fo_name = fo[0]
            mch = prog.match(fo_name)
            if mch is None:
                attr.append(fo_name)

        for at_name in attr:
            at = self.__getattribute__(at_name)
            if isinstance(at,str):
                #print(at_name, at, "string")
                strings[at_name] = at
            elif isinstance(at, Number):
                #print(at_name, at, "number")
                numeric[at_name] = at
            elif isinstance(at, list):
                #print(at_name, at, "list")
                lists[at_name] = at
            elif isinstance(at, bool):
                #print(at_name, at, "boolean")
                boolean[at_name] = at
            elif isinstance(at, tuple):
                #print(at_name, at, "tuple")
                tuples[at_name] = at
            elif isinstance(at, dict):
                #print(at_name, at, "dictionary")
                dictionaries[at_name] = at
            elif isinstance(at, numpy.ndarray):
                #print(at_name, at, "numpy.array")
                numdata[at_name] = at
            elif isinstance(at, Saveable):
                #print(at_name, at, "Saveable object")
                objects[at_name] = at
            else:
                #print(at_name, at, "unknown") 
                pass
        
        self._do_save(file=file,
                      strings=strings, 
                      numeric=numeric, 
                      boolean=boolean,
                      numdata=numdata,
                      objects=objects,
                      lists=lists,
                      tuples=tuples,
                      dictionaries=dictionaries)



    def load(self, file):
        """This method should be implemented in classes that inherit from here
        
        """        
        strings = {}
        numeric = {}
        boolean = {}
        numdata = {}
        objects = {}
        lists = {}
        tuples = {}
        dictionaries = {}        
        
        self._do_load(file=file,
                      strings=strings, 
                      numeric=numeric, 
                      boolean=boolean,
                      numdata=numdata,
                      objects=objects,
                      lists=lists,
                      dictionaries=dictionaries,
                      tuples=tuples)

                              
            
    def _do_save(self, file="", strings={}, numeric={},
                 boolean={}, numdata={}, objects={},
                 lists={}, dictionaries={}, tuples={}):
        """Performs the save to hdf5 file
        
        
        """
        if file == "":
            raise Exception("No file specified")
            
        file_openned = False
    
        # check if we should open a file
        if isinstance(file, str):
            file_openned = True
            f = h5py.File(file,"w")
            root = self._create_root_group(f,str(self.__class__))
                
        else:
            root = self._create_root_group(file,str(self.__class__))

        # do the save of all dictionaries
        with energy_units("int"):            
            self._save_classid(root, self.__class__)
            self._save_strings(root, strings)
            self._save_numeric(root, numeric)
            self._save_boolean(root, boolean)
            self._save_numdata(root, numdata)
            self._save_objects(root, objects)
            self._save_lists(root, lists)
            self._save_tuples(root, tuples)
            self._save_dictionaries(root, dictionaries)
    
    
        # if file openned here, close it
        if file_openned:
            f.close()
                
        
    def _do_load(self, file="", strings={}, numeric={},
                 boolean={}, numdata={}, objects={},
                 lists={}, dictionaries={}, tuples={}):
        """Loads data from hdf5 file
        
        
        """        
        if file == "":
            raise Exception("No file specified")
            
        file_openned = False
    
        # check if we should open a file
        if isinstance(file, str):
            file_openned = True
            f = h5py.File(file,"r")
            root = f[str(self.__class__)]
                
        else:
            root = file[str(self.__class__)]
            
        cid = self._load_classid(root)
        if cid != str(self.__class__):
            raise Exception("File contains an object of unexpected type")
        
        self._load_strings(root, strings)
        self._load_numeric(root, numeric)
        self._load_boolean(root, boolean)
        self._load_numdata(root, numdata)
        self._load_objects(root, objects)
        self._load_lists(root, lists)
        self._load_tuples(root, tuples)
        self._load_dictionaries(root, dictionaries)     
        
        for key in strings.keys():
            setattr(self,key,strings[key])
        for key in numeric.keys():
            setattr(self,key,numeric[key])
        for key in numdata.keys():
            setattr(self,key,numdata[key])
        for key in boolean.keys():
            setattr(self,key,boolean[key])            
            
        
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
        loc.attrs.create("classid",numpy.string_(str(classid)))
        print(classid)
    
    
    def _save_strings(self, loc, dictionary):
        """Saves a dictionary of strings under the group "strings"
        
        """
        if len(dictionary) == 0: return
        strs = self._create_root_group(loc,"strings")
        for key in dictionary.keys(): 
            strs.attrs.create(key, numpy.string_(dictionary[key]))
        print(dictionary.keys())
    
    def _save_numeric(self, loc, dictionary):
        """Saves a dictionary of numeric values under the group "numeric"
        
        """
        if len(dictionary) == 0: return
        strs = self._create_root_group(loc,"numeric")
        for key in dictionary.keys(): 
            strs.attrs.create(key, dictionary[key])
        print(dictionary.keys())
    
    def _save_boolean(self, loc, dictionary):
        """Saves a dictionary of boolean values under the group "boolean"
        
        """
        if len(dictionary) == 0: return
        strs = self._create_root_group(loc,"boolean")
        for key in dictionary.keys(): 
            strs.attrs.create(key, dictionary[key])
        print(dictionary.keys())
            
            
    def _save_numdata(self, loc, dictionary):
        """Saves a dictionary of numpy arrays under the group "numdata"
        
        """
        if len(dictionary) == 0: return
        strs = self._create_root_group(loc,"numdata")
        for key in dictionary.keys(): 
            strs.create_dataset(key, data=dictionary[key])
        print(dictionary.keys())
    
    def _save_objects(self, loc, dictionary):
        """Saves a dictionary of objects under the group "objects"
        
        """
        if len(dictionary) == 0: return
        strs = self._create_root_group(loc,"objects")
        for key in dictionary.keys(): 
            obj = dictionary[key]
            obj.save(strs)
        print(dictionary.keys())
    
    def _save_lists(self, loc, dictionary):
        """Saves a dictionary of lists under the group "lists"
        
        """
        if len(dictionary) == 0: return
        strs = self._create_root_group(loc,"lists")
        for key in dictionary.keys(): 
            clist = dictionary[key]
            self._save_a_list(strs, clist)
        print(dictionary.keys())

    def _save_tuples(self, loc, dictionary):
        """Saves a dictionary of lists under the group "lists"
        
        """
        if len(dictionary) == 0: return
        strs = self._create_root_group(loc,"tuples")
        for key in dictionary.keys(): 
            clist = dictionary[key]
            self._save_a_tuple(strs, clist)
        print(dictionary.keys())
        
    def _save_dictionaries(self, loc, dictionary):
        """Saves a dictionary of dictionaries under the group "dicts"
        
        """
        if len(dictionary) == 0: return
        strs = self._create_root_group(loc,"dicts")
        for key in dictionary.keys(): 
            cdict = dictionary[key]
            self._save_a_dictionary(strs, cdict)           
        print(dictionary.keys())



            
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
        return loc.attrs["classid"].decode("utf-8")
    
    
    def _load_strings(self, loc, dictionary):
        """Load a dictionary of strings from the group "strings"
        
        """
        try:
            strs = loc["strings"]
        except:
            return
        
        for key in strs.attrs.keys(): 
            dictionary[key] = strs.attrs[key].decode("utf-8")
    
    def _load_numeric(self, loc, dictionary):
        """Saves a dictionary of numeric values under the group "numeric"
        
        """
        try:
            strs = loc["numeric"]
        except:
            return
        
        for key in strs.attrs.keys(): 
            dictionary[key] = strs.attrs[key]
    
    def _load_boolean(self, loc, dictionary):
        """Saves a dictionary of boolean values under the group "boolean"
        
        """
        try:
            strs = loc["boolean"]
        except:
            return
        
        for key in strs.attrs.keys(): 
            dictionary[key] = strs.attrs[key]
            
            
    def _load_numdata(self, loc, dictionary):
        """Saves a dictionary of numpy arrays under the group "numdata"
        
        """
        try:
            strs = loc["numdata"]
        except:
            return
        
        for key in strs.keys(): 
            dictionary[key] = numpy.array(strs[key])
            
    
    def _load_objects(self, loc, dictionary):
        """Saves a dictionary of objects under the group "objects"
        
        """
        try:
            strs = loc["objects"]
        except:
            return
        
        #for key in dictionary.keys(): 
        #    obj = dictionary[key]
        #    obj.save(strs)
    
    def _load_lists(self, loc, dictionary):
        """Saves a dictionary of lists under the group "lists"
        
        """
        try:
            strs = loc["lists"]
        except:
            return
        
        for key in dictionary.keys(): 
            clist = dictionary[key]
            #self._save_a_list(strs, clist)
 
    def _load_dictionaries(self, loc, dictionary):
        """Saves a dictionary of dictionaries under the group "dicts"
        
        """
        try:
            strs = loc["dicts"]
        except:
            return
        
        for key in dictionary.keys(): 
            cdict = dictionary[key]
            #self._save_a_dictionary(strs, cdict)           

    def _load_tuples(self, loc, dictionary):
        """Saves a dictionary of dictionaries under the group "dicts"
        
        """
        try:
            strs = loc["tuples"]
        except:
            return
        
        for key in dictionary.keys(): 
            cdict = dictionary[key]
            #self._save_a_dictionary(strs, cdict)

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
        
           