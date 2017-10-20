# -*- coding: utf-8 -*-
import numpy
import h5py
import inspect
import fnmatch
import re
from numbers import Number

       
from .managers import energy_units

def _isattr(obj):
    ismethod = inspect.ismethod(obj)
    return not ismethod 

     
simple_types = ("numeric", "strings", "boolean")        

class Saveable:
    """Class implementing object persistence through saving to hdf5 files
    
    
    
    
    """
    
    def _before_save(self):
        pass
    
    def _after_save(self):
        pass
    
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
        
        #
        # Before save start-up
        #
        self._before_save()
         
        
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
        
        # members that are attributes
        info = inspect.getmembers(self, _isattr)
        
        # members that are properties
        propts = inspect.getmembers(self.__class__, 
                                    lambda o: isinstance(o, property))
        
        
        prop_names = []
        for p in propts:
            prop_names.append(p[0])
            
        attr = []
        for fo in info:
            fo_name = fo[0]
            mch = prog.match(fo_name)
            if (mch is None) and (fo_name not in prop_names):
                # attributes starting with _S__ are protected and will not be saved
                if not fo_name.startswith("_S__"):
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

        # 
        # After save clean-up
        #
        self._after_save()


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
    
        fullclassid = self._get_full_class_name(self)
    
        # check if we should open a file
        if isinstance(file, str):
            file_openned = True
            f = h5py.File(file,"w")
            root = self._create_root_group(f,fullclassid)
                
        else:
            root = self._create_root_group(file,fullclassid)

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
          
    def _get_full_class_name(self, obj):
        return obj.__class__.__module__+"."+obj.__class__.__name__
        
    def _do_load(self, file="", strings={}, numeric={},
                 boolean={}, numdata={}, objects={},
                 lists={}, dictionaries={}, tuples={}):
        """Loads data from hdf5 file
        
        
        """        
        if file == "":
            raise Exception("No file specified")
            
        file_openned = False
        thisclass = self._get_full_class_name(self)
        
        # check if we should open a file
        if isinstance(file, str):
            file_openned = True
            f = h5py.File(file,"r")
            root = f[thisclass]
                
        else:
            root = file[thisclass]
            
        cid = self._load_classid(root)
        
        
        if cid != thisclass:
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
        for key in lists.keys():
            setattr(self,key,lists[key]) 
        for key in objects.keys():
            setattr(self,key,objects[key])
        for key in dictionaries.keys():
            setattr(self,key,dictionaries[key])            
        
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
        fullclassid = classid.__module__+"."+classid.__name__
        loc.attrs.create("classid",numpy.string_(fullclassid))
        #print(classid)
    
    
    def _save_strings(self, loc, dictionary):
        """Saves a dictionary of strings under the group "strings"
        
        """
        if len(dictionary) == 0: return
        strs = self._create_root_group(loc,"strings")
        for key in dictionary.keys(): 
            strs.attrs.create(key, numpy.string_(dictionary[key]))
        #print(dictionary.keys())
    
    def _save_numeric(self, loc, dictionary):
        """Saves a dictionary of numeric values under the group "numeric"
        
        """
        if len(dictionary) == 0: return
        strs = self._create_root_group(loc,"numeric")
        for key in dictionary.keys(): 
            strs.attrs.create(key, dictionary[key])
        #print(dictionary.keys())
    
    def _save_boolean(self, loc, dictionary):
        """Saves a dictionary of boolean values under the group "boolean"
        
        """
        if len(dictionary) == 0: return
        strs = self._create_root_group(loc,"boolean")
        for key in dictionary.keys(): 
            strs.attrs.create(key, dictionary[key])
        #print(dictionary.keys())
            
            
    def _save_numdata(self, loc, dictionary):
        """Saves a dictionary of numpy arrays under the group "numdata"
        
        """
        if len(dictionary) == 0: return
        strs = self._create_root_group(loc,"numdata")
        for key in dictionary.keys(): 
            strs.create_dataset(key, data=dictionary[key])
        #print(dictionary.keys())
    
    def _save_objects(self, loc, dictionary):
        """Saves a dictionary of objects under the group "objects"
        
        """
        if len(dictionary) == 0: return
        strs = self._create_root_group(loc,"objects")
        for key in dictionary.keys(): 
            obj = dictionary[key]
            nroot = self._create_root_group(strs,key)
            obj.save(nroot)
        #print(dictionary.keys())
    
    def _save_lists(self, loc, dictionary):
        """Saves a dictionary of lists under the group "lists"
        
        """
        if len(dictionary) == 0: return
        strs = self._create_root_group(loc,"lists")
        for key in dictionary.keys(): 
            clist = dictionary[key]
            self._save_a_listuple(strs, key, "list", clist)
        #print(dictionary.keys())

    def _save_tuples(self, loc, dictionary):
        """Saves a dictionary of lists under the group "lists"
        
        """
        if len(dictionary) == 0: return
        strs = self._create_root_group(loc,"tuples")
        for key in dictionary.keys(): 
            clist = dictionary[key]
            self._save_a_listuple(strs, key, "tuple", clist)
        #print(dictionary.keys())
        
    def _save_dictionaries(self, loc, dictionary):
        """Saves a dictionary of dictionaries under the group "dicts"
        
        """
        if len(dictionary) == 0: return
        strs = self._create_root_group(loc,"dicts")
        for key in dictionary.keys(): 
            cdict = dictionary[key]
            self._save_a_dictionary(strs, key, cdict)           
        #print(dictionary.keys())


    def _get_type(self, obj):
        """Returns a simple type or None
        
        """
        htyp = None
        if isinstance(obj, Number):
            htyp = "numeric"
        elif isinstance(obj, str):
            htyp = "strings"
        elif isinstance(obj, bool):
            htyp = "boolean"
            
        return htyp
        

    def _is_homogeneous_simple(self, liple):
        """Returns True and type of the content of the list if the list
        contains objects of the same type and the type is one of the simple
        types: numeric, boolean or strings
        
        """
        
        obj = liple[0]
        htyp = self._get_type(obj)
        homo = True
        
        if htyp is None:
            homo = False
            return (homo, htyp)
        
        for obj in liple:
            if htyp != self._get_type(obj): 
                htyp = None
                homo = False
                break
            
        return (homo, htyp)


    def _save_a_listuple(self, loc, key, typ, ctuple):
        """ Saves a list or tuple of items
        
        All simple values (strings, numeric, boolean), numpy arrays, Savable
        objects, and also lists, tuples and dictionaries should be acceptable
        
        """
        
        root = self._create_root_group(loc, key)
        root.attrs.create("type", numpy.string_(typ))
        
        (homo, htyp) = self._is_homogeneous_simple(ctuple)
        
        root.attrs.create("homogeneous_simple", homo)
        
        if homo:
            if (htyp == "numeric") or (htyp == "boolean"):
                root.create_dataset(htyp, data=numpy.array(ctuple))
                
            elif htyp == "strings":
                # get the maximum length of the string
                slen = 0
                for item in ctuple:
                    nlen = len(item)
                    if nlen > slen:
                        slen = nlen

                dtp = "a"+str(slen)
                dset = root.create_dataset(htyp, (len(ctuple),), dtype=dtp)
                dset[:] = numpy.array(numpy.string_(ctuple))
            else:
                raise Exception("Attempt to save an empty list")
        
        else:
            root.attrs.create("length",len(ctuple))
            k = 0
            for item in ctuple:
                tupkey = str(k)
                if self._in_simple_types(item):
                    self._save_a_simple_type(root, tupkey, item)
                elif isinstance(item, Saveable):
                    objroot = self._create_root_group(root,tupkey)
                    objroot.attrs.create("type", numpy.string_("Saveable"))
                    item.save(objroot)
                    #print("Object saved at:",objroot)
                    #for key in objroot.keys():
                    #    print("                ",key)
                    
                elif isinstance(item, list):
                    self._save_a_listuple(root, tupkey, "list", item)
                elif isinstance(item, tuple):
                    self._save_a_listuple(root, tupkey, "tuple", item)
                elif isinstance(item, dict):
                    self._save_a_dictionary(root, tupkey, item)
                else:
                    # possibly warn that something was not saved
                    pass
                k += 1
    
    def _in_simple_types(self, item):
        if isinstance(item, Number):
            return True
        if isinstance(item, str):
            return True
        return False
    
    def _save_a_simple_type(self, loc, key, item):
        
        if isinstance(item,Number):
            grp = self._create_root_group(loc,key)
            grp.attrs.create("item",item)
            grp.attrs.create("type",numpy.string_("numeric"))
        elif isinstance(item,bool):
            grp = self._create_root_group(loc,key)
            grp.attrs.create("item",item)
            grp.attrs.create("type",numpy.string_("boolean"))
        elif isinstance(item,str):
            grp = self._create_root_group(loc,key)
            grp.attrs.create("item",numpy.string_(item))
            grp.attrs.create("type",numpy.string_("strings"))


    def _load_a_simple_type(self, loc):
        
        typ = loc.attrs["type"]
        typ = typ.decode("utf-8")
        if typ == "numeric":
            return loc.attrs["item"]
        if typ == "strings":
            ss = loc.attrs["item"]
            return ss.decode("utf-8")
        if typ == "boolean":
            return loc.attrs["item"]
        
        
    
    def _save_a_dictionary(self, loc, key, cdict):
        """Saves a dictionary of values
        
        All simple values (strings, numeric, boolean), numpy arrays, Savable
        objects, and also lists and dictionaries should be acceptable
        
        """
        root = self._create_root_group(loc, key)
        root.attrs.create("type", numpy.string_("dict"))
        
        for dkey in cdict.keys():
            item  = cdict[dkey]
            if self._in_simple_types(item):
                self._save_a_simple_type(root, dkey, item)
            elif isinstance(item, Saveable):
                objroot = self._create_root_group(root,dkey)
                objroot.attrs.create("type", numpy.string_("Saveable"))
                item.save(objroot)

            elif isinstance(item, list):
                self._save_a_listuple(root, dkey, "list", item)
            elif isinstance(item, tuple):
                self._save_a_listuple(root, dkey, "tuple", item)
            elif isinstance(item, dict):
                self._save_a_dictionary(root, dkey, item)
            else:
                # possibly warn that something was not saved
                pass
    
    
    
    #
    # loading methods
    #   
    
    def _load_classid(self, loc):
        """Loads a string as "classid" attribute
        
        """
        fullclassid = loc.attrs["classid"].decode("utf-8")
        return fullclassid
    
    
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
          
    def _get_class(self, kls):
        
        parts = kls.split('.')

        module = ".".join(parts[:-1])
        m = __import__( module )
        parts = parts[1:]

        for comp in parts:
            m = getattr(m, comp)            
        return m
    
    def _load_objects(self, loc, dictionary):
        """Saves a dictionary of objects under the group "objects"
        
        """
        try:
            strs = loc["objects"]
        except:
            return
        
        for key in strs.keys():
            
            objloc = strs[key]
            k = 0
            for kex in objloc.keys():
                objcls = self._get_class(kex)
                k += 1
            obj = objcls()
            obj.load(objloc)
            dictionary[key] = obj
    
    def _load_lists(self, loc, dictionary):
        """Saves a dictionary of lists under the group "lists"
        
        """
        try:
            strs = loc["lists"]
        except:
            return
        
        for key in  strs.keys():
            loc = strs[key]
            clist = self._load_a_listuple(loc)
            dictionary[key] = clist
            
 
    def _load_dictionaries(self, loc, dictionary):
        """Saves a dictionary of dictionaries under the group "dicts"
        
        """
        try:
            strs = loc["dicts"]
        except:
            return
        
        for key in strs.keys(): 
            loc = strs[key]
            cdict = self._load_a_dictionary(loc)
            dictionary[key] = cdict          

    def _load_tuples(self, loc, dictionary):
        """Saves a dictionary of dictionaries under the group "dicts"
        
        """
        try:
            strs = loc["tuples"]
        except:
            return
        
        for key in  strs.keys():
            loc = strs[key]
            clist = self._load_a_listuple(loc)
            dictionary[key] = clist
            

    def _load_a_listuple(self, loc):
        """Loads a list or tuple of values
        
        All simple values (strings, numeric, boolean), numpy arrays, Savable
        objects, and also lists and dictionaries should be acceptable
        
        """
        
        typ = loc.attrs["type"].decode("utf-8")
        clist = []
        homo = loc.attrs["homogeneous_simple"]
        
        if homo:

            if ("numeric" in loc.keys()):
                data = loc["numeric"]
            elif ("boolean" in loc.keys()):
                data = loc["boolean"]
            elif ("strings" in loc.keys()):
                data = []
                for dt in loc["strings"]:
                    data.append(dt.decode("utf-8"))
                    
            else:
                for key in loc.keys():
                    print(key)
                print("length: ", len(loc.keys()))
                raise Exception()
            
            
            for dt in data:
                clist.append(dt)
                
        else:
            
            N = loc.attrs["length"]
            for k in range(N):
                tupkey = str(k)
                itemloc = loc[tupkey]
                ltyp = itemloc.attrs["type"].decode("utf-8")
                if ltyp in simple_types:
                    item = self._load_a_simple_type(itemloc)
                    clist.append(item)
                elif ltyp ==  "Saveable":

                    for key in itemloc.keys():
                        classname = key

                    cls = self._get_class(classname)
                    obj = cls()
                    obj.load(itemloc)
                    clist.append(obj)

                elif ltyp == "list":
                    item = self._load_a_listuple(itemloc)
                    clist.append(item)
                elif ltyp == "tuple":
                    item = self._load_a_listuple(itemloc)
                    clist.append(item)
                elif ltyp == "dict":
                    item = self._load_a_dictionary(itemloc)
                    clist.append(item)
                else:
                    # possibly warn that something was not saved
                    pass
                
                
                k += 1 
                
        # if the type is tuple, we convert list to tuple
        if typ == "tuple":
            #print("convering to tuple")
            clist = tuple(clist)
            
        return clist
            
    
    def _load_a_dictionary(self, loc):
        """Loads a dictionary of values
        
        All simple values (strings, numeric, boolean), numpy arrays, Savable
        objects, and also lists and dictionaries should be acceptable
        
        """
        #typ = loc.attrs["type"].decode("utf-8")
        cdict = {}

        for dkey in loc.keys():
            itemloc = loc[dkey]
            ltyp = itemloc.attrs["type"].decode("utf-8")

            if ltyp in simple_types:
                item = self._load_a_simple_type(itemloc)
                cdict[dkey] = item
            elif ltyp ==  "Saveable":

                for key in itemloc.keys():
                    classname = key

                cls = self._get_class(classname)
                obj = cls()
                obj.load(itemloc)
                cdict[dkey] = obj

            elif ltyp == "list":
                item = self._load_a_listuple(itemloc)
                cdict[dkey] = item
            elif ltyp == "tuple":
                item = self._load_a_listuple(itemloc)
                cdict[dkey] = item
            elif ltyp == "dict":
                item = self._load_a_dictionary(itemloc)
                cdict[dkey] = item
            else:
                # possibly warn that something was not saved
                pass
    
        return cdict
    
    
class TSaveable(Saveable):
    
    def __init__(self):
        
        self.a = 1.0
        self.text = None
        self._txt = None
        self.b1 = False
        self.b2 = False
        self.b3 = False
        
        self.dat = None
        
    def set_a(self, a):
        self.a = a
        
    def set_str(self, text1, text2):
        self.text = text1
        self._txt = text2
        
    def set_bool(self, b1, b2, b3):
        self.b1 = b1
        self.b2 = b2
        self.b3 = b3
        
    def set_data(self, dat):
        self.dat = numpy.array(dat)
        
    def set_saveable(self, obj):
        self.obj = obj
        
        
    def set_liple(self, liple):
        """Sets a list or tuple 
        """
        self.liple = liple           