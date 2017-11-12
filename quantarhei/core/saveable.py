# -*- coding: utf-8 -*-
"""
    Module saveable

    Defines the class Saveable. Saveable object knows how to inspect
    itself and save its data into a hdf5 file.


"""
import inspect
import fnmatch
import re
import datetime
from numbers import Number

import h5py
import numpy

from .managers import energy_units
from .managers import Manager

def _isattr(obj):
    """Returns True if the object is an attribute of the class, i.e. it is
    not a method.

    """
    ismethod = inspect.ismethod(obj)
    return not ismethod


SIMPLE_TYPES = ("numeric", "strings", "boolean")

class Saveable:
    """Class implementing object persistence through saving to hdf5 files


    FIXME: 
    New Quantarhei issues
    
    1. Multiple saves of the same object should be a avoided
    If the same object is found twice, only a path to it should be saved
       
    2. Block cyclic references of the objects to themselves
    In the first stage, simply through an Exception if cyclic reference is
    encountered. Later, check how such references can be reconstructed
    from the saved information
       
    3. Save version of Quantarhei on the top level of the file
    Something like
    
        def save(self, loc, report_unsaved=False, parent=None):
    
        version = Manager().version
        if parent is None:
            self._save_version(version)
            
        pass
        

    """
    _S__saved_in_this_session_as = None
    _S__stack = []


    def _mark_as_saved(self, loc):
        self._S__saved_in_this_session_as=loc
        
        
    def _saved_in_this_session(self):
        return (self._S__saved_in_this_session_as is not None)

        
    def _save_pointer_instead(self, loc, key):
        """Saves a pointer to the location of this object 
        elsewhere in this file, rather than the object itself
        
        """
        root = self._create_root_group(loc, key)
        root.attrs.create("type",
            numpy.string_("_quantarhei_Saveable_pointer"))
        root.attrs.create("location", self._S__saved_in_this_session_as)

    def _before_save(self):
        """Operations to be done before saving the object

        """
        pass

    def _after_save(self):
        """Operations to be done after saving the object

        """
        pass

    def save(self, file, report_unsaved=False, stack=None, comment=None, 
             test=False):
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
        
        #
        # Check the stack
        #
        sid = id(self)  
        m = Manager()
        
        if stack is not None:
            
            if sid in stack:

                _save_cyclic(file, sid)
                self._after_save()
                return
            
            else:
                
                self._S__stack = stack

        else:
            
            self._S__stack = []
        
        self._S__stack.append(sid) 
        delete_save_dict = False
        
        # if we test or save to a file, the dictionary must be reinitialized
        if isinstance(file, str) or test:
            m.save_dict = {}
            delete_save_dict = True
        
        m.save_dict[sid] = file
        
        #
        # Analyze the class and save it
        #
        strings = {}
        numeric = {}
        boolean = {}
        numdata = {}
        objects = {}
        lists = {}
        tuples = {}
        dictionaries = {}

        attributes = {}

        attr = _get_udef_attributes(self)
        

        nb = 0
        for at_name in attr:
            atr = self.__getattribute__(at_name)
            if atr is None:
                pass # ignoring Nones
            elif isinstance(atr, str):
                strings[at_name] = atr
                # bool but be before Number because apparently bool is Number
            elif isinstance(atr, bool):
                boolean[at_name] = atr
            elif isinstance(atr, Number):
                numeric[at_name] = atr
            elif isinstance(atr, list):
                lists[at_name] = atr
            elif isinstance(atr, tuple):
                tuples[at_name] = atr
            elif isinstance(atr, dict):
                dictionaries[at_name] = atr
            elif isinstance(atr, numpy.ndarray):
                nb += atr.nbytes
                #print("Numdata: ", at_name, id(atr), atr.nbytes/(1024**2))
                numdata[at_name] = atr
            elif isinstance(atr, Saveable):
                objects[at_name] = atr
            else:
                if report_unsaved:
                    print("Attribute: ", at_name, " of ",
                          self.__class__.__name__, "not saved")
                    print("           ", atr)

            attributes["strings"] = strings
            attributes["numeric"] = numeric
            attributes["boolean"] = boolean
            attributes["numdata"] = numdata
            attributes["objects"] = objects
            attributes["lists"] = lists
            attributes["dictionaries"] = dictionaries
            attributes["tuples"] = tuples


        print("Data size: ", nb/(1024**2))
            
        _do_save(self, file=file,
                 attributes=attributes,
                 report_unsaved=report_unsaved,
                 comment=comment, test=test)

        self._S__stack.pop()
        
        #
        # After save clean-up
        #
        self._after_save()

        if delete_save_dict:
            m.save_dict = None

    def _before_load(self):
        """Operations to be done before loading the object

        """
        pass

    def _after_load(self):
        """Operations to be done after loading the object

        """
        pass


    def load(self, file, test=False):
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

        self._before_load()
        
        m = Manager()
        delete_save_dict = False
        
        if isinstance(file, str) or test:
            m.save_dict = {}
            loc = "/"
            delete_save_dict = True
        else:
            loc = file.name
            
        m.save_dict[loc] = self
        
        _do_load(self, file=file,
                 strings=strings,
                 numeric=numeric,
                 boolean=boolean,
                 numdata=numdata,
                 objects=objects,
                 lists=lists,
                 dictionaries=dictionaries,
                 tuples=tuples)

        for key in strings:
            setattr(self, key, strings[key])
        for key in numeric:
            setattr(self, key, numeric[key])
        for key in numdata:
            setattr(self, key, numdata[key])
        for key in boolean:
            setattr(self, key, boolean[key])
        for key in lists:
            setattr(self, key, lists[key])
        for key in objects:
            setattr(self, key, objects[key])
        for key in dictionaries:
            setattr(self, key, dictionaries[key])

        self._after_load()
        
        if delete_save_dict:
            m.save_dict = None

###############################################################################
#
# Utility functions
#
###############################################################################

""" Usage


from quantarhei.saveable import load

m = Molecule()
print(m)
m.save("molecule.hdf5")

m1 = load("molecule.hdf5")
print(m1)

from quantarhei.core.dump import dump_object, load_object

dump_object(m, "molecule.pkl")
m2 = load_object("molecule.pkl")



Loader function

"""

def get_class(fid):
    """Returns class of the saved object
    
    fid is assumed to be the location in the hdf5 file
    
    """

    try:
        k = 0
        for key in fid:
            typ = key
            k += 1
        if k != 1:
            raise Exception("Incorrect format of record")

    except KeyError:
        raise Exception("Object record not recognized")
        
    # get class from string
    cls = _get_class(typ)
    
    return cls
    
    
def load(file, test=False):
    """Loads a Saveable object from a file
    
    
    """
    
    file_opened = False
    
    # if file is a string, we need to open it
    if isinstance(file, str):
        # open file on a disk
        fid = h5py.File(file, "r")
        file_opened = True
    else: 
        fid = file
        
    cls = get_class(fid)
    
    if file_opened or test:
        test = True
        
    # get object from class and load it
    obj = cls()
    obj.load(fid, test=test)
        
    if file_opened:
        fid.close()
        
    # return the loaded object
    return obj

    
def print_info(file):
    """Prints the info of a file
    
    """
    info = read_info(file)
    print("File info on:", file)
    for key in info:
        if key is not "creating_file":
            print(key,": ", info[key])
            

def extract_source_file(file, filename=None):
    """Extracts creating file
    
    """
    info = read_info(file)
    
    fname = info["source_file_name"]
    ftext = info["source_file"]
    
    if filename is not None:
        fname = fname
    if fname is None:
        raise Exception("No file name specified")
        
    if ftext is None:
        print("Source file not recorded: no file created")
        return
        
    with open(fname,"w") as fid:
        print(fid, ftext)
        print("Source file saved as: ", fname)
    

def read_info(file):
    """Reads the information about the file and returns it as a dictionary
    
    """
    
    opened_file = False
    if isinstance(file,str):
        fid = h5py.File(file, "r")
        opened_file = True
    else:
        fid = file
        
    info = {}
    
    k = 0
    for key in fid:
        loc = fid[key]
        #print(key, loc)
        k += 1
    if k != 1:
        raise Exception("Incorrect record format")
        
    info["class"] = loc.attrs["classid"].decode("utf-8")
    info["version"] = loc.attrs["version"].decode("utf-8")
    info["timestamp"] = loc.attrs["timestamp"].decode("utf-8") # includes date
    info["format"] = loc.attrs["format"].decode("utf-8")
    info["filename"] = loc.attrs["filename"].decode("utf-8") # original file name
    info["comment"] = loc.attrs["comment"].decode("utf-8")
#    info["source_file"] = ... # here the whole python file which created this result can be stored
#    info["source_file_name"] = ...
    
    if opened_file:
        fid.close()
    
    return info
        




###############################################################################
#
# helper functions
#
###############################################################################


def _create_root_group(start, name):
    """Creates the root "directory" of the tree to save

    """
    return start.create_group(name)


def _save_cyclic(loc, sid):
    """Saves cyclic reference as a symbolic link
    
    """
    # previous location to which we will make a soft link
    ploc = Manager().save_dict[sid]
    
    # name of the location
    try:
        name = ploc.name
    except AttributeError:
        name = "/"
    #name = ploc.name
#    print(name)
    
    # create soft link a save it in stead of saving the object
    loc.attrs.create("soft_link", numpy.string_(name))    
    
    #raise Exception()
    

def _save_info(loc, info):
    """Saves a dictionary with info

    """    
    for rec in info:
        if rec != "source_file":
            loc.attrs.create(rec, numpy.string_(info[rec]))
        else:
            grp = _create_root_group(loc, "source_file")
            grp.create_dataset("file_content", data=numpy.string_(info[rec]))
            


def _in_simple_types(item):
    """Returns True if the item is of a simple type, i.e. one of the
    string or Number

    """
    if isinstance(item, Number):
        return True
    if isinstance(item, str):
        return True
    return False


def _get_class(kls):
    """Returns a class by its name

    """
    parts = kls.split('.')

    module_name = ".".join(parts[:-1])
    module = __import__(module_name)
    parts = parts[1:]

    for comp in parts:
        module = getattr(module, comp)
    return module


def _get_full_class_name(obj):
    """Returns a full class name of an object

    """
    return obj.__class__.__module__+"."+obj.__class__.__name__


def _get_type(obj):
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


def _is_homogeneous_simple(liple):
    """Returns True and type of the content of the list if the list
    contains objects of the same type and the type is one of the simple
    types: numeric, boolean or strings

    """
    try:
        obj = liple[0]
    except IndexError:
        obj = None

    htyp = _get_type(obj)
    homo = True

    if htyp is None:
        homo = False
        return (homo, htyp)

    for obj in liple:
        if htyp != _get_type(obj):
            htyp = None
            homo = False
            break

    return (homo, htyp)


def _get_udef_attributes(cls):
    """Returns all user defined attributes that are not methods,
    functions, or decorator created properties

    """
    # get all attributes to save
    prog = re.compile(fnmatch.translate("__*__"))

    # members that are properties
    prop_names = []
    for prop in inspect.getmembers(cls.__class__,
                                   lambda o: isinstance(o, property)):
        prop_names.append(prop[0])

    # members that are attributes
    attr = []
    for membs in inspect.getmembers(cls, _isattr):
        if (prog.match(membs[0]) is None) and \
           (membs[0] not in prop_names):
            # attributes starting with _S__ are protected
            # and will not be saved
            if not membs[0].startswith("_S__"):
                attr.append(membs[0])
    return attr


#
# saving functions
#


def _save_classid(loc, classid):
    """Saves a string as "classid" attribute

    """
    fullclassid = classid.__module__+"."+classid.__name__
    loc.attrs.create("classid", numpy.string_(fullclassid))


def _save_a_simple_type(loc, key, item):
    """Saves an object of simple type

    """
    if isinstance(item, Number):
        grp = _create_root_group(loc, key)
        grp.attrs.create("item", item)
        grp.attrs.create("type", numpy.string_("numeric"))
    elif isinstance(item, bool):
        grp = _create_root_group(loc, key)
        grp.attrs.create("item", item)
        grp.attrs.create("type", numpy.string_("boolean"))
    elif isinstance(item, str):
        grp = _create_root_group(loc, key)
        grp.attrs.create("item", numpy.string_(item))
        grp.attrs.create("type", numpy.string_("strings"))


def _save_numdata(loc, dictionary):
    """Saves a dictionary of numpy arrays under the group "numdata"

    """
    if not dictionary:
        return

    strs = _create_root_group(loc, "numdata")
    for key in dictionary.keys():
        strs.create_dataset(key, data=dictionary[key])


def _save_objects(loc, dictionary, stack=None):
    """Saves a dictionary of objects under the group "objects"

    """
    if not dictionary:
        return

    strs = _create_root_group(loc, "objects")
    for key in dictionary.keys():
        obj = dictionary[key]
        nroot = _create_root_group(strs, key)
        obj.save(nroot, stack=stack)


def _save_strings(loc, dictionary):
    """Saves a dictionary of strings under the group "strings"

    """
    if not dictionary:
        return

    strs = _create_root_group(loc, "strings")
    for key in dictionary.keys():
        strs.attrs.create(key, numpy.string_(dictionary[key]))


def _save_numeric(loc, dictionary):
    """Saves a dictionary of numeric values under the group "numeric"

    """
    if not dictionary:
        return

    strs = _create_root_group(loc, "numeric")
    for key in dictionary.keys():
        strs.attrs.create(key, dictionary[key])


def _save_boolean(loc, dictionary):
    """Saves a dictionary of boolean values under the group "boolean"

    """
    if not dictionary:
        return

    strs = _create_root_group(loc, "boolean")
    for key in dictionary.keys():
        strs.attrs.create(key, dictionary[key])


def _save_lists(loc, dictionary, report_unsaved=False, stack=None):
    """Saves a dictionary of lists under the group "lists"

    """
    if not dictionary:
        return

    strs = _create_root_group(loc, "lists")
    for key in dictionary.keys():
        clist = dictionary[key]
        _save_a_listuple(strs, key, "list", clist, report_unsaved, stack=stack)


def _save_tuples(loc, dictionary, report_unsaved=False, stack=None):
    """Saves a dictionary of lists under the group "lists"

    """
    if not dictionary:
        return

    strs = _create_root_group(loc, "tuples")
    for key in dictionary.keys():
        clist = dictionary[key]
        _save_a_listuple(strs, key, "tuple", clist, report_unsaved,
                         stack=stack)


def _save_dictionaries(loc, dictionary, report_unsaved=False, stack=None):
    """Saves a dictionary of dictionaries under the group "dicts"

    """
    if not dictionary:
        return

    strs = _create_root_group(loc, "dicts")
    for key in dictionary.keys():
        cdict = dictionary[key]
        _save_a_dictionary(strs, key, cdict, report_unsaved, stack=stack)


def _save_a_listuple(loc, key, typ, ctuple, report_unsaved=False, stack=None):
    """ Saves a list or tuple of items

    All simple values (strings, numeric, boolean), numpy arrays, Savable
    objects, and also lists, tuples and dictionaries should be acceptable

    """

    if stack is not None:
        stack.append(id(ctuple))
        Manager().save_dict[id(ctuple)] = loc   
        
    root = _create_root_group(loc, key)
    root.attrs.create("type", numpy.string_(typ))

    (homo, htyp) = _is_homogeneous_simple(ctuple)

    root.attrs.create("homogeneous_simple", homo)

    if homo:

        _save_homo_listuple(root, htyp, ctuple)

    else:
        root.attrs.create("length", len(ctuple))

        k = 0
        for item in ctuple:
            tupkey = str(k)

            if _in_simple_types(item):
                _save_a_simple_type(root, tupkey, item)
            elif isinstance(item, Saveable):
                objroot = _create_root_group(root, tupkey)
                objroot.attrs.create("type", numpy.string_("Saveable"))
                item.save(objroot, report_unsaved, stack=stack)
            elif isinstance(item, list):
                _save_a_listuple(root, tupkey, "list", item,
                                 report_unsaved, stack=stack)
            elif isinstance(item, tuple):
                _save_a_listuple(root, tupkey, "tuple", item,
                                 report_unsaved, stack=stack)
            elif isinstance(item, dict):
                _save_a_dictionary(root, tupkey, item, report_unsaved, stack=stack)
            elif item is None:
                pass # Ignoring Nones
            else:
                # possibly warn that something was not saved
                if report_unsaved:
                    print("Unsaved: ", item)
            k += 1

    if stack is not None:
        stack.pop()

def _save_homo_listuple(root, htyp, ctuple):
    """Saves homogenous list into numpy array

    """
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


def _save_a_dictionary(loc, key, cdict, report_unsaved=False, stack=None):
    """Saves a dictionary of values

    All simple values (strings, numeric, boolean), numpy arrays, Savable
    objects, and also lists and dictionaries should be acceptable

    """

    if stack is not None:
        stack.append(id(cdict))
        Manager().save_dict[id(cdict)] = loc

    root = _create_root_group(loc, key)
    root.attrs.create("type", numpy.string_("dict"))

    for dkey in cdict.keys():
        item = cdict[dkey]
        if _in_simple_types(item):
            _save_a_simple_type(root, dkey, item)
        elif isinstance(item, Saveable):
            objroot = _create_root_group(root, dkey)
            objroot.attrs.create("type", numpy.string_("Saveable"))
            item.save(objroot, stack=stack)

        elif isinstance(item, list):
            _save_a_listuple(root, dkey, "list", item, stack=stack)
        elif isinstance(item, tuple):
            _save_a_listuple(root, dkey, "tuple", item, stack=stack)
        elif isinstance(item, dict):
            _save_a_dictionary(root, dkey, item, stack=stack)
        elif item is None:
            pass # ignoring Nones
        else:
            # possibly warn that something was not saved
            if report_unsaved:
                print("Unsaved: ", item)
                
    if stack is not None:
        stack.pop()


def _do_save(cls, file="", attributes=None, report_unsaved=False, comment=None,
             test=False):
    """Performs the save to hdf5 file


    """
    if file == "":
        raise Exception("No file specified")

    file_openned = False

    strings = attributes["strings"]
    numeric = attributes["numeric"]
    boolean = attributes["boolean"]
    numdata = attributes["numdata"]
    objects = attributes["objects"]
    lists = attributes["lists"]
    dictionaries = attributes["dictionaries"]
    tuples = attributes["tuples"]

    
    # check if we should open a file
    fname = ""
    if isinstance(file, str):
        file_openned = True
        fid = h5py.File(file, "w")
        root = _create_root_group(fid, _get_full_class_name(cls))
        fname = file
    else:
        root = _create_root_group(file, _get_full_class_name(cls))

    if comment is not None:
        cmt = comment
    else:
        cmt = ""
        
    if file_openned or test:

        info = {"version":Manager().version,
                "comment":cmt,
                "filename":fname,
                "timestamp":'{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()),
                "format":"Quantarhei_hdf5.ver0.0.1"}
    
        
        _save_info(root, info)

    # do the save of all dictionaries
    with energy_units("int"):
        _save_classid(root, cls.__class__)
        _save_strings(root, strings)
        _save_numeric(root, numeric)
        _save_boolean(root, boolean)
        _save_numdata(root, numdata)
        _save_objects(root, objects, stack=cls._S__stack)
        _save_lists(root, lists, report_unsaved, stack=cls._S__stack)
        _save_tuples(root, tuples, report_unsaved, stack=cls._S__stack)
        _save_dictionaries(root, dictionaries, report_unsaved, stack=cls._S__stack)

    # if file openned here, close it
    if file_openned:
        fid.close()
        


#
# loading functions
#


def _load_classid(loc):
    """Loads a string as "classid" attribute

    """
    try:
        fullclassid = loc.attrs["classid"].decode("utf-8")
    except KeyError:
        fullclassid = ""
        
    return fullclassid


def _load_strings(loc, dictionary):
    """Load a dictionary of strings from the group "strings"

    """
    try:
        strs = loc["strings"]
    except KeyError:
        return

    for key in strs.attrs.keys():
        dictionary[key] = strs.attrs[key].decode("utf-8")


def _load_numeric(loc, dictionary):
    """Saves a dictionary of numeric values under the group "numeric"

    """
    try:
        strs = loc["numeric"]
    except KeyError:
        return

    for key in strs.attrs.keys():
        dictionary[key] = strs.attrs[key]


def _load_boolean(loc, dictionary):
    """Saves a dictionary of boolean values under the group "boolean"

    """
    try:
        strs = loc["boolean"]
    except KeyError:
        return

    for key in strs.attrs.keys():
        dictionary[key] = bool(strs.attrs[key])


def _load_numdata(loc, dictionary):
    """Saves a dictionary of numpy arrays under the group "numdata"

    """
    try:
        strs = loc["numdata"]
    except KeyError:
        return

    for key in strs.keys():
        dictionary[key] = numpy.array(strs[key])


def _load_objects(loc, dictionary, test=False):
    """Saves a dictionary of objects under the group "objects"

    """
    try:
        strs = loc["objects"]
    except KeyError:
        return

    for key in strs.keys():

        objloc = strs[key]
        k = 0
        for kex in objloc.keys():
            objcls = _get_class(kex)
            k += 1
        if k != 1:
            try:
                sftl = objloc.attrs["soft_link"].decode("utf-8")
            except KeyError:
                raise Exception("Empty object")

            obj = Manager().save_dict[sftl]
            
        else:
            obj = objcls()
            # we save the location of a newly created object
            # and save it under the current location

            Manager().save_dict[objloc.name] = obj
            obj.load(objloc, test=test)
            
        dictionary[key] = obj


def _load_lists(loc, dictionary, test=False):
    """Saves a dictionary of lists under the group "lists"

    """
    try:
        strs = loc["lists"]
    except KeyError:
        return

    for key in  strs.keys():
        loc = strs[key]
        clist = _load_a_listuple(loc, test=test)
        dictionary[key] = clist


def _load_dictionaries(loc, dictionary, test=False):
    """Saves a dictionary of dictionaries under the group "dicts"

    """
    try:
        strs = loc["dicts"]
    except KeyError:
        return

    for key in strs.keys():
        loc = strs[key]
        cdict = _load_a_dictionary(loc, test=test)
        dictionary[key] = cdict


def _load_tuples(loc, dictionary, test=False):
    """Saves a dictionary of dictionaries under the group "dicts"

    """
    try:
        strs = loc["tuples"]
    except KeyError:
        return

    for key in  strs.keys():
        loc = strs[key]
        clist = _load_a_listuple(loc, test=test)
        dictionary[key] = clist


def _load_a_simple_type(loc):
    """Loads any of the simple types

    """

    typ = loc.attrs["type"]
    typ = typ.decode("utf-8")
    if typ == "numeric":
        return loc.attrs["item"]
    if typ == "strings":
        sstr = loc.attrs["item"]
        return sstr.decode("utf-8")
    if typ == "boolean":
        return loc.attrs["item"]


def _load_a_listuple(loc, test=False):
    """Loads a list or tuple of values

    All simple values (strings, numeric, boolean), numpy arrays, Savable
    objects, and also lists and dictionaries should be acceptable

    """

    typ = loc.attrs["type"].decode("utf-8")
    homo = loc.attrs["homogeneous_simple"]

    clist = []

    if homo:

        _load_homo_listuple(loc, clist)

    else:

        N = loc.attrs["length"]
        for k in range(N):
            tupkey = str(k)

            read_item = False
            try:
                itemloc = loc[tupkey]
                read_item = True
            except KeyError:
                #print(tupkey, "does not exist")
                pass

            if read_item:

                _read_item(itemloc, clist, test=test)

            else:

                # if the object in the list does not exist, it is replaced
                # by None
                clist.append(None)

            k += 1

    # if the type is tuple, we convert list to tuple
    if typ == "tuple":
        #print("convering to tuple")
        clist = tuple(clist)

    return clist


def _load_homo_listuple(loc, clist):
    """Loads a homogenous list from a numpy array

    """
    if "numeric" in loc.keys():
        data = loc["numeric"]
    elif "boolean" in loc.keys():
        data = loc["boolean"]
    elif "strings" in loc.keys():
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


def _read_item(itemloc, clist, test=False):
    """Reads one item from a location

    """
    ltyp = itemloc.attrs["type"].decode("utf-8")
    if ltyp in SIMPLE_TYPES:
        item = _load_a_simple_type(itemloc)
        clist.append(item)
    elif ltyp == "Saveable":

        k = 0
        for key in itemloc.keys():
            classname = key
            k += 1
            
        if k != 1:
            slink = itemloc.attrs["soft_link"].decode("utf-8")
            obj = Manager().save_dict[slink]
        else:
            cls = _get_class(classname)
            obj = cls()
            Manager().save_dict[itemloc.name] = obj
            obj.load(itemloc, test=test)
        
        clist.append(obj)

    elif ltyp == "list":
        item = _load_a_listuple(itemloc, test=test)
        clist.append(item)
    elif ltyp == "tuple":
        item = _load_a_listuple(itemloc, test=test)
        clist.append(item)
    elif ltyp == "dict":
        item = _load_a_dictionary(itemloc, test=test)
        clist.append(item)
    else:
        # possibly warn that something was not saved
        pass


def _load_a_dictionary(loc, test=False):
    """Loads a dictionary of values

    All simple values (strings, numeric, boolean), numpy arrays, Savable
    objects, and also lists and dictionaries should be acceptable

    """
    cdict = {}

    for dkey in loc.keys():

        read_item = False
        try:
            itemloc = loc[dkey]
            read_item = True
        except KeyError:
            #print(dkey, "does not exist")
            pass

        if read_item:

            ltyp = itemloc.attrs["type"].decode("utf-8")

            if ltyp in SIMPLE_TYPES:
                item = _load_a_simple_type(itemloc)
                cdict[dkey] = item
            elif ltyp == "Saveable":

                k = 0
                for key in itemloc.keys():
                    classname = key
                    k += 1
                    
                if k != 1:
                    slink = itemloc.attrs["soft_link"].decode("utf-8")
                    obj = Manager().save_dict[slink]
                
                else:

                    cls = _get_class(classname)
                    obj = cls()
                    Manager().save_dict[itemloc.name] = obj
                    obj.load(itemloc, test=test)
                cdict[dkey] = obj

            elif ltyp == "list":
                item = _load_a_listuple(itemloc, test=test)
                cdict[dkey] = item
            elif ltyp == "tuple":
                item = _load_a_listuple(itemloc, test=test)
                cdict[dkey] = item
            elif ltyp == "dict":
                item = _load_a_dictionary(itemloc, test=test)
                cdict[dkey] = item
            else:
                # possibly warn that something was not saved
                pass

        else:

            cdict[dkey] = None

    return cdict


def _do_load(cls, file="", strings=None, numeric=None,
             boolean=None, numdata=None, objects=None,
             lists=None, dictionaries=None, tuples=None,
             test=False):
    """Loads data from hdf5 file


    """
    if file == "":
        raise Exception("No file specified")

    file_openned = False
    thisclass = _get_full_class_name(cls)

    # check if we should open a file
    if isinstance(file, str):
        file_openned = True
        fid = h5py.File(file, "r")
        root = fid[thisclass]

    else:
           
        root = file[thisclass]
            
    cid = _load_classid(root)

    if cid != thisclass:
        raise Exception("File contains an object of unexpected type")

    _load_strings(root, strings)
    _load_numeric(root, numeric)
    _load_boolean(root, boolean)
    _load_numdata(root, numdata)
    _load_objects(root, objects, test=test)
    _load_lists(root, lists, test=test)
    _load_tuples(root, tuples, test=test)
    _load_dictionaries(root, dictionaries, test=test)

    # if file openned here, close it
    if file_openned:
        fid.close()
