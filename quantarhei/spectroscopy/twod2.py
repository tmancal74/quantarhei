# -*- coding: utf-8 -*-
"""Two-dimensional Fourier Transform Spectrum and its calculation


"""
from functools import partial
import numbers
import matplotlib.pyplot as plt  
import numpy

from ..core.frequency import FrequencyAxis
from .. import COMPLEX
from .. import signal_REPH, signal_NONR, signal_TOTL, signal_DC
from .. import TWOD_SIGNALS
from .. import part_REAL, part_IMAGINARY, part_ABS
from ..core.saveable import Saveable
from ..core.datasaveable import DataSaveable
from ..core.dfunction import DFunction
from ..core.valueaxis import ValueAxis
from ..utils.types import check_numpy_array
from .twod import TwoDSpectrum

# FIXME: Check these names

#
#  Pathway types
#
_ptypes = ["R1g", "R2g", "R3g", "R4g", "R1fs", "R2fs", "R3fs", "R4fs"]

#
# Processes --- GSB, SE, ESA and DC
#
_processes = dict(GSB=[_ptypes[0], _ptypes[1]], SE=[_ptypes[2], _ptypes[3]],
                  ESA=[_ptypes[4], _ptypes[5]], DC=[_ptypes[6], _ptypes[7]])

#
# Types of signals --- rephasing (REPH), non-rephasing (NONR) 
#                      and double coherence (DC)
#
_signals = {signal_REPH:[_ptypes[1], _ptypes[2], _ptypes[4]],
            signal_NONR:[_ptypes[0], _ptypes[3], _ptypes[5]],
            signal_DC:[_ptypes[6], _ptypes[7]]}

_total = signal_TOTL

#
# Storage resolutions
#
_resolutions = ["off", "signals", "processes", "types", "pathways"]


def _resolution2number(res):
    """Converts resolution string to number
    
    Parameters
    ----------
    
    res : string
        Resolution string. Here is conversion table    
        | string    | integer |
        | pathways  |   4     |
        | types     |   3     |
        | processes |   2     |
        | signals   |   1     |
        | off       |   0     |
    
    
    >>> _resolution2number("pathways")
    4
    
    >>> _resolution2number("types")
    3
    
    >>> _resolution2number("processes")
    2
    
    >>> _resolution2number("signals")
    1
    
    >>> _resolution2number("off")
    0
    
    
    
    """
    if res in _resolutions:
        return _resolutions.index(res)
    else:
        raise Exception("Unknow resolution level in TwoDSpectrum")


def _get_type_and_tag(obj, storage):

    if obj.current_dtype not in _ptypes:
        # check the current_type attribute
        raise Exception("Wrong pathways type")

    try:
        # get the dictionary of pathways with a give type
        piece = storage[obj.current_dtype]
    except IndexError:
        # if it does not exist, create it
        storage[obj.current_dtype] = {}
        piece = storage[obj.current_dtype]

    if obj.current_tag in piece.keys():
        # if the tag exists raise Exception
        raise Exception("Tag already exists")
        
    return piece


def _pathways_to_processes(obj, process):

    if obj.storage_initialized:
        data = numpy.zeros((obj.xaxis.length,
                            obj.yaxis.length),
                            dtype=COMPLEX)
    else:
        return None #numpy.zeros((1,1), dtype=COMPLEX)
        
    if process not in _processes:
        raise Exception("Unknown process: "+process)
        
    else:
        
        # types corresponding to process
        types = _processes[process]
        for typ in types:
            # pways corresponding to a given type
            try:
                pways = obj._d__data[typ]
            except:
                pways = []
            # sum those pathways
            for tag in pways:
                data += pways[tag]
        
    return data


def _pathways_to_signals(obj, signal):

    if obj.storage_initialized:
        data = numpy.zeros((obj.xaxis.length,
                            obj.yaxis.length),
                            dtype=COMPLEX)
    else:
        data = None #numpy.zeros((1,1), dtype=COMPLEX)
        
    if signal not in _signals:
        raise Exception("Unknown signal: "+signal)
        
    else:
        
        # types corresponding to signal
        types = _signals[signal]
        for typ in types:
            # pways corresponding to a given type
            try:
                pways = obj._d__data[typ]
            except:
                pways = []
            # sum those pathways
            for tag in pways:
                data += pways[tag]
        
    return data        


def _pathways_to_total(obj):
    
    if obj.storage_initialized:
        data = numpy.zeros((obj.xaxis.length,
                            obj.yaxis.length),
                            dtype=COMPLEX)
    else:
        data = None #numpy.zeros((1,1), dtype=COMPLEX)

    for signal in _signals:
        
        data += _pathways_to_signals(obj, signal)

    return data        


def _types_to_processes(obj, process):
    """Sums pathways of different types into a specified process spectrum
    
    
    """
    if obj.storage_initialized:
        data = numpy.zeros((obj.xaxis.length,
                            obj.yaxis.length),
                            dtype=COMPLEX)
    else:
        data = None #numpy.zeros((1,1), dtype=COMPLEX)


    types = _processes[process]
    for dtype in types:
        try:
            ddata = obj._d__data[dtype]
        except KeyError:
            # set to None if dtype not present
            ddata = None
        except AttributeError:
            # no data
            ddata = None
            
        if ddata is not None:
            if data is not None:
                data += ddata
            else:
                data = ddata

    return data


def _types_to_signals(obj, signal):
    """Sums pathways of different types into a specified signal spectrum
    
    
    """
    if obj.storage_initialized:
        data = numpy.zeros((obj.xaxis.length,
                            obj.yaxis.length),
                            dtype=COMPLEX)
    else:
        data = None #numpy.zeros((1,1), dtype=COMPLEX)
        
    types = _signals[signal]
    
    for dtype in types:
        try:
            ddata = obj._d__data[dtype]
        except KeyError:
            # set to None if dtype not present
            ddata = None
        except AttributeError:
            # no data
            ddata = None
            
        if ddata is not None:
            if data is not None:
                data += ddata
            else:
                data = ddata

    return data


def _signals_to_total(obj):
    """Sums spectra corresponding to different signals into the total spectrum


    """    
    if obj.storage_initialized:
        data = numpy.zeros((obj.xaxis.length,
                            obj.yaxis.length),
                            dtype=COMPLEX)
    
        for signal in _signals:
            
            try:
                data += obj._d__data[signal]
            except KeyError:
                pass
                # likely the corresponding spectrum is not defined, this may be
            
    else:
        data = None #numpy.zeros((1,1), dtype=COMPLEX)

    return data
            

def _processes_to_total(obj):
    """Sums spectra corresponding to different processes into the total spectrum


    """
    
    if obj.storage_initialized:
        data = numpy.zeros((obj.xaxis.length,
                            obj.yaxis.length),
                            dtype=COMPLEX)
    
        for process in _processes:
            
            try:
                data += obj._d__data[process]
            except KeyError:
                pass
                # likely the corresponding spectrum is not defined, this my be            
            
    else:
        data = None #numpy.zeros((1,1), dtype=COMPLEX)

    return data


def _types_to_total(obj):
    """Sums all pathways into the total spectrum
    
    
    """
    if obj.storage_initialized:
        data = numpy.zeros((obj.xaxis.length,
                            obj.yaxis.length),
                            dtype=COMPLEX)
    
        # sum over all processes
        for process in _processes:
            
            data += _types_to_processes(obj, process)
            
    else:
        data = None #numpy.zeros((1,1), dtype=COMPLEX)

    return data


def twodspectrum_dictionary(name, dtype):
    """Defines operations of setting and retrieving (getting) data 
    
    Operations on the storage of two-dimensional spectral data
    
    
    """
    
    storage_name = '_'+name

    @property
    def prop(self):
        
        storage = getattr(self, storage_name)
        
        #
        # with pathway resolution => type and tag has to be specified
        #
        if self.storage_resolution == "pathways":
             
            if self.current_dtype in _ptypes:

                if self.storage_initialized:
    
                    try:
                        piece = storage[self.current_dtype]
                    except:
                        return numpy.zeros((self.xaxis.length,
                                            self.yaxis.length),
                                            dtype=COMPLEX)

                else:
                    return None #numpy.zeros((1,1), dtype=COMPLEX)
                                
                #
                # return as pathway
                #
                if self.current_tag is not None:
                    return piece[self.current_tag]
            
                #
                # return as type
                #
                else:
                    # tag not specified so we add up all pathways
                    # of a given type
                    k_i = 0
                    for tag in piece:
                        dat = piece[tag]
                        if k_i == 0:
                            data = dat.copy()
                        else:
                            data += dat
                        k_i += 1
                            
                    return data
                    
            elif self.current_dtype in _processes:
                
                #
                # return as process
                #
                return _pathways_to_processes(self, self.current_dtype)
                
            elif self.current_dtype in _signals:
                
                #
                # return as signals
                #
                return _pathways_to_signals(self, self.current_dtype)
                
            elif self.current_dtype == _total:
                
                #
                # return total spectrum
                #
                return _pathways_to_total(self)
            
            else:
                raise Exception("Inappropriate data type: "
                                +self.current_dtype+" for storage resolution "
                                +self.storage_resolution)
 
        #
        # Resolution = "types"
        #
        elif self.storage_resolution == "types":
            
            if self.current_dtype in _ptypes:
                try: 
                    ret = storage[self.current_dtype]
                except KeyError:
                    ret = None
                    
                return ret

            elif self.current_dtype in _processes:
                
                # return as process
                return _types_to_processes(self, self.current_dtype)
            
            elif self.current_dtype in _signals:
                #
                # return as signals
                #
                return _types_to_signals(self, self.current_dtype)
                
            elif self.current_dtype == _total:
                
                #
                # return total spectrum
                #
                return _types_to_total(self)

            else:
                raise Exception("Inappropriate data type: "
                                +self.current_dtype+" for storage resolution "
                                +self.storage_resolution)

        #
        # Resolution = "processes"
        #
        elif self.storage_resolution == "processes":
            
            if self.current_dtype in _processes:
                try:
                    ret = storage[self.current_dtype]
                except KeyError:
                    ret = None
                    
                return ret

            elif self.current_dtype == _total:
                
                #
                # Return total spectrum
                #
                return _processes_to_total(self)

            else:
                raise Exception("Inappropriate data type: "
                                +self.current_dtype+" for storage resolution "
                                +self.storage_resolution)

        #
        # Resolution = "signals"
        #
        elif self.storage_resolution == "signals":

            if self.current_dtype in _signals:
                try:
                    ret = storage[self.current_dtype]
                except KeyError:
                    ret = None
                    
                return ret
            
            elif self.current_dtype == _total:
                
                #
                # return total spectrum
                #
                return _signals_to_total(self)                

            else:
                raise Exception("Inappropriate data type: "
                                +self.current_dtype+" for storage resolution "
                                +self.storage_resolution)
                
        elif self.storage_resolution == "off":
            
            if self.current_dtype == _total:
                try:
                    ret = storage[_total]
                except KeyError:
                    ret = None
                    
                return ret

            else:
                raise Exception("Inappropriate data type: "
                                +self.current_dtype+" for storage resolution "
                                +self.storage_resolution)

        else:
            raise Exception("not implemented")              

        
    @prop.setter
    def prop(self, value):
        
        ini = getattr(self, "storage_initialized")
        if not ini:
            setattr(self, storage_name, {})
            setattr(self, "storage_initialized", True)
        
        if isinstance(value, numpy.ndarray):

            storage = getattr(self, storage_name)

            #
            # with pathway resolution => type and tag has to be specified
            #
            if self.storage_resolution == "pathways":

                if self.current_dtype not in _ptypes:
                    # check the current_type attribute
                    raise Exception("Wrong pathways type")
            
                try:
                    # get the dictionary of pathways with a give type
                    piece = storage[self.current_dtype]
                except KeyError:
                    # if it does not exist, create it
                    storage[self.current_dtype] = {}
                    piece = storage[self.current_dtype]
            
                if self.current_tag in piece.keys():
                    # if the tag exists raise Exception
                    raise Exception("Tag "+self.current_tag+" already exists")
                        
                    if value.shape != (self.xaxis.length, self.yaxis.length):
                        # if the data shape is not consistent, raise Exception
                        raise Exception("Data not consistent "+
                                        "with spectrum axes")

                piece[self.current_tag] = value

            #
            # pathway type resolution
            #
            elif self.storage_resolution == "types":
                
                if self.current_dtype not in _ptypes:
                    # check the current_type attribute
                    raise Exception("Wrong pathways type: "+self.current_dtype)
                    
                storage[self.current_dtype] = value  

            #
            # signal type resolution
            #
            elif self.storage_resolution == "signals":

                if self.current_dtype not in _signals:
                    # check the current_type attribute
                    raise Exception("Wrong signal type: "+self.current_dtype)
                    
                storage[self.current_dtype] = value  

            #
            # signal type resolution
            #
            elif self.storage_resolution == "processes":

                if self.current_dtype not in _processes:
                    # check the current_type attribute
                    raise Exception("Wrong process type: "+self.current_dtype)
                    
                storage[self.current_dtype] = value

            #
            # No resolution
            #
            elif self.storage_resolution == "off":

                if self.current_dtype != _total:
                    # check the current_type attribute
                    raise Exception("Wrong data type: "++self.current_dtype)
                    
                storage[_total] = value
                    
                                        
        else:
            raise TypeError('{} must contain \
                            values of type {})'.format(name, dtype), dtype)
            
    return prop

#
# Storage type for 2D spectra
#
TwoDSpectrumDataArray = partial(twodspectrum_dictionary,
                                dtype=numbers.Complex)     


def twod_data_wrapper(dtype):
    """Handles access to the d__data property of TwoDSpectrum 
    
    """
    storage_name = "d__data"
    
    @property
    def prop(self):
        return getattr(self,storage_name)
    
    @prop.setter
    def prop(self, value):
        if self._allow_data_writing:
            try:
                vl = check_numpy_array(value) 
                setattr(self,storage_name,vl)
            except:
                raise TypeError(
                'Value must be either a list or numpy.array')
        else:  
            raise Exception("`data` property is for reading only")
        
    return prop


DataWrapper = partial(twod_data_wrapper, dtype=numbers.Complex)


class TwoDSpectrumBase(DataSaveable):
    """Basic class of a two-dimensional spectrum
    
    
    """
    
    # spectral types
    stypes = TWOD_SIGNALS #["rephasing", "non-rephasing", "nonrephasing", "total"]
    
    # spectral parts
    sparts = [part_REAL, part_IMAGINARY] #["real", "imaginary"]
    
    # to keep Liouville pathways separate?
    keep_pathways = False
    
    dtypes = TWOD_SIGNALS #["Tot", "Reph", "Nonr"]
    
    # to keep stypes separate?
    keep_stypes = True

    #
    # Storage of 2D data
    #
    d__data = TwoDSpectrumDataArray("d__data")
    data = DataWrapper()

    _allow_data_writing = False
    
    def __init__(self):
        super().__init__()
        
        self.xaxis = None
        self.yaxis = None
            
        #
        # The followinh attributes will be removed 
        ##########################################
        self.reph2D = None
        self.nonr2D = None
        #self.data = None
        
        self.dtype = None
        ##########################################
        
        
        # initially, the highest possible resolution is set
        self.storage_resolution = "pathways"
        self.storage_initialized = False
        
        # if no dtype is specified, we return total spectrum
        self.current_dtype = signal_TOTL #"total"  
        self.current_tag = None
        self.address_length = 1


    def set_data_writable(self):
        """Lifts the protection of the data property 
        
        """
        if not self.storage_initialized:
            self.set_resolution("off")
            
        self._allow_data_writing = True

    def set_data_protected(self):
        """Sets back the protection of the data property 
        
        """
        self._allow_data_writing = False
        
    def set_axis_1(self, axis):
        """Sets the x-axis of te spectrum (omega_1 axis)
        
        """
        self.xaxis = axis


    def set_axis_3(self, axis):
        """Sets the y-axis of te spectrum (omega_3 axis)
        
        """
        self.yaxis = axis


#    def initialize_storage(self, storage_resolution):
#        self.storage_resolution = storage_resolution
        

    # FIXME: to be deprecated
    def set_data_type(self, dtype=signal_TOTL): #"Tot"):
        """Sets the data type for this 2D spectrum
        
        Parameters
        ----------
        
        dtype : string
           Specifies the type of data stored in this TwoDSpectrum object 
        """
        
        if dtype in self.dtypes:
            self.dtype = dtype
        else:
            raise Exception("Unknown data type for TwoDSpectrum object")

    # FIXME: to be replaced       
    def set_data(self, data, dtype=signal_TOTL): #"Tot"):
        """Sets the data of the 2D spectrum
        
        Sets the object data depending on the specified type and stores the
        type
        
        
        Parameters
        ----------
        
        data : 2D array
            Data of the spectrum, float or complex
            
        dtype : string
            Type of the data stored. Three values are allowed: if quantarhei
            is import as `qr` that they are qr.signal_REPH
            for rephasing spectra, qr.signal_NONR for non-rephasing spectra,
            and qr.signal_TOTL for total spectrum, which is the sum of both
            
        """
        if self.dtype is None:
            if dtype in self.dtypes:
                self.dtype = dtype
        else:
            if dtype != self.dtype:
                raise Exception("Incorrect data type in TwoDSpectrum")
                
        if dtype == signal_TOTL:
            self.data = data
            
        elif dtype == signal_REPH:
            self.reph2D = data
            
        elif dtype == signal_NONR:        
            self.nonr2D = data
            
        else:
            
            raise Exception("Unknow type of data: "+dtype)

    # FIMXE: to be relaced
    def add_data(self, data, dtype=signal_TOTL): #"Tot"):
        
        if dtype is None:
            if dtype in self.dtypes:
                self.dtype = dtype
        else:
            if dtype != self.dtype:
                raise Exception("Incorrect data type in TwoDSpectrum")

        if dtype == signal_TOTL:
            
            if self.data is None:
                self.data = numpy.zeros(data.shape, dtype=data.dtype)
            self.data += data
            
        elif dtype == signal_REPH:
            if self.reph2D is None:
                self.reph2D = numpy.zeros(data.shape, dtype=data.dtype)
            self.reph2D += data                
            
        elif dtype == signal_NONR:
            if self.nonr2D is None:
                self.nonr2D = numpy.zeros(data.shape, dtype=data.dtype)                
            self.nonr2D += data

        else:
            
            raise Exception("Unknow type of data: "+dtype)


    def set_data_flag(self, flag):
        """Sets a flag by which date will be retrieved and save to `data`
        attribute
        
        
        """
        
        if isinstance(flag, list):
            self.current_dtype = flag[0]
            try:
                self.current_tag = flag[1]
            except IndexError:
                raise Exception("flag in form of a list must have"
                                +" two elements")
            self.address_length = 2
        else:
            self.current_dtype = flag
            self.current_tag = None
            self.address_length = 1


    def get_all_tags(self):
        """Retunrs tags of the pathways stored under `pathways` resolution
        
        The tags are returned in a two membered list with its type, i.e.
        as [type, tag]
        
        """
        
        if self.storage_resolution != "pathways":
            return []
        
        tags = []
        for typ in _ptypes:
            try:
                pdict = self._d__data[typ]
            except KeyError:
                pdict = {}
            keys = pdict.keys()
            for key in keys:
                tags.append([typ, key])
                
        return tags


    # FIXME: maybe set_storage_resolution ?
    def set_resolution(self, resolution):
        """Sets the storage resolution attribute of TwoDSpectrum


        Parameters
        ----------
        
        resolution : string
            Resolution in which data are stored in TwoDSpectrum object.
            Values are one of the strings: `pathways` - stores individual
            Liouville pathways, `types` - stores sets of pathways corresponding
            to different shapes of Feynman diagrams, `processes` - stores
            data corresponding to processes, such as stimulated emission and
            ground state bleach, `signals` - stores only rephasing and
            non-rephasing signals separately, `off` - stores only total
            spectrum.
            
        
        Examples
        --------
        
        Initial resolution is the highest one, i.e. 'pathways'
        
        >>> spect1 = TwoDResponse()
        >>> spect2 = TwoDResponse()
        >>> spect1.storage_resolution
        'pathways'
        
        We can set only decreasing resolution
        
        >>> spect1.set_resolution("types")
        >>> spect2.set_resolution("types")
        >>> spect1.storage_resolution
        'types'

        "types" can be converted either to "processes" or "signals"
        
        >>> spect1.set_resolution("processes")
        >>> spect1.storage_resolution
        'processes'

        "processes" cannot be converted to "signals"

        >>> spect1.set_resolution("signals")
        Traceback (most recent call last):
            ...
        Exception: Cannot convert resolution for level 2 to level 1
                
        >>> spect2.set_resolution("signals")
        >>> spect2.storage_resolution
        'signals'
        
        If we set increasing resolution we get an exception
        
        >>> spect1.set_resolution("types")
        Traceback (most recent call last):
            ...
        Exception: Cannot convert from lower to higher resolution

        From "signals" and "types" you can convert to no resolution, i.e.
        to total spectrum where storage resolution is 'off'

        >>> spect1.set_resolution("off")
        >>> spect1.storage_resolution
        'off'
        
        >>> spect2.set_resolution("off")
        >>> spect2.storage_resolution
        'off'
        
        """
        if resolution in _resolutions:
            res_old = _resolution2number(self.storage_resolution)
            res_new = _resolution2number(resolution)
            if res_old < res_new:
                raise Exception("Cannot convert from lower"+
                                " to higher resolution")
            elif res_old > res_new:
                # recalculate data towards lower resolution
                self._convert_resolution(res_old, res_new)
            
            #self.storage_resolution = resolution
        else:
            raise Exception("Unknown resolution: "+resolution)


    # FIXME: maybe get_storage_resolution
    def get_resolution(self):
        """Returns storage resolution
        
        """
        return self.storage_resolution


    def _convert_resolution(self, old, new):
        """Converts storage from one level of resolution to another
        
        """
        
        _conversion_paths = {4:{3:[4,3], 2:[4,3,2], 1:[4,3,1], 0:[4,3,2,0]}, 
                 3:{2:[3,2], 1:[3,1], 0:[3,2,0]},
                 2:{0:[2,0]},
                 1:{0:[1,0]}}
        
        try:
            path = _conversion_paths[old][new]
        except KeyError:
            raise Exception("Cannot convert resolution for level "+str(old)+
                            " to level "+str(new))
        
        for step in path:
            if step == old:
                start = old
            else:
                end = step
                self._convert_res_elementary(start, end)
                self.storage_resolution = _resolutions[end]
                start = end

        
    def _convert_res_elementary(self, old, new):
        """Performs simple storage conversions involving only one step conversion
        
        """
        
        # convert "pathways" to "types"
        if (old == 4) and (new == 3):
            
            storage = {}
            data_present = True
            
            for dtype in _ptypes:
                # get a dictionary of pathways of a given type
                try:
                    pdict = self._d__data[dtype]
                except KeyError:
                    # ignore if some are absent
                    pdict = {}                    
                except AttributeError:
                    # no data
                    pdict = {}
                    data_present = False

                # sum all data
                if data_present:
                    data = numpy.zeros((self.xaxis.length, self.yaxis.length),
                                       dtype=COMPLEX)
                else:
                    data = numpy.zeros((1, 1),
                                       dtype=COMPLEX)
                    
                for key in pdict.keys():
                    data += pdict[key]
                
                storage[dtype] = data
           
            self._d__data = storage
                  
        # convert "types" to "processes"
        elif (old == 3) and (new == 2):
            
            storage = {}
            
            for process in _processes.keys():
                
                data = _types_to_processes(self, process)                    
                storage[process] = data
           
            self._d__data = storage
            
        # convert "types" to "signals"
        elif (old == 3) and (new == 1):

            storage = {}
            
            for signal in _signals.keys():
                
                data = _types_to_signals(self, signal) 
                storage[signal] = data
           
            self._d__data = storage
            
        # converts "signals" to "off"
        elif (old == 1) and (new == 0):
            storage = {}
            
            data = _signals_to_total(self)
            storage[_total] = data
            
            self._d__data = storage
        
        # converts "processes" to "off"
        elif (old == 2) and (new == 0):
            storage = {}
            
            data = _processes_to_total(self)
            storage[_total] = data
            
            self._d__data = storage
            
        else:
            raise Exception("Cannot convert resolution for level "+str(old)+
                            " to level "+str(new))


    # FIXME: this will become the main add_data method
    def _add_data(self, data, resolution=None, dtype=_total, tag=None):
        """Adds data to this 2D spectrum
        
        This method is used when partial data are stored in this spectrum
        object, and it is expected that more data will come. To set spectrum
        data in one shot, you can use `set_data` method.
        
        
        Parameters
        ----------
        
        data : array
            Numpy array compatible in dimensions with the axes of the spectrum
            
        resolution : string or None
            Resolution of adding data. If the data correspond to individual
            `pathways` (resolution="pathways"), a `types` of pathways
            (resolution="types") such as "R1g", "R2g" etc., 
            `process` (resolution="processes") such as "GSB", "ESA" etc.,
            or `signals` such as "Reph" or "Nonr" (resolution="signals"). 
            One can also store a complete spectrum under the resolution="off".
            
        dtype : string
            Type of data; under resolution="pathway", dtype specifies the 
            character of the pathway (such as "R1g", "R2g", etc.)
            
        tag : string
            Used in the resolution="pathway". It provides a unique tag to 
            identify the pathway
            
        """
        if not self.storage_initialized:
            self._d__data = {}
            self.storage_initialized =  True
            if resolution is not None:
                self.storage_resolution = resolution
            
        if resolution is None:
            resolution = self.storage_resolution
        else:
            res1 = _resolution2number(resolution)
            res2 = _resolution2number(self.storage_resolution)
            if res1 <= res2:
            
                pass
            
            else:
                raise Exception("This TwoDSpectrum does not have enough "
                                +"resolution to add data with resolution = "
                                +resolution)

        if resolution == "pathways":
            if dtype in _ptypes:
                if tag is not None:
                    self.set_data_flag([dtype, tag])
                    try:
                        odata = self.d__data
                    except:
                        odata = None
                        
                    if odata is None:    
                        self.d__data = data
                    else:
                        self.d__data = odata + data
                else:
                    raise Exception("Tag for Liouville pathway not specified")
            else:
                raise Exception("Storage resolution 'pathways': "+
                                "Unknown type of Liouville pathway: "+dtype)
                
        elif resolution == "types":
            if dtype in _ptypes:
                if tag is not None:
                    raise Exception("Tag specified for storage resolutios"+
                                    " 'types'. Tag would be ignored and"+
                                    " information lost")
                self.set_data_flag(dtype)
                try:
                    odata = self.d__data
                except:
                    odata = None
                    
                if odata is None:
                    self.d__data = data
                else:
                    self.d__data = odata + data
            
            else:
                raise Exception("Storage resolution 'types': "+
                                "Unknown type of Liouville pathway: "+dtype)


        elif resolution == "processes":
            
            if dtype in _processes:
                if tag is not None:
                    raise Exception("Tag specified for storage resolutios"+
                                    " 'processes'. Tag would be ignored and"+
                                    " information lost")
                self.set_data_flag(dtype)
                try:
                    odata = self.d__data
                except:
                    odata = None
                    
                if odata is None:
                    self.d__data = data
                else:
                    self.d__data = odata + data
            
            else:
                raise Exception("Storage resolution 'processes': "+
                                "Unknown type of signal: "+dtype)
        
        elif resolution == "signals":
            
            if dtype in _signals:
                if tag is not None:
                    raise Exception("Tag specified for storage resolutios"+
                                    " 'signals'. Tag would be ignored and"+
                                    " information lost")
                self.set_data_flag(dtype)
                try:
                    odata = self.d__data
                except:
                    odata = None
                    
                if odata is None:
                    self.d__data = data
                else:
                    self.d__data = odata + data
            
            else:
                raise Exception("Storage resolution 'signals': "+
                                "Unknown type of signal: "+dtype)
        
        elif resolution ==  "off":
            
            if dtype == _total:
                if tag is not None:
                    raise Exception("Tag specified for storage resolutios"+
                                    " 'off'. Tag would be ignored and"+
                                    " information lost")
                self.set_data_flag(dtype)
                try:
                    odata = self.d__data
                except:
                    odata = None
                    
                if odata is None:
                    self.d__data = data
                else:
                    self.d__data = odata + data
                
            else:
                raise Exception("Storage has no resolution. "+
                                "Only total spectrum can be saved. "
                                +"Used data type: "+dtype)
                

    # FIXME: implement
    def _set_data(self, data, dtype=_total, tag=None):
        """Sets the 2D spectrum data
        
        Depending on the storage resolution inferred from the combination
        of the data type and tag, we store the data, and set the corresponding
        storage resolution flag.
        
        Parameters
        ----------
        
        data : numpy array
            Data of the spectrum
            
        dtype : string
            Data type: "total", "REPH", "NONR", "R1g", "R2g", etc.
            
        tag : string
            If a type of pathway is specified ("R1g", "R2g", etc.), tag which
            is not None allows storage of individual pathways
            
        """
        
        pass
    
    
    def get_all_data(self):
        """Returns a dictionary of all data stored in the object
        
        The data are returned as a dictionary with keys corresponding
        to current storage mode. If the storage resolution is `pathway`,
        keys are strings constucted as TYPE+"_"+str(TAG), where TYPE is the
        type string and TAG is the tag of the pathway.
        
        It is expected that most often the keys will not be even used
        
        """
        data_dict = {}
        if self.storage_resolution == "pathways":
            for typ in _ptypes:
                try:
                    piece = self._d__data[typ]
                except KeyError:
                    piece = []
                for tag in piece:
                    key = typ+"_"+str(tag)
                    data_dict[key] = piece[tag]
        
        elif self.storage_resolution == "types":
            for typ in _ptypes:
                data_ex = True
                try:
                    data = self._d__data[typ]
                except KeyError:
                    data_ex = False
                if data_ex:
                    key = typ
                    data_dict[key] = data
                
        elif self.storage_resolution == "signals":
            for typ in _signals:
                data_ex = True
                try:
                    data = self._d__data[typ]
                except KeyError:
                    data_ex = False
                if data_ex:
                    key = typ
                    data_dict[key] = data                

        elif self.storage_resolution == "processes":
            for typ in _processes:
                data_ex = True
                try:
                    data = self._d__data[typ]
                except KeyError:
                    data_ex = False
                if data_ex:
                    key = typ
                    data_dict[key] = data             
                
        elif self.storage_resolution == "off":
            data_ex = True
            try:
                data = self._d__data[_total]
            except:
                data_ex = False
            if data_ex:
                key = _total
                data_dict[key] = data
            
        else:
            raise Exception("Unknown storage resolution: "
                            +self.storage_resolution)

        return data_dict


class TwoDResponse(TwoDSpectrumBase, Saveable):
    """This class represents a single 2D spectrum
    
    Methods
    -------
    
    plot(fig=None, window=None, stype=_total, spart="real",
         vmax=None, vmin_ratio=0.5, 
         colorbar=True, colorbar_loc="right",
         cmap=None, Npos_contours=10,
         show_states=None,
         text_loc=[0.05,0.9], fontsize="20", label=None)
    
    
    """
    
    def __init__(self, keep_pathways=False, keep_stypes=True):
        self.keep_pathways = keep_pathways
        self.keep_stypes = keep_stypes
        self.t2 = -1.0
        super().__init__()


    def set_t2(self, t2):
        """Sets the t2 (waiting time) of the spectrum
        
        
        """
        self.t2 = t2


    def get_t2(self):
        """Returns the t2 (waiting time) of the spectrum
        
        """
        return self.t2
    
   
    def get_value_at(self, x, y):
        """Returns value of the spectrum at a given coordinate
        
        """
        
        legacy = False

        if legacy:
            if self.dtype is None:
                raise Exception("Data type not set")
        else:
            if self.current_dtype is None:
                raise Exception("No data type specified")

        (ix, dist) = self.xaxis.locate(x)
        (iy, dist) = self.yaxis.locate(y)    

        if legacy:
            
            if self.dtype == "Tot":
                return self.data[iy,ix]
                #return numpy.real(self.reph2D[ix,iy]+self.nonr2D[ix,iy])
            elif self.dtype == "Reph":
                return self.reph2D[iy,ix]
            elif self.dtype == "Nonr":
                return self.nonr2D[iy,ix]
        
        else:
            
            
            # FIXME: try linear interpolation
            return self.d__data[iy, ix]


    def get_area_integral(self, area, dpart=part_REAL):
        """Returns an integral of a given area in the 2D spectrum
        
        """
        def integral_square(x1, x2, y1, y2, data, dx, dy):
            (n1, n2) = data.shape
            data.reshape(n1*n2)
            return numpy.sum(data)*dy*dy
        
        area_shape = area[0]
        x1 = area[1][0]
        x2 = area[1][1]
        y1 = area[1][2]
        y2 = area[1][3]
        
        dx = self.xaxis.step
        dy = self.yaxis.step
        
        (nx1, derr) = self.xaxis.locate(x1)
        (nx2, derr) = self.xaxis.locate(x2)
        (ny1, derr) = self.yaxis.locate(y1)
        (ny2, derr) = self.yaxis.locate(y2)
        
        if area_shape == "square":
            int_fce = integral_square
        else:
            raise Exception("Unknown area type: "+area_shape)   
            
        data = self.data[nx1:nx2, ny1:ny2]
        
        if dpart == part_REAL:
            return int_fce(x1, x2, y1, y2, numpy.real(data), dx, dy)
        elif dpart == part_IMAGINARY:
            return int_fce(x1, x2, y1, y2, numpy.imag(data), dx, dy)
        elif dpart == part_ABS:
            return int_fce(x1, x2, y1, y2, numpy.abs(data))
        else:
            raise Exception("Unknown data part")


    def get_cut_along_x(self, y0):
        """Returns a DFunction with the cut of the spectrum along the x axis
        
        """
        (iy, dist) = self.yaxis.locate(y0)
        
        ax = self.xaxis
        vals = numpy.zeros(ax.length, dtype=self.d__data.dtype)
        for ii in range(ax.length):
            vals[ii] = self.d__data[iy, ii]
    
        return DFunction(ax, vals)
    

    def get_cut_along_y(self, x0):
        """Returns a DFunction with the cut of the spectrum along the y axis
        
        """
        (ix, dist) = self.xaxis.locate(x0)
        
        ay = self.yaxis
        vals = numpy.zeros(ay.length, dtype=self.d__data.dtype)
        for ii in range(ay.length):
            vals[ii] = self.d__data[ii, ix]
    
        return DFunction(ay, vals)   
    

    def get_cut_along_line(self, point1, point2, which_step=None, step=None):
        """Returns a cut along a line specified by two points
        
        """
        vx1 = point1[0]
        vy1 = point1[1]
        vx2 = point2[0]
        vy2 = point2[0]
        
        (x1, dist) = self.xaxis.locate(vx1)
        (y1, dist) = self.yaxis.locate(vy1)
        (x2, dist) = self.xaxis.locate(vx2)
        (y2, dist) = self.yaxis.locate(vy2)
        
        length = numpy.sqrt((vx1-vx2)**2 + (vy1-vy2)**2)
        
        if which_step is None:
            wstep = "x"
        else:
            wstep = which_step
            
        if step is None:
            if wstep == "x":
                dx = self.xaxis.step
            elif wstep == "y":
                dx = self.yaxis.step
            else:
                raise Exception("Step along the cut line is not defined")
                
        else:
            dx = step
             
        Nstep = int(length/dx)
        
        axis = ValueAxis(0.0, Nstep+1, dx)
            
        vals = numpy.zeros(Nstep+1, dtype=self.d__data.dtype)
        ii = 0
        for val in axis.data:
            vx1 = self.xaxis.data[x1]
            vx2 = self.xaxis.data[x2]
            vy1 = self.yaxis.data[y1]
            vy2 = self.yaxis.data[y2]
            x = vx1 + val*(vx2-vx1)/(Nstep*dx)
            y = vy1 + val*(vy2-vy1)/(Nstep*dx)
            vals[ii] = self.get_value_at(x,y)
            ii += 1
            
        return DFunction(axis, vals)    
        
    
    def get_diagonal_cut(self):
        """Returns cut of the spectrum along the diagonal 
        
        
        """
        
        point1 = [self.xaxis.min, self.yaxis.min]
        point2 = [self.xaxis.max, self.yaxis.max]
        
        fce = self.get_cut_along_line(point1, point2, which_step="x")
        
        fce.axis.data += point1[0]
        
        return fce
    

    def get_anti_diagonal_cut(self, point):
        
        pass


    def get_max_value(self):
        """Maximum value of the real part of the spectrum
        
        
        """
        legacy = False
        
        if legacy:
            return numpy.amax(numpy.real(self.reph2D+self.nonr2D))
        
        else:
            return numpy.amax(numpy.real(self.d__data))


    def get_min_value(self):
        """Minimum value of the real part of the spectrum
        
        
        """
        legacy = False
        
        if legacy:
            return numpy.amin(numpy.real(self.reph2D+self.nonr2D))
        
        else:
            return numpy.amin(numpy.real(self.d__data))
            
    
    #FIXME: introduce new storage scheme
    def devide_by(self, val):
        """Devides the total spectrum by a value
        
        
        Parameters
        ----------
        
        val : float, int
            Value by which we devide the spectrum
        
        """
        
        legacy = False
        
        if legacy:
            self.reph2D = self.reph2D/val
            self.nonr2D = self.nonr2D/val
        else:
            data_dict = self.get_all_data()
            for data_key in data_dict:
                data = data_dict[data_key]
                data[:,:] = data/val



    def get_PumpProbeSpectrum(self):
        """Returns a PumpProbeSpectrum corresponding to the 2D spectrum
        
        """
        #from .pumpprobe import PumpProbeSpectrumCalculator
        from . import pumpprobe as pp
        #from ..core.time import TimeAxis
        #fake_t = TimeAxis(0,1,1.0)
        #ppc = PumpProbeSpectrumCalculator(fake_t, fake_t, fake_t)
        #return ppc.calculate_from_2D(self)
        return pp.calculate_from_2D(self)
    
    
    
    def get_TwoDSpectrum(self, dtype=None):
        """Returns a 2D spectrum based on this response
        
        """
        if dtype is None:
            dtype = signal_TOTL
        twod = TwoDSpectrum()
        twod.set_axis_1(self.xaxis.copy())
        twod.set_axis_3(self.yaxis.copy())
        
        twod.set_t2(self.t2)
        
        twod.set_data_type(dtype)
        self.set_data_flag(dtype)
        twod.set_data(self.d__data[:,:])

        return twod


    def plot(self, fig=None, window=None, 
             stype=_total, spart=part_REAL,
             vmax=None, vmin_ratio=0.5, 
             colorbar=True, colorbar_loc="right",
             cmap=None, Npos_contours=10,
             show_states=None,
             text_loc=[0.05,0.9], fontsize="20", label=None,
             xlabel=None,
             ylabel=None,
             axis_label_font=None,
             show=False, savefig=None):
        """Plots the 2D spectrum
        
        Parameters
        ----------
        
        fig : matplotlib.figure
            Figure into which plotting will be done. This is used e.g. when
            making a movie using moview writter (may be obsolete).
            If fig is None, we create a new figure object.
            
        window : list
            Specifies the plotted window in current energy units. When axes
            are x and y, the window is specified as window=[x_min,x_max,y_min,y_max]
            
        stype : {qr.signal_TOTL, "rephasing", "non-rephasing"}
            type of the spectrum 
            
        spart : {"real", "imaginary", "abs"}
            part of the spectrum to be plotted
            
        vmax : float
            max of the plotting range in the z-direction. If vmax is None,
            maximum of the real part of the spectrum is used to determine
            the values of `vmax`
            
            
            
            
        """
        legacy = False
        
        # 
        # What type of spectra to plot
        #
        if stype == signal_TOTL:
            
            if not legacy:
                self.set_data_flag(signal_TOTL)
                spect2D = self.d__data
                
            else:
                if (self.reph2D is not None) and (self.nonr2D is not None):
                    spect2D = self.reph2D + self.nonr2D 
                elif self.reph2D is not None:
                    spect2D = self.reph2D 
                elif self.nonr2D is not None:
                    spect2D = self.nonr2D
                
        elif stype == signal_REPH:
            if not legacy:
                self.set_data_flag(signal_REPH)
                spect2D = self.d__data
            else:
                spect2D = self.reph2D
                
        elif stype == signal_NONR:
            if not legacy:
                self.set_data_flag(signal_NONR)
                spect2D = self.d__data
            else:
                spect2D = self.nonr2D     
                
        else:
            raise Exception("Undefined spectrum type "+stype)
        
        
        #
        # What part of the spectrum to plot
        #
        if spart == part_REAL:
            spect2D = numpy.real(spect2D)
        elif spart == part_IMAGINARY:
            spect2D = numpy.imag(spect2D)
        elif spart == part_ABS:
            spect2D = numpy.abs(spect2D)
        else:
            raise Exception("Undefined part of the spectrum: "+spart)
         
            
        if window is not None: 
            axis = window
            w1_min = axis[0]
            w1_max = axis[1]
            w3_min = axis[2]
            w3_max = axis[3]

            (i1_min, dist) = self.xaxis.locate(w1_min)
            (i1_max, dist) = self.xaxis.locate(w1_max)

            (i3_min, dist) = self.yaxis.locate(w3_min)
            (i3_max, dist) = self.yaxis.locate(w3_max)   
            
        else:
            i1_min = 0
            i1_max = self.xaxis.length
            i3_min = 0
            i3_max = self.yaxis.length
            
    
        #
        # Plotting with given units on axes
        #
  
        realout = spect2D[i3_min:i3_max,i1_min:i1_max]
    
        #
        #  How to treat the figures
        #
        if fig is None:
            fig, ax = plt.subplots(1,1)
        else:
            fig.clear()
            fig.add_subplot(1,1,1)
            ax = fig.axes[0]
            
        #
        # Color map
        #
        if cmap is None:
            cmap = plt.cm.rainbow
            
            
        #
        # Actual plotting
        #
        if vmax is None:
            vmax = numpy.amax(realout)

        vmin = numpy.amin(realout)
        if vmin < -vmax*vmin_ratio:
            vmax = -vmin
        else:
            vmin = -vmax*vmin_ratio
        
        Npos = Npos_contours
        poslevels = [i*vmax/Npos for i in range(1, Npos)]
        neglevels = [-i*vmax/Npos for i in range(Npos,1,-1)]
        
        levo = self.xaxis.data[i1_min]
        prvo = self.xaxis.data[i1_max-1]
        dole = self.yaxis.data[i3_min]
        hore = self.yaxis.data[i3_max-1]
        
        cm = plt.imshow(realout, extent=[self.xaxis.data[i1_min],
                                    self.xaxis.data[i1_max-1],
                                    self.yaxis.data[i3_min],
                                    self.yaxis.data[i3_max-1]],
                   origin='lower', vmax=vmax, vmin=vmin,
                   interpolation='bilinear', cmap=cmap)  

        #
        # Label
        #
        pos = text_loc
        if label is not None:
            label = label    
            ax.text((prvo-levo)*pos[0]+levo,
                (hore-dole)*pos[1]+dole,
                label,
                fontsize=str(fontsize))

        #
        # axis labels
        #
        if axis_label_font is not None:
            font = axis_label_font
        else:
            font={'size':20}

        if xlabel is None:
            xl = ""
        if ylabel is None:
            yl = ""
            
        if xlabel is not None:
            xl = r'$\omega$ [fs$^{-1}$]'

        if isinstance(self.xaxis, FrequencyAxis):
            units = self.xaxis.unit_repr_latex()
            xl = r'$\omega_{1}$ ['+units+']'
            yl = r'$\omega_{3}$ ['+units+']'
#        if isinstance(self.axis, TimeAxis):
#            xl = r'$t$ [fs]'
#            yl = r'$f(t)$'

        if xlabel is not None:
            xl = xlabel
        if ylabel is not None:
            yl = ylabel

        if xl is not None:
            plt.xlabel(xl, **font)
        if xl is not None:
            plt.ylabel(yl, **font)  


        
        #
        # Contours
        #
        
        # positive contours are always plotted
        try:
            plt.contour(self.xaxis.data[i1_min:i1_max],
                     self.yaxis.data[i3_min:i3_max],
                     realout, levels=poslevels, colors="k")
                     #linewidth=1)
        except:
            print("No positive contours found; not plotted")
              
        # other contours only if we do not plot absolute values
        if spart != "abs":
            # zero contour
            try:
                plt.contour(self.xaxis.data[i1_min:i1_max],
                         self.yaxis.data[i3_min:i3_max],
                         realout, levels=[0],colors="b")
                         #linewidth=1)
            except:
                print("Zero contour not found; not plotting")
        
        
            # negative contours
            try:
                plt.contour(self.xaxis.data[i1_min:i1_max],
                         self.yaxis.data[i3_min:i3_max],
                         realout, levels=neglevels,colors="k")
                         #linewidth=1) 
            except:
                print("Negative contour not found; not plotting")
        
        #
        # Color bar presence
        #
        if colorbar:
            plt.clim(vmin=vmin,vmax=vmax)
            fig.colorbar(cm)
            
        #
        # Plot lines denoting positions of selected transitions
        #
        if show_states is not None:
            for en in show_states:  
                plt.plot([en,en],[dole,hore],'--k',linewidth=1.0)
                plt.plot([levo,prvo],[en,en],'--k',linewidth=1.0)
            
        #
        # Should the spectra be showed now?
        #
        if show:
            self.show()
            
        #
        # Saving the figure
        #
        if savefig:
            self.savefig(savefig)

            
    def show(self):
        """Show the plot of 2D spectrum
        
        By default, plots are not shown. It is waited until explicit show()
        is called
        
        """
        
        plt.show()
        

    def savefig(self, filename):
        """Saves the fige of the plot into a file
        
        """
        plt.savefig(filename)
            

    def trim_to(self, window=None):
        """Trims the 2D spectrum to a specified region
        
        Parameters
        ----------
        
        window : array
            Spectral window to which the present 2D spectrum 
            should be trimmed. The window is specified as
            an array [w1_min, w1_max, w3_min, w3_max], where 
            w1_min is the lower bound of the w1 axis (the x-axis),
            w1_max is the upper bound of the w1 axis and similarly
            for the w3 axis (y-axis) of the spectrum.
            
        
        """
        
        if window is not None:

            axis = window
            w1_min = axis[0]
            w1_max = axis[1]
            w3_min = axis[2]
            w3_max = axis[3]

            (i1_min, dist) = self.xaxis.locate(w1_min)
            (i1_max, dist) = self.xaxis.locate(w1_max)

            (i3_min, dist) = self.yaxis.locate(w3_min)
            (i3_max, dist) = self.yaxis.locate(w3_max)    
            
            # create minimal off-set
            i1_min -=1
            i1_max +=1
            i3_min -=1
            i3_max +=1
            
            # reconstruct xaxis
            start_1 = self.xaxis.data[i1_min]
            length_1 = i1_max - i1_min
            step_1 = self.xaxis.step
            atype = self.xaxis.atype
                
            xaxis = FrequencyAxis(start_1,length_1,step_1, 
                                  atype=atype,
                                  time_start=self.xaxis.time_start)
            self.xaxis = xaxis
            
            # reconstruct yaxis
            start_3 = self.yaxis.data[i3_min]
            length_3 = i3_max - i3_min
            step_3 = self.yaxis.step
            yaxis = FrequencyAxis(start_3,length_3,step_3, 
                                  atype=self.yaxis.atype,
                                  time_start=self.yaxis.time_start)                
            self.yaxis = yaxis            
                
            dtype_saved = self.current_dtype
            
            if self.storage_resolution == "pathways":
                for typ in _ptypes:
                    piece_ex = True
                    try:
                        piece = self._d__data[typ]
                    except:
                        piece_ex = False
                    if piece_ex:
                        for tag in piece.keys():
                            data_ex = True
                            try:
                                data = piece[tag]
                            except KeyError:
                                data_ex = False
                            if data_ex:
                                ndata = data[i1_min:i1_max,i3_min:i3_max]
                                piece[tag] = ndata                       
                    
            elif self.storage_resolution == "types":
                
                for typ in _ptypes:
                    self.set_data_flag(typ)
                    #data_ex = True
                    #try:
                    data = self.d__data
                    #except KeyError:
                    #    data_ex = False
                    if data is not None:
                        ndata = data[i1_min:i1_max,i3_min:i3_max]
                        self.d__data = ndata
                    
            elif self.storage_resolution == "signals":
                
                for typ in _signals:
                    self.set_data_flag(typ)
                    #data_ex = True
                    #try:
                    data = self.d__data
                    #except KeyError:
                    #    data_ex = False
                        
                    if data is not None:
                        ndata = data[i1_min:i1_max,i3_min:i3_max]
                        self.d__data = ndata
                    
            elif self.storage_resolution == "processes":
                
                for typ in _processes:
                    self.set_data_flag(typ)
                    #data_ex = True
                    #try:
                    data = self.d__data
                    #except KeyError:
                    #    data_ex = False
                    if data is not None:
                        ndata = data[i1_min:i1_max,i3_min:i3_max]
                        self.d__data = ndata
                        
            elif self.storage_resolution == _total:
                self.set_data_flag(_total)
                #data_ex = True
                #try:
                data = self.d__data
                #except KeyError:
                #    data_ex = False
                if data is not None:
                    ndata = data[i1_min:i1_max,i3_min:i3_max]
                    self.d__data = ndata
                
            self.set_data_flag(dtype_saved)
                
        else:
            # some automatic trimming in the future
            pass
