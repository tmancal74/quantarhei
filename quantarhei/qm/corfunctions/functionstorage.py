# -*- coding: utf-8 -*-

from itertools import combinations

import numpy
from ...core.time import TimeAxis


class FunctionStorage:
    """Data storage for discrete representation of correlation functions
    and lineshape functions.


    Parameters
    ----------

    N: int
        Number of functions to be stored.
    
    
    """
    
    import numpy 
    
    def __init__(self, N=1, timeaxis=None, dtype=numpy.complex64,
                 show_config=False, config=None):

        #
        # Storage configuration
        #
        # This dictionary describes what g(t) values will be stored and how
        #
        if config is None:
            self.config = {"response_times":{
                            "t1":{"reset":False,"integral":{"s1"},"axis":0},
                            "t2":{"reset":True, "integral":None},
                            "t3":{"reset":False,"integral":{"s3"},"axis":1}}} #,
                            #"t4":{"reset":False,"integral":None,"axis":2}}}

        else:
            self.config = config

        #######################################################################
        #
        # Analyze the configuration 
        #
        #######################################################################
        
        #
        # generate the list of argument combinations
        #
        values = self.config["response_times"]
        Nf = len(values)  # Total number of elements
        
        # Generate all M-tuples (1 ≤ M ≤ N) while preserving order
        all_tuples = [list(combinations(values, M)) for M in range(1, Nf + 1)]

        # Flatten and print results
        self.time_combination_tuples = []
        self.time_combination_strings = []
        for tuples in all_tuples:
            for t in tuples:
                self.time_combination_tuples.append(t)
                tstr = ""
                kk = 0
                for x in t:
                    if kk > 0:
                        tstr += "+"
                    tstr += x
                    kk += 1
                self.time_combination_strings.append(tstr)
                
        if show_config:
            print("\nStorage configuration")
            print("Time argument combinations: ", self.time_combination_strings)

        #
        # creating time_index
        #
        time_index = {}

        # find all reset times 
        reset_times = []
        for rst_key in self.config["response_times"]:
            if self.config["response_times"][rst_key]["reset"]:
                reset_times.append(rst_key)
                time_index[rst_key] = -1
                
        # find all 1D times
        oned_times = []
        for rst_key in self.config["response_times"]:
            if not self.config["response_times"][rst_key]["reset"]: 
                axis = self.config["response_times"][rst_key]["axis"]
                oned_times.append(rst_key)
                time_index[rst_key] = axis

        # is the number of 1D times matching the number of submitted time axes?
        if isinstance(timeaxis, TimeAxis):
            tlen = 1
        else:
            tlen = len(timeaxis)
        if len(oned_times) != tlen:
            print("Storage configuration:")
            print(self.config)
            raise Exception("The number of time axes has"
                           +" to be consistent with storage conguration.")

        self.oned_times = oned_times

        # find all combinations with 2 or more times and calculate their dimension
        for tpl, tstr in zip(self.time_combination_tuples, self.time_combination_strings):
            if len(tpl) > 1:
                dim = []
                for tm in tpl:
                    #if time_index[tm] >= 0:     # we add even -1 to time_index
                    dim.append(time_index[tm])
                if len(dim) == 1:
                    time_index[tstr] = dim[0]
                elif len(dim) > 1:
                    time_index[tstr] = tuple(dim)
                else:
                    raise Exception("Inconsistent configuration")

        # order time_index        
        # Sorting criteria:
        # - First, sort by type: Integers should come before tuples
        # - Second, for integers: Sort numerically
        # - Third, for tuples: Sort by length
        sorted_items = sorted(time_index.items(), 
                    key=lambda x: (isinstance(x[1], tuple), x[1] if isinstance(x[1], int) else len(x[1])))
        
        # Convert back to a dictionary (optional)
        time_index = dict(sorted_items)  

        if show_config:
            print("Property 'time_index': ", time_index) 
            #print(self.time_combination_tuples)
            
        ##########################################################################
        #
        ##########################################################################
        
        # Number of unique functions stored
        self.N = N 

        # The max value of the index identifying the stored functions for the user
        self.Nmax = self.N - 1

        # One or more time axes on which the functions should be represented
        if timeaxis is None:
            raise Exception("At least one time axis has to be specified")
        self.timeaxis = timeaxis
        
        # Storage dimensions; they correspond to the submitted time axes
        if isinstance(timeaxis, TimeAxis):
            self.Ndim = 1
            dim = numpy.zeros(self.Ndim, dtype=numpy.int32)
            dim[0] = timeaxis.length
        else:
            self.Ndim = len(timeaxis)
            dim = numpy.zeros(self.Ndim, dtype=numpy.int32)
            for ii, ta in enumerate(timeaxis):
                if not isinstance(ta, TimeAxis):
                    raise Exception("'timeaxis' argument must be a list/tuple of Quantarhei TimeAxis objects.")
                dim[ii] = ta.length
        self.dim = dim

        # Data type of the storage (default is complex64)
        self.dtype = dtype

        # Dictionary of time arguments and indices of the 'dim' argument defining the corresponding
        # length of the discrete representation. The mapping of the time arguments on the storage integer
        # indices is done in a natural order.
        #
        # The default setting means that the value stored under the index 0 is the one of g(t2) and 
        # its representation has a size of 1 element (this is the meaning of -1). The value stored
        # under index 6 is the one of g(t1+t2+t3) with a length of dim[0]*dim[1] (this what the values 
        # in the tuple define).
        #self.time_index = {"t2":-1,"t1":0,"t1+t2":0, "t3":1,"t2+t3":1,"t1+t3":(0,1),"t1+t2+t3":(0,1)}
        self.time_index = time_index
        
        reshapes = dict()
        for dms in self.time_index:
            val = self.time_index[dms]
            if isinstance(val, int):
                if val == -1:
                    reshapes[dms] = ""
                else:
                    reshapes[dms] = [self.dim[val]]

            else:
                rlist = []
                dcount = 0
                kl = 0
                for kk, dims in enumerate(val):

                    if dims >= 0:

                        rlist.append(self.dim[dims])
                        dcount += 1
                        kl += 1
                
                reshapes[dms] = rlist

        self.reshapes_dic = reshapes
        self.reshapes = []
        for key in self.reshapes_dic:
            self.reshapes.append(self.reshapes_dic[key])
            

        # The number of time arguments
        self.Nt = len(self.time_index)

        # mapping of times to indices
        self.time_mapping = {}
        kk = 0
        for tm in self.time_index:
            self.time_mapping[tm] = kk
            kk += 1

        #
        # Extract the sizes from 'dim' argument
        #
        # self._N stores the sizes of individual index representations
        self._N = numpy.zeros(self.Nt, dtype=numpy.int32)

        # loop over time indices
        kk = 0
        for ts in self.time_index:
            ival = self.time_index[ts]
            if isinstance(ival, tuple):
                self._N[kk] = 1
                for sz in ival:
                    if sz >= 0:
                        self._N[kk] *= dim[sz]  
            elif isinstance(ival,int):
                if ival == -1:
                    self._N[kk] = 1
                else:
                    self._N[kk] = dim[ival]
            kk += 1
                
        # one function takes one dimensional array of length 'data_stride'                
        self.data_stride = 0
        for kk in range(self.Nt):
            self.data_stride += self._N[kk]  
        
        # total number of values 
        self.data_dim = self.N*self.data_stride

        # data size in MB
        #print(self.data_dim)
        #print(numpy.dtype(self.dtype))
        self.data_size = \
            self.data_dim*numpy.dtype(self.dtype).itemsize/(1024*1024)
        self.size_units = "MB"

        #
        # Here the data are stored
        #
        self.data = numpy.zeros(self.data_dim, dtype=self.dtype)
        
        # 2D view of the data for fast access
        self._data2d = self.data.reshape((self.N, self.data_stride))

        #
        # Here the g(t) functions are stored
        #
        self.funcs = {}

        # The number of stored functions
        self.Nf = 0

        #
        # These arrays store the starts and ends of the g(t) arrays in a give stride
        #
        self.start = numpy.zeros(self.Nt, dtype=numpy.int32)
        self.end = numpy.zeros(self.Nt, dtype=numpy.int32)

        # filling the starts
        for kk in range(self.Nt):
            if kk == 0:
                self.start[kk] = 0
            else:
                self.start[kk] = self.start[kk-1] + self._N[kk-1]

        # filling the ends
        for ii in range(self.Nt):
            self.end[ii] = self.start[ii] + self._N[ii]

        #
        # As an alternative, one can use the string to access the data
        #
        self.start_dic = {}
        self.end_dic = {}
        kk = 0
        for tm in self.time_mapping:
            self.start_dic[tm] = self.start[kk]
            self.end_dic[tm] = self.end[kk]
            kk += 1

        #
        # Mapping between the actual stored functions and the index representing
        # the functions for the outside world.
        #
        self.mapping = numpy.zeros(self.Nmax+1,dtype=numpy.int32)
        for kk in range(self.N):
            self.mapping[kk] = -1

        if show_config:
            print("\n")


    def __str__(self):
        """String representation of the storage

        """
        return ("Correlation function storage\n"
               +" dimension: "+str(self.dim)+"\n"
               +" data size: "+str(self.data_size)+" "+self.size_units)


    def __setitem__(self, index, value):
        """Set the value(s) of the storage

        """
        if isinstance(index, tuple):
            
            if len(index) == 2:
                i, j = index
                start = i*self.data_stride
                
                if isinstance(j,int):
                    sta = start + self.start[j]
                    end = start + self.end[j]
                    
                else:
                    sta = start + self.start_dic[j]
                    end = start + self.end_dic[j]
                
                self.data[sta:end] = numpy.array(value, dtype=self.dtype) 
                
            else:
                raise Exception()
                
        else:
            raise Exception()
            

    def __getitem__(self, index):
        """Returns the stored function 


        """
        
        if isinstance(index, tuple):
            
            if len(index) == 2:
                
                i, j = index
                
                # if i is an integer, we return one array
                if isinstance(i, int):
                
                    start = self.mapping[i]*self.data_stride
                    
                    if isinstance(j,int):
                        sta = start + self.start[j]
                        end = start + self.end[j]
                        
                    else:
                        sta = start + self.start_dic[j]
                        end = start + self.end_dic[j]
                    
                    if sta + 1 == end:
                        return self.data[sta]
                    else:
                        if isinstance(j,int):
                            tpl = self.reshapes[j]
                        else:
                            tpl = self.reshapes_dic[j]
                            
                        return self.data[sta:end].reshape(tpl)
                
                # if i is a slice :, we return a view an all arrays
                elif isinstance(i, slice) and i == slice(None):
                    
                    if isinstance(j,int):
                        sta = self.start[j]
                        end = self.end[j]
                        
                    else:
                        sta = self.start_dic[j]
                        end = self.end_dic[j]
                        
                    if sta + 1 == end:
                        return self._data2d[i, sta]
                    else:
                        if isinstance(j,int):
                            tpl = [self._data2d.shape[0]]+self.reshapes[j]
                        else:
                            tpl = [self._data2d.shape[0]]+self.reshapes_dic[j]
                            
                        return self._data2d[i, sta:end].reshape(tpl)                   
            
            else:
                raise Exception()
                
        else:
            raise Exception()                    


    def set_goft(self, N, func=None):
        """ Sets the values for a given stored function. 


        """
        #
        # The argument N must be an integer or a list of indices by which the function should
        # represented for the outside world
        #
        if not isinstance(N, (int, list, tuple, numpy.ndarray)):
            raise Exception("Argument N has to be an integer or an array (list, tuple) of integers")

        #
        # Check the consistency of the submitted time axes
        #
        self.ta = [] # This is a synonym for the self.timeaxis
        if isinstance(self.timeaxis, (tuple, list)):
            self.ta = self.timeaxis
            #self.t1a = self.timeaxis[0]
            #self.t3a = self.timeaxis[1]
        else:
            self.ta = [self.timeaxis]
            #self.t1a = self.timeaxis
            #self.t3a = self.timeaxis
        
        # for an integer
        if isinstance(N, int):
            self._check_and_make_space(N)
            
            # FIXME: mapping has to be merged with the previous mappings
            if self.mapping[N] > 0:
                raise Exception("Function already set for "+str(N))
            # set the function to the first available position
            if self.Nf > self.N:
                raise Exception("The storage full. Unable to add additional functions")
                
            if func is not None:
                pos = self._fce_already_stored(func)
                if pos == -1:
                    self.funcs[self.Nf] = func
                    self.mapping[N] = self.Nf
                    self.Nf += 1
                else:
                    self.mapping[N] = pos
                    
            else:
                raise Exception("Function not submitted")
                
        else:
            # In this case N is a list of mappings to the submitted function
            added = False
            kk = 0
            for nn in N:
                self._check_and_make_space(nn)
                    
                # if the position is free
                if self.mapping[nn] < 0:

                    # check that the same function was not submitted before
                    pos = self._fce_already_stored(func)
                    if pos == -1:
                        if kk == 0:
                            self.funcs[self.Nf] = func
                            added = True
                        self.mapping[nn] = self.Nf
                        kk += 1
                    else:
                        self.mapping[nn] = pos
                        
                else:
                    raise Exception("Function already set for "+str(nn))

            # we added just one function
            if added:
                self.Nf += 1


    def _check_and_make_space(self, nn):
        """Check if the required position is outside the allocated mapping
           and if so, make more space.

        """
        if nn > self.Nmax:                
            # create more space
            mapping = numpy.zeros(nn+1,dtype=numpy.int32)
            mapping[:] = -1
            # copy earlier values
            #print(nn, self.Nmax, mapping.shape, self.mapping.shape)
            mapping[:self.Nmax+1] = self.mapping
            # set new Nmax
            self.Nmax = nn
            # make larger variables the object properties
            self.mapping = mapping   

    
    def _fce_already_stored(self, func):
        """If the function was already stored in the object, the function returns its position,
           otherwise it returns -1.

        """
        ret = -1
        if len(self.funcs) == 0:
            return ret
            
        for key in self.funcs:
            fc = self.funcs[key]
            if fc is func:
                ret = key
        
        return ret
        

    # FIXME: allow multiple reset times
    def create_data(self, reset=dict(t2=0.0)):
        """We create data for all submitted functions


        """
        tt_matrix = []
        t2 = reset["t2"]

        t2m = numpy.array([t2])
        _colon_ = slice(None, None, None)

        for tm in self.time_index:

            # list of time axes from the time_index
            dms = self.time_index[tm]
            
            # handling single variables
            if isinstance(dms, int):
                if dms == -1:
                    tt_matrix.append(t2m)
                else:
                    tt_matrix.append(self.ta[dms].data)  

            # handling tuples
            elif isinstance(dms, (tuple, list)):

                # from tuples we need to eliminate -1 
                # to get the final dimension of the stored representation
                tt_dim = []
                for elem in self.time_index[tm]:
                    if elem >= 0:
                        tt_dim.append(self.dim[elem])

                # tt_dim is the final dimension and we allocate space for time variable 

                tt = numpy.zeros(tt_dim, dtype=numpy.float32)

                #  now let's go axis by axis
                dim_needed = 0
                for ldm in dms:
                    # if we find -1, we add the corresponding single value (now it is t2)
                    # We do it like this always, no matter how big is 'tt' (constant can be broadcasted 
                    # automatically)
                    if ldm == -1:
                        tt += t2
  
                    # if it is not -1
                    else:
                        # but it is the only dimension available, we have a 1D case
                        if len(tt_dim) == 1:
                            tt += self.ta[ldm].data

                        # otherwise we have multi-dimensional case, 
                        # and we have to combine time axes into multi-D matrix
                        # We have to be aware of which axis we are - this is given by 'ldm' variable
                        else:
                            #
                            # Here we have to figure out how to construct the index for broadcasting
                            #
                            index = [None]*len(tt_dim)
                            #index[ldm] = _colon_
                            index[dim_needed] = _colon_
                            index = tuple(index)
                            
                            tt += self.ta[ldm].data.__getitem__(index)

                            dim_needed += 1
                            
                tt_matrix.append(tt)
        
        #
        # loop over all stored functions
        #   
        for key in self.funcs:
            
            func = self.funcs[key]
            #start = key*self.data_stride

            for kk, tm in zip(range(self.Nt), self.time_index):
                self.__setitem__((key, tm),
                                 func(tt_matrix[kk]).reshape(self._N[kk]))


    def effective_size(self):
            """ Effective size of the storage. Some stored functions have the same functional form
    
            """
            return len(self.mapping)


    def get_number_of_functions(self):
        """ Returns the number of different stored functions
        
        """
        return numpy.max(self.mapping)+1


    def get_number_of_sites(self):
        """ Returns the number of assigned sites
        
        """
        return (self.mapping >= 0).sum()
        

    def get_mapping_matrix(self):
        """ Returns a matrix that maps site index on the function index
        

        M[n,i] equals 1 if the site n has the correlation function stored
        at the position i in the storage

        Returns
        -------
        Mapping matrix
        

        """

        Nstore = self.get_number_of_functions()
        Nsites = self.get_number_of_sites()

        if Nsites != len(self.mapping):
            raise Exception("Mapping has to be complete."
                        +" All of its elements have to be set without gaps.")
            
        mtrx = numpy.zeros((Nsites,Nstore), dtype=numpy.int32)
        for ii in range(Nsites):
            jj_ind = self.mapping[ii]
            mtrx[ii,jj_ind] = 1.0
        
        return mtrx

        
    def get_reorganization_energies(self):
        """Returns the estimate of the reoganization energies of the stored functions
        
        """
        index = (slice(None,None,None),"t1")
        igg = -numpy.imag(self.__getitem__(index))
        tal = self.ta[0].length - 1
        print("t1_max = ", self.ta[0].data[tal])
        lam = igg[:,tal]/self.ta[0].data[tal]
        
        return lam
            


class FastFunctionStorage(FunctionStorage):
    """ Function storage with no overhead retrieval

    """
    
    def __getitem__(self, index):
        i, j = index
        
        # if the index is integer
        if isinstance(i, int):
            start = self.mapping[i]*self.data_stride
            sta = start + self.start[j]
            end = start + self.end[j]
            if sta + 1== end:
                return self.data[sta]
            return self.data[sta:end]
        
        # i is a slice :
        elif isinstance(i, slice) and i == slice(None):
            sta = self.start[j]
            end = self.end[j]
            if sta + 1 == end:
                return self._data2d[i, sta]
            return self._data2d[i, sta:end]
        
        else:
            raise Exception("Integer or slice : required as a first index")
            




class SingleGoft(FunctionStorage):
    """ Storage of a single lineshape function 

    """
    def __init__(self, timeaxis=None, dtype=numpy.complex64,
                 show_config=False, config=None):
        super().__init__(N=1, timeaxis=timeaxis,
                         dtype=dtype, show_config=show_config, config=config)


    def __getitem__(self, index):
        """ Returns the stored function 


        """  
        j = index
        start = 0
        
        if isinstance(j,int):
            sta = start + self.start[j]
            end = start + self.end[j]
            
        else:
            sta = start + self.start_dic[j]
            end = start + self.end_dic[j]
            
        return self.data[sta:end]
