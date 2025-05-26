# -*- coding: utf-8 -*-

from .triangle import triangle

class unique_list:
    """Helper list class that stores only unique objects
    
    Every time an object is stored, it is checked if it is already
    present. Only the index of the object is stored. One can ask
    for the number of such unique obejcts and for their list.
    
    Parameters
    ----------

    N : int
        Number of elements of the list. Default is 0.
        
        
    Examples
    --------
    
    >>> A = unique_list()
    >>> A.get_number_of_unique_elements()
    0
    
    Prealocating the size of the list does not affect the number of elements
    >>> A = unique_list(6)
    >>> A.get_number_of_unique_elements()
    0
    
    But one can start setting the elements
    
    >>> a = "x"
    >>> b = "y"
    >>> A.set_element(4, a)
    >>> A.set_element(3, b)
    >>> A.get_number_of_unique_elements()
    2
    
    When a possition is accessed to which the element was not set, 
    Exception occurs.
    
    >>> c = A.get_element(2)
    Traceback (most recent call last):
        ...
    Exception: Element not set
    
    
    Setting an element in an empty list raises an Exception
    
    >>> A = unique_list()
    >>> A.set_element(0,"z")
    Traceback (most recent call last):
        ...
    IndexError: list index out of range
   
    """
    
    def __init__(self,N=0):
        if N == 0:
            self._storage = []
            self._indices = []
        else:
            self._storage = []
            self._indices = [-1]*N
        
        
    def append(self,obj):
        """Appends a new object to the list

        Examples
        --------
    
        >>> A = unique_list()
        >>> A.append("x")
        >>> print(A.get_number_of_unique_elements())
        1
        
        """
        try:
            ind = self._storage.index(obj)
            self._indices.append(ind)
        except:
            self._storage.append(obj)
            ind = self._storage.index(obj)
            self._indices.append(ind)
            
    def set_element(self,i,obj):
        """Sets an object to the specified position 
        
        """
        # FIXME: when replacing an element, the list may get corrupted
        
        # old index of the object which sits under i
        oldind = self._indices[i]
        try:
            ind = self._storage.index(obj)
            self._indices[i] = ind
        except:
            self._storage.append(obj)
            ind = self._storage.index(obj)
            self._indices[i]= ind
            
        
    def get_element(self,i):
        """Returns an element residing on a given position
        
        """
        k = self._indices[i]
        if k == -1:
            raise Exception("Element not set")
        else:
            return self._storage[k]


    def get_number_of_unique_elements(self):
        """Returns a number of unique elements stored 
        
        Examples
        --------

        >>> ul = unique_list()
        >>> a = "x"
        >>> b = "y"
        >>> ul.append(a)
        >>> ul.append(b)
        >>> ul.append(b)
        >>> ul.append(a)
        >>> ul.get_number_of_unique_elements()
        2

        """
        return len(self._storage)

    #@deprecated        
    def get_unique_element(self):
        return self._storage
        
    def get_unique_elements(self):
        """Returns a list of elements stored 
        
        Examples
        --------

        >>> ul = unique_list()
        >>> a = "x"
        >>> b = "y"
        >>> ul.append(a)
        >>> ul.append(b)
        >>> ul.append(b)
        >>> ul.append(a)
        >>> lst = ul.get_unique_elements()
        >>> print(lst)
        ['x', 'y']

        """
        return self._storage
        
        
class unique_triangle_array:
    
    def __init__(self,N):
        self.N = N
        self.triangle = triangle(self.N)
        self._storgage = unique_list(N=N)
        
    def set_element(self,i,j,obj):
        I = self.triangle.locate(i,j)
        self._storgage.set_element(I,obj)
        
    def get_element(self,i,j):
        I = self.triangle.locate(i,j)
        return self._storgage.get_element(I)
        
    def get_number_of_unique_elements(self):
        return len(self._storgage._storage)

    #@deprecated        
    def get_unique_element(self):
        return self._storage._storage

    def get_unique_elements(self):
        return self._storage._storage

        
class unique_array:
    
    def __init__(self,N,M):
        self.N = N
        self.M = M
        self._storage = unique_list(N=self.N*self.M)
        
        
    def get_element(self,i,j):
        if ((i<self.N) and (j < self.M)) and ((i>=0) and (j>=0)) :
            return self._storage.get_element((self.M-1)*i+j)
        else:
            raise Exception("Index out of range")
        
    def set_element(self,i,j,obj):
        if ((i<self.N) and (j < self.M)) and ((i>=0) and (j>=0)) :
            self._storage.set_element((self.M-1)*i+j,obj)
        else:
            raise Exception("Index out of range")
            
    def get_number_of_unique_elements(self):
        return len(self._storgage._storage)       

    #@deprecated        
    def get_unique_element(self):
        return self._storage._storage
    
    def get_unique_elements(self):
        return self._storage._storage
        