# -*- coding: utf-8 -*-

from .triangle import triangle

class unique_list:
    """Stores only unique objects
    
    
    """
    
    def __init__(self,N=0):
        if N == 0:
            self._storage = []
            self._indices = []
        else:
            self._storage = []
            self._indices = [-1]*N
        
        
    def append(self,obj):
        try:
            ind = self._storage.index(obj)
            self._indices.append(ind)
        except:
            self._storage.append(obj)
            ind = self._storage.index(obj)
            self._indices.append(ind)
            
    def set_element(self,i,obj):
        try:
            ind = self._storage.index(obj)
            self._indices[i] = ind
        except:
            self._storage.append(obj)
            ind = self._storage.index(obj)
            self._indices[i]= ind
            
        
    def get_element(self,i):
        return self._storage[self._indices[i]]
        
    def get_number_of_unique_elements(self):
        return len(self._storage)
        
    def get_unique_element(self):
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
        
    def get_unique_element(self):
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
        
    def get_unique_element(self):
        return self._storage._storage
        