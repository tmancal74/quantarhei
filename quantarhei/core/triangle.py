# -*- coding: utf-8 -*-
from .saveable import Saveable

class triangle(Saveable):
    """Class representing a symmetric matrix by a linear list

    """
    def __init__(self,N=0):
        self.N = N
    
    def get_empty_list(self):
        return self.get_list(init=None)
        
    def get_list(self,init=None):
        return [init]*(((self.N**2)-self.N)//2+self.N)

    
    def locate(self,i,j,transpose=True,report_transpose=False):
        
        if (((i>=self.N) or (j>=self.N)) or ((i < 0) or (j < 0))):
            raise Exception("Index out of range")
            
        trans = False
        if (j > i):
            if (transpose):
                trans = True
                pom = i
                i = j
                j = pom
            else:
                raise Exception("Index out of range (transposed is in range)")
            
        I = 0
        for m in range(1,i+1):
            I -= (m-1)    
        I += (j-i)
        I += i*self.N
        
        if report_transpose:
            return trans,I
        else:
            return I
            
    def indices(self,I):
        pass
    
    
#FIXME: define triangle_list class which can store using triangle