# -*- coding: utf-8 -*-
import dill
from .managers import Manager

class Parcel:
    
    def set_content(self, obj):
        """Set the content of the parcel
        """
        self.content = obj
        
        self.class_name = "class name"
        self.qrversion = Manager().version
        
        
    def save(self, filename):
        """Saves the parcel to a file
        
        """
        
        with open(filename, "wb") as f:
            dill.dump(self, f)
            
        
def load_parcel(filename):
    pass


def check_parcel(filename):
    pass
