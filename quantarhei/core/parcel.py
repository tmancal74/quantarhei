# -*- coding: utf-8 -*-

import dill as pickle

    
from .managers import Manager

class Parcel:
    
    def set_content(self, obj):
        """Set the content of the parcel
        """
        self.content = obj
        
        self.class_name = "{0}.{1}".format(obj.__class__.__module__,
                                           obj.__class__.__name__)
        self.qrversion = Manager().version
        
        self.comment = ""

    
    def set_comment(self, comm):
        if comm is not None:
            self.comment = comm
            
        
    def save(self, filename):
        """Saves the parcel to a file
        
        """
        
        with open(filename, "wb") as f:
            pickle.dump(self, f)

      
def save_parcel(obj, filename, comment=None):
    """Saves a given object as a parcel 
    
    """
    p = Parcel()
    p.set_content(obj)
    p.set_comment(comment)
    
    p.save(filename)

        
def load_parcel(filename):
    """Loads the object saved as parcel
    
    """
    with open(filename, "rb") as f:
        obj = pickle.load(f)
        
    if isinstance(obj, Parcel):
        return obj.content
    else:
        raise Exception("Only Quantarhei Parcels can be loaded")
        


def check_parcel(filename):
    """Checks the content of a Quantarhei parcel
    
    """
    with open(filename, "rb") as f:
        obj = pickle.load(f)
        
    if isinstance(obj, Parcel):
        return dict(class_name=obj.class_name, qrversion=obj.qrversion,
                    comment=obj.comment)
    else:
        raise Exception("The file does not represent a Quantarhei parcel")    
