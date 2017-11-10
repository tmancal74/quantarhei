# -*- coding: utf-8 -*-
"""
    This script copies the examples under tests to
    the standard example directory and processes
    certain setting variables. Under the tests, the
    examples do not plot and save data. In examples
    directory they do.


"""
import os
from shutil import copyfile

def get_example_files(loc):
    """Returns files containing examples
    
    """
    onlyfiles = [f for f in os.listdir(loc) 
                 if os.path.isfile(os.path.join(loc, f))]
    return onlyfiles
    
def get_prefix(file):
    """Returns the prefix of the file
    
    """
    parts = file.split(sep="_")
    return parts[0]

    
def get_number(file):
    """Returns the number part of the file name
    
    """
    parts = file.split(sep="_")
    return parts[1]


def get_root(file):
    """Returns the root part of the filename
    
    """
    parts = file.split(sep="_")
    ret = ""
    for k in range(len(parts)):
        part = parts[k]
        if k > 1:
            ret += "_"+part
    return ret

def copy_file(from_loc, file, to, as_file, fltr=None):
    """Copies the file into "examples" directory
    
    """
    if fltr is None:
        src = os.path.join(from_loc, file)
        dst = os.path.join(to, as_file)
        copyfile(src, dst)
    
def fltr_fce(text):
    """Changes predefined variables in examples
    
    """
    pass
    


exloc = os.path.join("..", "quantarhei", "wizard", "examples")

files = get_example_files(exloc)

print("Copying example files from ", exloc)
for file in files:
    pref = get_prefix(file)
    if pref == "ex":
        nr = get_number(file)
        root = get_root(file)
    
        demo_file = "demo_"+nr+root
    
        print(file,"-->", demo_file)
        copy_file(from_loc=exloc, file=file,
                  to=".", as_file=demo_file) #, fltr=fltr_fce)
    
    