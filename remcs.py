# -*- coding: utf-8 -*-

""" Removes *.c files created during cythonization

"""

from setuptools import find_packages

import os
import glob

def process_packages(packages):
    packg = list()
    for pckg in packages:
        p = pckg.replace(".", "/")+"/*.c"
        packg.append(p)
    
    for p in packg:

        fls = glob.glob(p)
        for fl in fls:
            print(fl)
            os.remove(fl) 
            
    return packg


if __name__ == "__main__":
    process_packages(find_packages(exclude=['quantarhei.implementations.*',
                                   'quantarhei.implementations',
                                   'tests.*','tests' ,'docs']))