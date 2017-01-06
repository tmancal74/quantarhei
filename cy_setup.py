# -*- coding: utf-8 -*-

from setuptools import find_packages
from distutils.core import setup
from Cython.Build import cythonize


packages = find_packages(exclude=['quantarhei.implementations.*',
                                  'quantarhei.implementations',
                                  'tests.*','tests' ,'docs'])

packg = list()
for pckg in packages:
    p = pckg+"/*.py"
    packg.append(p)
    
for p in packg:
    print(p)
    
setup(ext_modules=cythonize(packg))
