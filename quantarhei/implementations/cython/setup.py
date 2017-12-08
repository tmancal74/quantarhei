# -*- coding: utf-8 -*-

import numpy

from distutils.core import setup
from distutils.extension import Extension
#from Cython.Build import cythonize
from Cython.Distutils import build_ext

ext_modules = [Extension("loopit", ["loopit.pyx"])]

setup(
    name = "My hello app",
    cmdclass = {'build_ext': build_ext},
    include_dirs = [numpy.get_include()],
    ext_modules = ext_modules, #cythonize('loopit.pyx'), 
    )