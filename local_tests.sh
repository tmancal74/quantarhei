#! /bin/sh

#
# This script runs a test of the local source code of Quantarhei package.
# To test an installed version of Quantarhei run
#
#   > paver
#
# 
#
#
#


#
# Make sure now quantarhei is installed
#
pip uninstall quantarhei
#pip uninstall aceto

# point to the local version of the package
export PYTHONPATH=`pwd`

#
# Run tests 
#
#nosetests -vs tests/unit/core/test_valueaxis.py
#nosetests -vs tests/unit/core/test_saveable.py
#nosetests -vs tests/unit/core/time_test.py
#nosetests -vs tests/unit/core/frequency_test.py
#nosetests -vs tests/unit/core/test_dfunction.py
#nosetests -vs tests/unit/qm/corfunctions/correlationfunctions_test.py
#nosetests -vs tests/unit/builders/test_molecules.py
#nosetests -vs tests/unit/builders/test_aggregates.py
nosetests -vs tests/unit/spectroscopy/abs_test.py

#python examples/demo_load.py


paver


