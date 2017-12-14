##############################################################
#
# Quantarhei package Makefile
#
# This file defines several tasks useful to Quantarhei
# developers 
#
##############################################################

#
# Set this to required version or override from command line
# Default is the current development version 
#
VERSION=0.0.34
TASK=

ANACONDA_BIN=anaconda3/bin
PIP=${HOME}/${ANACONDA_BIN}/pip
PYTHON=${HOME}/${ANACONDA_BIN}/python

all:

###########################
# Quantarhei installation #
###########################
install: inst


inst: sdist
	${PIP} install `ls dist/quantarhei-${VERSION}*`


################
# Distribution #
################
sdist:
	${PYTHON} setup.py sdist


##################
# Uninstallation #
##################
uninstall: uninst


uninst: 
	${PIP} uninstall -y quantarhei


##################
# Upload to pypi #
##################
upload: sdist
	twine upload `ls dist/quantarhei-${VERSION}*`


############
# Clean-up #
############
clean:
	rm -rf dist
	rm -rf quantarhei.egg-info
	rm -rf result_images
	rm -rf qrconf.py quantarhei/qrconf.py



################################
# Reinstallation with clean-up #
################################
reinst: clean uninst inst


#############################################################
# Local tests: this will reinstall Quantarhei and run tests #
#############################################################
local_tests: reinst
	paver  ${TASK}


#############################################
# Test of the install version of quantarhei #
#############################################
tests: 
	paver
	

####################
# Test of plotting #
####################
plot_tests: 
	paver matplotlib_tests


pylint:
	paver pylint

help:
	@echo "inst, reinst, local_tests, plot_tests, tests, sdist, clean"


sync:
	git fetch upstream
	git checkout master
	git merge upstream/master


