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
VERSION=0.0.39
TASK=

ANACONDA_BIN=anaconda3/bin


PIP=${HOME}/${ANACONDA_BIN}/pip
PYTHON=${HOME}/${ANACONDA_BIN}/python

REPOSITORY=https://github.com/tmancal74/quantarhei

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
	rm -rf coverage.xml


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
test: 
	paver
	

####################
# Test of plotting #
####################
plot_tests: 
	paver matplotlib_tests

####################
# Update examples  #
####################
update_examples:
	cd examples; python admin/make_demos.py


pylint:
	paver pylint


help:
	@echo ""
	@echo "Quantarhei Makefile Help"
	@echo "========================"
	@echo ""
	@echo "Essential tasks: "
	@echo "----------------"
	@echo "inst        ... install quantarhei from this source code"
	@echo "reinst      ... uninstall and install from this source code"
	@echo "local_tests ... uninstall, install and run tests"
	@echo "plot_tests  ... run tests of plotting"
	@echo "test        ... run tests"
	@echo "sdist       ... create distribution"
	@echo "clean       ... clean the repository"
	@echo "examples    ... updates examples"
	@echo ""
	@echo "Git tasks: "
	@echo "----------"
	@echo "git_add_upstream, git_update_master"
	@echo ""

############################################
#  Helper tasks for managing pull-requests
############################################

#
# update from master branch of the quantarhei's main repository
#
git_update_master:
	git fetch upstream
	git checkout master
	git merge upstream/master


#
# connect a forked local repository to the main quantarhei repository	
#
git_add_upstream:
	git remote add upstream ${REPOSITORY}


