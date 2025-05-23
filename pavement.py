"""
    Quantarhei Package Tests
    ========================
    
    To run the tests, the following packages have to be installed
    
    paver (pip install paver)
    nose (pip install nose)
    aloe (pip install aloe)
    behave (pip install behave)
    pylint (pip install pylint)
    
    About our tests
    ---------------
    
    Quantarhei is tested against unit tests run by the `nose` tool, against
    acceptance tests run by `aloe` tool and against doc tests (run by `nose`).
    We also calculate coverege using the `coverage` tool and we plan to 
    check the code by `pylint`.
    
    
    The ultimit goal of testing in 100 % coverage and perfect code style.
    
    
    The following tasks are available. Run them with `paver` as
    
    .. code:: bash
    
        $ paver TASK_NAME
        
    To test the whole suit, just run `paver` without any task.

    .. code:: bash

        $ paver
        
        
        
    Doctest Tasks
    -------------
    
    doc_tests
    doc_tests_v
    doc_tests_vs
    doc_tests_cov
    doc_tests_cov_v
    doc_tests_cov_vs
        
    Runs test doc tests. _cov means `` option, _v means
    with `-v` option for verbose output and _vs means with `-vs` option, i.e.
    verbose and printing standard output of the tests.
        
        
    Unit tests
    ----------
    
    unit_tests_vs
    unit_tests_cov_vs


    Acceptance tests
    ----------------
    
    aloe_tests_vs
    aloe_tests_cov_vs
    
    
    Tests to be run during development
    ----------------------------------
    
    doc_dev
    unit_dev
    
    Edit the list of files or directories that you want to test during
    development.
    

"""
import contextlib
import os
import subprocess
import platform


from paver.tasks import task
from paver.tasks import needs
from paver.easy import sh


""" The project version has to be set here and then in the following files:
    
    Makefile
    docs/sphinx/conf.py
    quantarhei/core/managers.py
    setup.py
  
    Ideally it would be all updated base on the present file by calling:
    
    > paver version-update
    
"""
version = "0.0.67"


sys_name = platform.system()



pip = 'pip'
python = 'python'

# 
# Commands for deleting files and directories silently and without error codes
#
if sys_name == "Darwin" or sys_name == "Linux":
    deldir = 'rm -r -f '
    delfile = 'rm -r -f '
elif sys_name == "Windows":
    deldir = 'rmdir /s /q '
    delfile = 'del /s /q '
else:
    raise Exception("Unknown system")

#
# The Mother of all repositories
#
repository = 'https://github.com/tmancal74/quantarhei'


#
# look for location of `behave`
#
if sys_name != "Windows":
    p = subprocess.Popen('which behave', shell=True,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    ii = 0
    for line in p.stdout.readlines():
        behave_bin = str(line.decode())
        behave_bin = behave_bin.strip()
        ii += 1
    if ii != 1:
        raise Exception("Don't know where `behave` is")
else:
    behave_bin = "behave"
    

#
# Context manager for getting into subdirectories
#
@contextlib.contextmanager
def cd(path):
   old_path = os.getcwd()
   os.chdir(path)
   try:
       yield
   finally:
       os.chdir(old_path)
       
def rm_rf(path):
    """Removal of files and directories, recursively and silently
    
    """
    pass
    
#
###############################################################################
#
#      Standard developer tasks
#
###############################################################################
#      
###############################################################################
#  Creates distribution
###############################################################################   
@task
def sdist():
    sh(python+" setup.py sdist")

###############################################################################
# Installs Quantarhei from the source
###############################################################################
@needs('sdist')
@task
def inst():
    
    sh(pip+' install dist/quantarhei-'+version+'.tar.gz')
    
#
#  The same as above
#
@needs('inst')
@task
def install():
    pass

###############################################################################
#  Uninstall local installation of Quantarhei
###############################################################################
@needs('uninst')
@task
def uninstall():
    pass

#
# The same as above  
#
@task
def uninst():
	sh(pip+' uninstall -y quantarhei')
    

###############################################################################
#  Upload to pypi 
###############################################################################
@needs('sdist')
@task
def upload():
	sh('twine upload dist/quantarhei-'+version+'.tar.gz')


###############################################################################
#  Clean-up 
###############################################################################
@task
def clean():
    try:
        sh(deldir+'dist')
    except:
        print("Directory not present - no problem")
    try:    
        sh(delfile+'quantarhei.egg-info')
    except:
        print("File not present - no problem")
    try:
        sh(deldir+'result_images')
    except:
        print("Directory not present - no problem")
    try:
        sh(delfile+'qrconf.py quantarhei/qrconf.py')
    except:
        print("File not present - no problem")
    try:        
        sh(delfile+'test.log')
    except:
        print("File not present - no problem")


###############################################################################
# Reinstallation with clean-up 
###############################################################################
@needs('clean','uninst', 'inst')
@task
def reinst():
    pass

###############################################################################
# Local tests: this will reinstall Quantarhei and run tests 
###############################################################################
@needs('reinst')
@task
def local_tests():
	sh('paver')

###############################################################################
# Test of plotting 
###############################################################################
@needs('matplotlib_tests')
@task
def plot_tests(): 
	pass 

###############################################################################
# Update examples  
###############################################################################
@task
def update_examples():
	sh('cd examples; python admin/make_demos.py')

###############################################################################
#
###############################################################################
@task
def tasks():
	print("")
	print("Quantarhei Paver Tasks")
	print("======================")
	print("")
	print("Essential tasks: ")
	print("----------------")
	print("inst        ... install quantarhei from this source code")
	print("reinst      ... uninstall and install from this source code")
	print("local_tests ... uninstall, install and run tests")
	print("plot_tests  ... run tests of plotting")
	print("test        ... run tests")
	print("sdist       ... create distribution")
	print("clean       ... clean the repository")
	print("examples    ... updates examples")
	print("")
	print("Git tasks: ")
	print("----------")
	print("git_add_upstream, git_update_master")
	print("")

###############################################################################
#  Helper tasks for managing pull-requests
###############################################################################

#
# update from master branch of the quantarhei's main repository
#
@task
def git_update_master():
	sh('git fetch upstream')
	sh('git checkout master')
	sh('git merge upstream/master')


###############################################################################
# connect a forked local repository to the main quantarhei repository	
###############################################################################
@task
def git_add_upstream():
	sh('git remote add upstream '+repository)



###############################################################################
# 
#     Graphical output testing
#
###############################################################################

#
# TASK: matplotlib_tests
#
@task
def matplotlib_tests():
    """Tests matplotlib output
    
    This test seems not to work on all platforms. 
    
    !!!Currently it is not run in automatic build testing!!!
    
    """
    sh('cd tests/matplotlib; ./tests.py; cd ..')

#
# TASK: remove .coverage
#
@task
def coverage_erase():
    """ Remove coverage from previous build

    """
    sh('coverage erase')


#
# TASK: combine coverage reports
#
@task
def coverage_combine():
    """ Combines all coverage reports into one
    
    """
    sh('coverage combine --append')
    

###############################################################################
#    
# Tests of Quantarhei installed examples
#
###############################################################################

cover_flags = ' '
covr = '  ' 
nose_test_dir = 'tests/unit/wizard/examples/ '

#
# TASK: examples
#
@task
def examples():
    """Tests if Quantarhei examples run correctly
    
    """
    sh('coverage run -m nose  -vs '+covr+cover_flags+nose_test_dir)

    
###############################################################################
#
# Unit tests
#
###############################################################################

nose_test_dir = 'tests/unit'

#
# TASK: unit_tests_vs
#
@task
def unit_tests_vs():
    """Performs verbose units tests printing output but without coverage 
    
    """
    sh('pytest -vs --cov=quantarhei --cov-report=term --cov-report=xml --cov-report=html '+cover_flags+nose_test_dir)

#
# TASK: unit_tests_v
#
@task
def unit_tests_v():
    """Performs verbose units tests capturing output but without coverage 
    
    """
    sh('pytest -v --cov=quantarhei --cov-report=term --cov-report=xml --cov-report=html '+cover_flags+nose_test_dir)

#
# TASK: unit_tests_cov_vs
#
@task
def unit_tests_cov_vs():
    """Performs verbose units tests without capturing output and with coverage 
    
    """
    sh('pytest -vs --cov=quantarhei --cov-report=term --cov-report=xml --cov-report=html '+covr+cover_flags+nose_test_dir)    

#
# TASK: unit_tests_cov_v
#
@task
def unit_tests_cov_v():
    """Performs verbose units tests capturing output and with coverage 
    
    """
    sh('pytest -v --cov=quantarhei --cov-report=term --cov-report=xml --cov-report=html '+covr+cover_flags+nose_test_dir)
 


###############################################################################
#
#   Doc tests
#
###############################################################################

# Doc tests with coverage
cmdd = 'pytest --doctest-modules'

# Directories to test
dirs = ['core', 'builders', 'qm/corfunctions', 'spectroscopy',
        'qm/liouvillespace',  'functions', 'qm/hilbertspace', 'qm/propagators']

# Specific modules to test
mods = ['qm/propagators/poppropagator.py']

#
# TASK: doc_tests_cov_vs
#
@task
def doc_tests_cov_vs():    
    """Verbose doc test with coverage
    
    This test runs with coverage in a verbose mode. Output is not captured
    but printed to the screen.
 
    Run it as
    
    .. code:: bash
        
        paver doc_tests_cov_vs
        
    """
    tcmd = cmdd+covr+' -vs '+cover_flags
    _run_tests(tcmd, mods, dirs)
    
    
#
# TASK: doc_tests_cov_vs
#    
@task
def doc_tests_cov_v():
    """Verbose doc test with coverage
    
    This test runs with coverage in a verbose mode. Output is captured.
 
    Run it as
    
    .. code:: bash
        
        paver doc_tests_cov_v
        
    """
    tcmd = cmdd+covr+' -v '+cover_flags
    _run_tests(tcmd, mods, dirs)

#
# TASK: doc_tests_vs
#
@task
def doc_tests_vs():    
    """Verbose doc test without coverage
    
    This test runs without coverage in a verbose mode. Output is not captured
    but printed to the screen.
 
    Run it as
    
    .. code:: bash
        
        paver doc_tests_vs
        
    """
    tcmd = cmdd+' -vs '+cover_flags
    _run_tests(tcmd, mods, dirs)

#
# TASK: doc_tests_vs
#
@task
def doc_tests_v():    
    """Verbose doc test without coverage, capturing the output
    
    This test runs without coverage in a verbose mode. Output is captured.
 
    Run it as
    
    .. code:: bash
        
        paver doc_tests_v
        
    """
    tcmd = cmdd+' -v '+cover_flags
    _run_tests(tcmd, mods, dirs) 


def _run_tests(tcmd, dirs, mods):
    """Runs tests on directories and modules
    
    """
    
    # run on directories
    for _dir in dirs:
        sh(tcmd+"quantarhei/"+_dir)

    # run on modules
    for _mod in mods:
        sh(tcmd+"quantarhei/"+_mod)
    
    
    
###############################################################################
#
#   Aloe tests
#
###############################################################################

aloe_test_dir = 'tests/bdd'

@task
def aloe_tests_vs():
    sh("aloe -vs -a !in_development "+cover_flags+aloe_test_dir)

@task
def aloe_tests_v():
    sh("aloe -v -a !in_development "+cover_flags+aloe_test_dir)
    
@task
def aloe_tests_cov_vs():
    sh("aloe  -vs -a !in_development "+cover_flags+aloe_test_dir)

@task
def aloe_tests_cov_v():
    sh('aloe  -v -a !in_development '+cover_flags+aloe_test_dir)


###############################################################################
#
#    Behave tests
#
###############################################################################
@task
def behave():
    path = os.path.join('.','tests', 'behave', 'features')
    sh("coverage run "+behave_bin+" "+path)

###############################################################################
#
#   Documentation build test
#
###############################################################################
    
#
# TASK: html
#
@task
def html():
    sh('cd docs/sphinx/; make html')

#
# TASK: html_clean
#   
@task
def html_clean():
    sh('cd docs/sphinx/; make clean')

    
###############################################################################
#
#   Pylint tests
#
###############################################################################
class Runner:  
    """Class to run pylint tests
    
    """
    def __init__(self,path):
        self.path = path
        
    def set_path(self,path):
        self.path = path
        
    def un(self,file):
        rcfile = 'tests/pylint/pylintrc'
        sh('pylint --rcfile='+rcfile+' '+self.path+'/'+file)
    
@task
def pylint():
    """Runs pylint tests on specified files
    
    """

    path = 'quantarhei/core/'
    r = Runner(path)
#    r.un('valueaxis.py')
#    r.un('time.py')
#    r.un('frequency.py')
#    r.un('dfunction.py')

    path = 'quantarhei/qm/corfunctions'
    r.set_path(path)
#    r.un('correlationfunctions.py')
#    r.un('spectraldensities.py')  
#    r.un('cfmatrix.py')

    path = 'quantarhei/qm/hilbertspace'
    r.set_path(path)
#    r.un('statevector.py')
#    r.un('aggregate_test.py')

    path = 'quantarhei/qm/liouvillespace'
    r.set_path(path)
#    r.un('lindbladform.py')
    
    path = 'quantarhei/builders'
    r.set_path(path)
#    r.un('pdb.py')    

    path = 'quantarhei/core'
    r.set_path(path)
    #r.un('matrixdata.py')
    #r.un('saveable.py')    
    
    path = 'quantarhei/utils/'
    r.set_path(path)
#    r.un('logging.py')
    
    path = 'quantarhei/scripts/'
    r.set_path(path)
    #r.un('ghenerate.py')
    r.un('qrhei.py')


###############################################################################
#
#   Development tasks
#
#   Edit this if you need to test only particular part
#   of the code which you develop. Then maybe ignore this file during commits
#
###############################################################################
  
#
# Targets for doc tests
#       
doc_devs = ['quantarhei/spectroscopy/twod2.py', 
            'quantarhei/qm/liouvillespace',
            'quantarhei/functions']

#
# Targets for unit tests
#
unit_devs = ['tests/unit/qm/liouvillespace/test_evolutionsuperoperator.py',
             'tests/unit/qm/propagators/rdmpropagator_test.py']


cmd = 'pytest'
doct = ''

#
# TASK: doc_dev
#
@task
def doc_dev():
   """Task to run doctests on files under development
   
   """
   tcmd = cmd+doct+' -vs ' 
   for _dev in doc_devs:
       sh(tcmd+_dev)

#
# TASK: unit_dev
#
@task
def unit_dev():
   """Task to run unit tests on files under development
   
   """
   tcmd = cmd+' -vs ' 
   for _dev in unit_devs:
       sh(tcmd+_dev)
       

###############################################################################
#  
#  Codecov (upload of coverage results to codecov.io server)
#
###############################################################################
   
#
# TASK: codecov
#
@task
def codecov():
    """Uploads the results of `coverage` command on the codecov.io server
    
    """
    sh('codecov --token=1d501d6b-b32d-44c2-8253-8b12ef3c8cc9')


###############################################################################
#
#   Collections of tests including default test task
#
###############################################################################

# Removing matplotlib tests because it does not work on remote travis CI
#@needs('matplotlib_tests',
#       'unit_tests_cov_vs',
#       'doc_tests_cov_vs',
#       'aloe_tests_cov_vs',
#       'pylint')
      
#
# This is called when paver is run without any task
#
#@needs('unit_tests_cov_v',
#       'doc_tests_cov_v',
#       'aloe_tests_cov_v',
#       'behave')
@needs('test')
@task
def default():
    """Default paver task
    
    """
    pass

#
# This is called when paver is run without any task
#

#@needs('coverage_erase',
#       'unit_tests_cov_v',
#       'doc_tests_cov_v',
#       'aloe_tests_cov_v',
#       'behave',
#       'coverage_combine')
@needs('coverage_erase',
       'unit_tests_cov_v',
       'doc_tests_cov_v',
       'behave')
@task
def test():
    """Default paver task
    
    """
    pass

#
# Verbose version of the default tests
#
#@needs('coverage_erase',
#       'unit_tests_cov_vs',
#       'doc_tests_cov_vs',
#       'aloe_tests_cov_vs')
@needs('coverage_erase',
       'unit_tests_cov_vs',
       'doc_tests_cov_vs')
@task
def verbose():
    """Run the default tests but with printing output
    
    """
    pass


