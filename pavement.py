"""
    Quantarhei Package Tests
    ========================
    
    To run the tests, the following packages have to be installed
    
    paver (pip install paver)
    nose (pip install nose)
    aloe (pip install aloe)
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
        
    Runs test doc tests. _cov means `--with-coverage` option, _v means
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
from paver.tasks import task
from paver.tasks import needs
from paver.easy import sh


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


###############################################################################
#    
# Tests of Quantarhei installed examples
#
###############################################################################

#
# TASK: examples
#
@task
def examples():
    """Tests if Quantarhei examples run correctly
    
    """
    sh('nosetests  -vs tests/unit/wizard/examples/')

    
###############################################################################
#
# Unit tests
#
###############################################################################

#
# TASK: unit_tests_vs
#
@task
def unit_tests_vs():
    """Performs verbose units tests printing output but without coverage 
    
    """
    sh('nosetests -vs tests/unit')

#
# TASK: unit_tests_v
#
@task
def unit_tests_v():
    """Performs verbose units tests capturing output but without coverage 
    
    """
    sh('nosetests -v tests/unit')

#
# TASK: unit_tests_cov_vs
#
@task
def unit_tests_cov_vs():
    """Performs verbose units tests without capturing output and with coverage 
    
    """
    sh('nosetests  --with-coverage -vs tests/unit')    

#
# TASK: unit_tests_cov_v
#
@task
def unit_tests_cov_v():
    """Performs verbose units tests capturing output and with coverage 
    
    """
    sh('nosetests  --with-coverage -v tests/unit')     


###############################################################################
#
#   Doc tests
#
###############################################################################

# Doc tests with coverage
cmdd = 'nosetests --with-doctest'
covr = ' --with-coverage'

# Directories to test
dirs = ['core', 'builders', 'qm/corfunctions', 'spectroscopy',
        'qm/liouvillespace']

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
    tcmd = cmdd+covr+' -vs '
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
    tcmd = cmdd+covr+' -v '
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
    tcmd = cmdd+' -vs '
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
    tcmd = cmdd+' -v '
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
@task
def aloe_tests_vs():
    sh("aloe -vs -a !in_development tests/bdd")

@task
def aloe_tests_v():
    sh("aloe -v -a !in_development tests/bdd")
    
@task
def aloe_tests_cov_vs():
    sh('aloe --with-coverage -vs -a !in_development tests/bdd')

@task
def aloe_tests_cov_v():
    sh('aloe --with-coverage -v -a !in_development tests/bdd')


###############################################################################
#
#    Behave tests
#
###############################################################################
@task
def behave():
    sh("cd tests/behave/features; coverage run $(which behave)")


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
    r.un('ghenerate.py')


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


cmd = 'nosetests'
doct = ' --with-doctest'

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
    sh('codecov --token=34f2053d-7aa9-4b0c-b776-92d852e597ca')


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
@needs('unit_tests_cov_v',
       'doc_tests_cov_v',
       'aloe_tests_cov_v',
       'behave')
@task
def default():
    """Default paver task
    
    """
    pass

#
# Verbose version of the default tests
#
@needs('unit_tests_cov_vs',
       'doc_tests_cov_vs',
       'aloe_tests_cov_vs')
@task
def verbose():
    """Run the default tests but with printing output
    
    """
    pass