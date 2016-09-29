from paver.tasks import task
from paver.tasks import needs
from paver.easy import sh

""" 
    Particular tasks
"""

@task
def unit_tests_vs():
    sh('nosetests -vs tests/unit')
@task
def doc_tests_vs():    
    sh('nosetests --with-doctest -vs quantarhei/core')

@task
def unit_tests_v():
    sh('nosetests -v tests/unit')
    
@task
def doc_tests_v():    
    sh('nosetests --with-doctest -vs quantarhei/core')
    
@task
def aloe_tests_vs():
    sh("aloe -vs -a !in_development tests/bdd")

@task
def aloe_tests_v():
    sh('aloe -v -a !in_development tests/bdd')
    

"""
    Optional tasks
"""    
@task
def pylint():
    sh('pylint --rcfile=pylint/pylintrc quantarhei/core/valueaxis.py')
    sh('pylint --rcfile=pylint/pylintrc quantarhei/core/time.py')
    sh('pylint --rcfile=pylint/pylintrc quantarhei/core/frequency.py')
    
"""
    Default 
"""

@needs('unit_tests_v','doc_tests_v','aloe_tests_v')
@task
def default():
    pass

@needs('unit_tests_vs','doc_tests_vs','aloe_tests_vs')
@task
def vs():
    pass

@needs('unit_tests_v','aloe_tests_v')
@task
def nodoc():
    pass

@needs('unit_tests_vs','aloe_tests_vs')
@task
def nodoc_vs():
    pass

@needs('unit_tests')
@task
def windows():
    """ On windows, aloe tool does not work. We do only unit tests"""
    pass

