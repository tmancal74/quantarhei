from paver.tasks import task
from paver.tasks import needs
from paver.easy import sh

@task
def unit_tests_vs():
    sh('nosetests -vs tests/unit')

@task
def aloe_tests_vs():
    sh("aloe -vs -a !in_development tests/bdd")

@task
def unit_tests_v():
    sh('nosetests -v tests/unit')

@task
def aloe_tests_v():
    sh('aloe -v -a !in_development tests/bdd')

@needs('unit_tests_v','aloe_tests_v')
@task
def default():
    pass


@needs('unit_tests')
@task
def windows():
    """ On windows, aloe tool does not work. We do only unit tests"""
    pass

