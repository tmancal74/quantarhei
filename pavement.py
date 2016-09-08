from paver.tasks import task
from paver.tasks import needs
from paver.easy import sh

@task
def unit_tests():
    sh('nosetests -vs tests/unit')

@task
def aloe_tests():
    sh("aloe -vs -a !in_development tests/bdd")


@needs('unit_tests','aloe_tests')
@task
def default():
    pass


@needs('unit_tests')
@task
def windows():
    """ On windows, aloe tool does not work. We do only unit tests"""
    pass

