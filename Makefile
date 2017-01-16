all: install 

VERSION = 0.0.11.dev1

uninst:
	pip uninstall quantarhei

install:
	python setup.py sdist
	pip install dist/quantarhei-${VERSION}.tar.gz

reinst: uninst install
	@echo "Reinstallation"
    
