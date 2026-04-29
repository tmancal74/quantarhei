.PHONY: sync lint lint-fix test unit doctest behave all clean-coverage

all: sync lint test


#
# Environment setup
#

sync:
	uv sync --extra dev


#
# Lint tests and fixes
#

lint:
	uv run ruff format --check .
	uv run pre-commit run --all-files

lint-fix:
	uv run ruff format .
	uv run ruff check quantarhei/ --fix
	uv run ruff check examples/ --fix

pre-commit: lint

# 
# Package tests
#

unit:
	uv run pytest -v \
		--cov=quantarhei \
		--cov-report=term \
		--cov-report=xml:coverage1.xml \
		--junit-xml=results-unit.xml \
		tests/unit

doctest:
	uv run pytest --doctest-modules -v \
		--cov=quantarhei \
		--cov-append \
		--cov-report=xml:coverage2.xml \
		--junit-xml=results-doctest.xml \
		quantarhei/core \
		quantarhei/builders \
		quantarhei/qm/corfunctions \
		quantarhei/spectroscopy \
		quantarhei/qm/liouvillespace \
		quantarhei/functions \
		quantarhei/qm/hilbertspace \
		quantarhei/qm/propagators \
		quantarhei/qm/propagators/poppropagator.py

behave:
	uv run coverage run --data-file=.coverage.behave \
		$$(uv run which behave) tests/behave/features
	uv run coverage combine --data-file=.coverage.behave
	uv run coverage xml --data-file=.coverage.behave -o coverage3.xml

test: unit doctest behave


#
# After test clean-up
#

clean-coverage:
	rm -f .coverage .coverage.* coverage*.xml
	rm -rf htmlcov

clean: clean-coverage
	rm -f results-*.xml
	rm -f test.log

