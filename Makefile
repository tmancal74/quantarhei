.PHONY: sync lint test unit doctest behave all clean-coverage

sync:
	uv sync --extra dev

lint:
	uv run ruff format --check .
	uv run pre-commit run --all-files

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
		src/quantarhei/core \
		src/quantarhei/builders \
		src/quantarhei/qm/corfunctions \
		src/quantarhei/spectroscopy \
		src/quantarhei/qm/liouvillespace \
		src/quantarhei/functions \
		src/quantarhei/qm/hilbertspace \
		src/quantarhei/qm/propagators \
		src/quantarhei/qm/propagators/poppropagator.py

behave:
	uv run coverage run --data-file=.coverage.behave \
		$$(uv run which behave) tests/behave/features
	uv run coverage combine --data-file=.coverage.behave
	uv run coverage xml --data-file=.coverage.behave -o coverage3.xml

test: unit doctest behave

all: sync lint test

clean-coverage:
	rm -f .coverage .coverage.* coverage*.xml
	rm -rf htmlcov