.. _write-tests: 

How to Write Tests for Quantarhei
=================================

Every contribution to Quantarhei has to come with appropriate tests. Because
of Quantarhei's development history, less than a half of current code is
subject to testing at every build. This will change in near future.

Running Tests
-------------

To run the tests on Quantarhei current installation, just invoke the `paver`
command (`paver` package is listed as one of the testing dependencies):

.. code:: bash

    $ paver
    
Equivalently, you can use the predefined task *test* the Makefile as:

.. code:: bash

    $ make test
    
To run the same tests, but with the output of the scipts printed on the screen
(normally such output is camptured and only the information about tests is
printed), you can run:

.. code:: bash

    $ paver verbose
    
This runs the same set of tests, except for the capturing option removed.
    
Sometimes you want to make sure that tests are run against your latest code
(it is assumed here that tests run against the installed copy of Quantarhei).
For this purpose we have a Makefile task *local_tests*, which rebuilds and
reinstalls Quantarhei, before running tests:

.. code:: bash

    $ make local_tests
    
Appart from these predefined tests, you can also run other versions of the
tests or just partial tests, e.g. on the files that you develop. Consult the
`pavement.py` file in the root directory of Quantarhei package.

Writting Doc Tests
------------------

Writting Unit Tests
-------------------

Writting Acceptance Tests
-------------------------

