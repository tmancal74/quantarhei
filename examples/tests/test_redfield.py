# -*- coding: utf-8 -*-

from tests.qm.liouvillespace.test_redfield import TestRedfield


t1 = TestRedfield()
t1.setUp(verbose=True)

t1.test_comparison_of_rates()

