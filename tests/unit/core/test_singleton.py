import threading
import unittest

from quantarhei.core.singleton import Singleton


class _TestSingleton(metaclass=Singleton):
    def __init__(self):
        self.value = 0


class TestSingletonThreadSafety(unittest.TestCase):
    def setUp(self):
        Singleton._instances.pop(_TestSingleton, None)

    def test_concurrent_instantiation_returns_same_instance(self):
        instances = []
        barrier = threading.Barrier(10)

        def create():
            barrier.wait()
            instances.append(_TestSingleton())

        threads = [threading.Thread(target=create) for _ in range(10)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        self.assertEqual(len(set(id(i) for i in instances)), 1)

    def test_constructor_runs_exactly_once_under_concurrency(self):
        call_count = []
        barrier = threading.Barrier(10)

        original_init = _TestSingleton.__init__

        def counting_init(self):
            call_count.append(1)
            original_init(self)

        _TestSingleton.__init__ = counting_init
        try:

            def create():
                barrier.wait()
                _TestSingleton()

            threads = [threading.Thread(target=create) for _ in range(10)]
            for t in threads:
                t.start()
            for t in threads:
                t.join()

            self.assertEqual(len(call_count), 1)
        finally:
            _TestSingleton.__init__ = original_init
