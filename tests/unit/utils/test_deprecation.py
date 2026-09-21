import warnings

from quantarhei.utils.deprecation import deprecated


class TestDeprecated:
    def test_emits_deprecation_warning(self):
        @deprecated("Use bar() instead.")
        def foo():
            return 42

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = foo()

        assert result == 42
        assert len(w) == 1
        assert w[0].category is DeprecationWarning
        assert "foo is deprecated" in str(w[0].message)
        assert "Use bar() instead." in str(w[0].message)

    def test_preserves_function_metadata(self):
        @deprecated("Gone.")
        def documented_func():
            """My docstring."""
            return 1

        assert documented_func.__name__ == "documented_func"
        assert documented_func.__doc__ == "My docstring."

    def test_passes_arguments(self):
        @deprecated("Use add() instead.")
        def legacy_add(a, b, offset=0):
            return a + b + offset

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            assert legacy_add(2, 3, offset=10) == 15

    def test_works_on_methods(self):
        class MyClass:
            @deprecated("Use new_method().")
            def old_method(self):
                return "old"

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = MyClass().old_method()

        assert result == "old"
        assert "old_method is deprecated" in str(w[0].message)
