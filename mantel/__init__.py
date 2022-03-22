from ._test import compute_correlations, mantel_test_from_correlations, mantel_test, test

try:
    from ._version import __version__
except ImportError:
    __version__ = "???"
