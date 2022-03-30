from ._test import test
from ._plot import plot

try:
    from ._version import __version__
except ImportError:
    __version__ = "???"
