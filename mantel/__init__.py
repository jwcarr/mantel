from ._test import compute_correlations, test
from ._plot import plot_correlations
from ._mantel import Mantel

try:
    from ._version import __version__
except ImportError:
    __version__ = "???"
