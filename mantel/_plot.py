from ._utils import \
    probability_from_sample as _probability_from_sample, \
    confidence_interval as _confidence_interval
import numpy as np
from scipy import stats
try:
    import matplotlib.pyplot as plt
    _has_matplotlib = True
except ImportError:
    _has_matplotlib = False


def plot(correlations, plot,
                      tail = "two-tail",
                      significance_level = 0.05,
                      gaussian_background_color = 'blue',
                      gaussian_background_alpha = 0.1,
                      gaussian_color = 'blue',
                      gaussian_alpha = 0.3,
                      gaussian_curve_color = 'blue',
                      gaussian_curve_alpha = 0.3,
                      hist_fill_color = 'orange',
                      hist_edge_color = 'green',
                      hist_alpha = 0.7,
                      acceptance_color = 'green',
                      rejection_color = 'red'):
    """
    Plot the correlations previously computed (see
    compute_correlations() function) and on the given matplotlib axes,
    as well as the theoretical normal distribution with its confidence
    interval highlighted.

    Parameters
    ----------
    correlations : array of floats
            The correlations computed by the Mantel test.
    plot : matplotlib.axes.Axes
            The matplotlib figure to draw on.
    tail : str, optional
            Which tail to test in the calculation of the empirical p-value; either
            'upper', 'lower', or 'two-tail' (default: 'two-tail').
    significance_level : float
            Significance level of the null hypothesis of the Mantel
            test (default: 5%)
    gaussian_background_color: str
            Color used for painting the area under the normal
            distribution curve (default: 'blue').

    gaussian_background_alpha: float (between [0, 1])
            Opacity (alpha channel) of the area under the normal
            distribution curve (default: 0.1).

    gaussian_color: str
            Color used for painting the area of the confidence
            interval under the normal distribution curve (default:
            'blue').
    gaussian_alpha: float (between [0, 1])
            Opacity (alpha channel) of the area of the confidence
            interval under the normal distribution curve (default:
            0.3).
    gaussian_curve_color: str
            Color used for painting the normal distribution curve and
            the confidence interval limits (default: 'blue').
    gaussian_curve_alpha: float (between [0, 1])
            Opacity (alpha channel) of the normal distribution curve
            and the confidence interval limits (default: 0.3).
    hist_fill_color: str
            Color used for filling the correlations histogram bars
            (default: 'orange').
    hist_alpha: float (between [0, 1])
            Opacity (alpha channel) of the correlations histogram bars
            (default: 0.7).
    hist_edge_color: str
            Color used for drawing the correlations histogram bar
            edges (default: 'green').
    acceptance_color: str
            Color used for drawing the vertical line and the label of
            the veridical correlation if the null hypothesis is
            accepted according to the significance level value
            (default: 'green').
    rejection_color: str
            Color used for drawing the vertical line and the label of
            the veridical correlation if the null hypothesis is
            rejected according to the significance level value
            (default: 'red').
    Returns
    -------
    min_correlation : float
            Minimal correlation value to accept the null hypothesis.
    max_correlation : float
            Maximal correlation value to accept the null hypothesis.
    """
    if not _has_matplotlib:
        raise Exception("In order to produce histograms, you need to install the 'matplotlib' library first.")

    r = correlations[0]
    m = np.mean(correlations)
    s = np.std(correlations)
    p = _probability_from_sample(value = r,
                                 sample = correlations,
                                 tail = tail)

    (min_correlation, max_correlation) = _confidence_interval(
        mean = m, std = s,
        significance_level = significance_level,
        tail = tail
    )

    x = np.linspace(min_correlation, max_correlation, 100)
    y = stats.norm.pdf(x, m, s)

    lower = -5 * s + m
    upper = 5 * s + m
    x_all = np.linspace(lower, upper, 100)
    y_all = stats.norm.pdf(x_all, m, s)

    plot.fill_between(x_all, y_all, 0, color=gaussian_background_color, alpha=gaussian_background_alpha)
    plot.fill_between(x, y, 0, color=gaussian_color, alpha=gaussian_alpha)
    plot.plot(x_all, y_all, color=gaussian_curve_color, alpha=gaussian_curve_alpha)
    plot.vlines(x=[min_correlation, max_correlation], ymin=0, ymax=[stats.norm.pdf(min_correlation, m, s), stats.norm.pdf(max_correlation, m, s)], linestyle='-', color=gaussian_curve_color, alpha=gaussian_curve_alpha)

    plot.hist(correlations, bins=20, range=(lower, upper), density=True, histtype='stepfilled', color=hist_fill_color, edgecolor=hist_edge_color, alpha=hist_alpha)

    plot.set_xlim(left=lower, right=upper)
    plot.set_xlabel("correlation coefficients")

    plot.set_ylabel("Density")

    threshold_color = acceptance_color if min_correlation <= r <= max_correlation else rejection_color
    plot.axvline(x=r, linestyle=':', color=threshold_color)
    plot.annotate('{:.2f}'.format(r), xy=(r, 0.9), color=threshold_color)
    return (min_correlation, max_correlation)
