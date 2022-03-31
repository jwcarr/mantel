import numpy as np
from scipy import stats
from ._test import MantelResult

try:
    import matplotlib.pyplot as plt
    from matplotlib.axes import Axes
except ImportError:
    plt = None


def plot(
    result,
    axis=None,
    alpha=0.05,
    hist_color="lightgray",
    gaussian_color="black",
    acceptance_color="black",
    rejection_color="black",
):
    """
    Plot a histogram of the sample correlations from a MantelResult object on
    a Matplotlib axis. If no axis is provided, a new plot is created and
    returned. The theoretical normal distribution (and confidence interval)
    is also shown, along with the veridical correlation.

    Parameters
    ----------
    result : MantelResult
        The reuslt object output from mantel.test()
    axis : matplotlib.axes.Axes
        The matplotlib figure to draw on.
    alpha : float
        Significance level for rejecting the null hypothesis (default: 5%)
    hist_color: str
        Color used for the histogram bars (default: 'lightgray').
    gaussian_color: str
        Color used for the normal distribution curve and the confidence
        interval limits (default: 'black').
    acceptance_color: str
        Color used for drawing the vertical line and the label of the
        veridical correlation if the null hypothesis is rejected according to
        the significance level value (default: 'black').
    rejection_color: str
        Color used for drawing the vertical line and the label of the
        veridical correlation if the null hypothesis cannot be rejected
        according to the significance level value (default: 'black').

    Returns
    -------
    Matplotlib Figure and Axis objects
    """
    if plt is None:
        raise ImportError("Matplotlib is required for plotting")
    if not isinstance(result, MantelResult):
        raise TypeError("result shoud be a MantelResult object")
    if axis is None:
        fig, axis = plt.subplots()
    elif isinstance(axis, Axes):
        fig = axis.get_figure()
    else:
        raise TypeError("axis should be a Matplotlib Axis object")

    lower = max(-5 * result.std + result.mean, -1)
    upper = min(5 * result.std + result.mean, 1)

    n_bins = 100 if result.perms >= 1000 else 20
    axis.hist(
        result.correlations,
        bins=n_bins,
        range=(lower, upper),
        density=True,
        histtype="stepfilled",
        color=hist_color,
    )

    x = np.linspace(lower, upper, 200)
    y = stats.norm.pdf(x, result.mean, result.std)
    axis.plot(x, y, color=gaussian_color)

    min_corr, max_corr = confidence_interval(result, alpha)
    axis.vlines(
        x=[min_corr, max_corr],
        ymin=0,
        ymax=[
            stats.norm.pdf(min_corr, result.mean, result.std),
            stats.norm.pdf(max_corr, result.mean, result.std),
        ],
        linestyle="-",
        color=gaussian_color,
    )

    threshold_color = (
        rejection_color if min_corr <= result.r <= max_corr else acceptance_color
    )
    axis.axvline(x=result.r, linestyle=":", color=threshold_color)
    axis.annotate("{:.2f}".format(result.r), xy=(result.r, 0.9), color=threshold_color)

    axis.set_xlim(lower, upper)
    axis.set_xlabel("Correlation coefficient")
    axis.set_yticks([])
    return fig, axis


def confidence_interval(result, alpha=0.05):
    """
    Return the confidence interval in a normal distribution(characterized by
    its given mean and its given standard deviation) for rejecting the null
    hypothesis with some given significance level and some given tail
    method.

    Parameters
    ----------
    result : float
        MantelResult object
    alpha : float (between 0 and 1), optional
        The probability of rejecting the null hypothesis when it is true
        (default: 0.05).

    Returns
    -------
    lower_bound : float
        Lower bound of the confidence interval.
    upper_bound : float
        Upper bound of the confidence interval.
    """
    if not (0 < alpha < 1):
        raise ValueError("The alpha level must be in the interval (0, 1)")
    if result.tail == "upper":
        lower_bound = stats.norm.ppf(1e-10, result.mean, result.std)
        upper_bound = stats.norm.ppf(1 - alpha, result.mean, result.std)
    elif result.tail == "lower":
        lower_bound = stats.norm.ppf(alpha, result.mean, result.std)
        upper_bound = stats.norm.ppf(1 - 1e-10, result.mean, result.std)
    elif result.tail == "two-tail":
        lower_bound = stats.norm.ppf(alpha / 2, result.mean, result.std)
        upper_bound = stats.norm.ppf(1 - alpha / 2, result.mean, result.std)
    return lower_bound, upper_bound
