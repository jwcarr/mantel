import numpy as np
from scipy import stats

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


def plot(
    result,
    axis,
    significance_level=0.05,
    gaussian_background_color="blue",
    gaussian_background_alpha=0.1,
    gaussian_color="blue",
    gaussian_alpha=0.3,
    gaussian_curve_color="blue",
    gaussian_curve_alpha=0.3,
    hist_fill_color="orange",
    acceptance_color="green",
    rejection_color="red",
):
    """
    Plot the correlations previously computed (see
    compute_correlations() function) and on the given matplotlib axes,
    as well as the theoretical normal distribution with its confidence
    interval highlighted.

    Parameters
    ----------
    results : MantelResult
            The reuslt object output from mantel.test()
    axis : matplotlib.axes.Axes
            The matplotlib figure to draw on.
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
    """
    if plt is None:
        raise ImportError("Matplotlib is required for plotting")

    min_corr, max_corr = confidence_interval(result, significance_level)

    x = np.linspace(min_corr, max_corr, 100)
    y = stats.norm.pdf(x, result.mean, result.std)

    lower = -5 * result.std + result.mean
    upper = 5 * result.std + result.mean
    x_all = np.linspace(lower, upper, 100)
    y_all = stats.norm.pdf(x_all, result.mean, result.std)

    axis.fill_between(
        x_all,
        y_all,
        0,
        color=gaussian_background_color,
        alpha=gaussian_background_alpha,
    )
    axis.fill_between(x, y, 0, color=gaussian_color, alpha=gaussian_alpha)
    axis.plot(x_all, y_all, color=gaussian_curve_color, alpha=gaussian_curve_alpha)
    axis.vlines(
        x=[min_corr, max_corr],
        ymin=0,
        ymax=[
            stats.norm.pdf(min_corr, result.mean, result.std),
            stats.norm.pdf(max_corr, result.mean, result.std),
        ],
        linestyle="-",
        color=gaussian_curve_color,
        alpha=gaussian_curve_alpha,
    )
    axis.hist(
        result.correlations,
        bins=20,
        range=(lower, upper),
        density=True,
        histtype="stepfilled",
        color=hist_fill_color,
    )

    axis.set_xlim(left=lower, right=upper)
    axis.set_xlabel("Correlation coefficient")
    axis.set_yticks([])

    threshold_color = (
        acceptance_color
        if min_corr <= result.r <= max_corr
        else rejection_color
    )
    axis.axvline(x=result.r, linestyle=":", color=threshold_color)
    axis.annotate("{:.2f}".format(result.r), xy=(result.r, 0.9), color=threshold_color)


def confidence_interval(result, significance_level=0.05):
    """
    Return the confidence interval in a normal distribution
    (characterized by its given mean and its given standard deviation)
    for accepting the null hypothesis with some given significance
    level and some given tail method.
    Parameters
    ----------
    mean : float, optional
            The mean of the normal distribution (default: 0).
    std : float, optional
            The standard deviation of the normal distribution
            (default: 1).
    significance_level : float (between 0 and 1), optional
            The probability of rejecting the null hypothesis when it
            is true (default: 0.05).
    tail : str, optional
            Which tail to test in the calculation of the confidence
            interval; either 'upper', 'lower', or 'two-tail' (default:
            'two-tail').
    Returns
    -------
    lower_bound : float
            Lower bound of the confidence interval.
    upper_bound : float
            Upper bound of the confidence interval.
    """
    if not (0 < significance_level < 1):
        raise ValueError("The significance level must be in the range ]0, 1[")
    if result.tail == "upper":
        lower_bound = stats.norm.ppf(1e-10, result.mean, result.std)
        upper_bound = stats.norm.ppf(1 - significance_level, result.mean, result.std)
    elif result.tail == "lower":
        lower_bound = stats.norm.ppf(significance_level, result.mean, result.std)
        upper_bound = stats.norm.ppf(1 - 1e-10, result.mean, result.std)
    elif result.tail == "two-tail":
        lower_bound = stats.norm.ppf(significance_level / 2, result.mean, result.std)
        upper_bound = stats.norm.ppf(1 - significance_level / 2, result.mean, result.std)
    return (lower_bound, upper_bound)
