from itertools import permutations
import numpy as np
from scipy import spatial, stats
try:
    import matplotlib.pyplot as plt
    _has_matplotlib = True
except ImportError:
    _has_matplotlib = False


def compute_correlations(X, Y, perms=10000, method="pearson", ignore_nans=False):
    """
    Takes two distance matrices (either redundant matrices or
    condensed vectors) and computes the correlation between symmetric
    permutations two distance matrices.

    Parameters
    ----------
    X : array_like
            First distance matrix (condensed or redundant).
    Y : array_like
            Second distance matrix (condensed or redundant), where the
            order of elements corresponds to the order of elements in
            the first matrix.
    perms : int, optional
            The number of permutations to perform (default: 10000). A
            larger number gives more reliable results but takes longer
            to run. If the number of possible permutations is smaller,
            all permutations will be tested. This can be forced by
            setting perms to 0.
    method : str, optional
            Type of correlation coefficient to use; either 'pearson'
            or 'spearman' (default: 'pearson').
    ignore_nans : bool, optional
            Ignore NaN values in the Y matrix (default: False). This
            can be useful if you have missing values in one of the
            matrices.

    Returns
    -------
    correlations : array of floats
            Computed correlation coefficients between symmetric
            permutations two distance matrices.  The first element of
            the array is the veridical correlation.  The array size is
            the number of computed permutations (see 'perms'
            parameter).
    """

    # Ensure that X and Y are represented as Numpy arrays.
    X = np.asarray(X)
    Y = np.asarray(Y)

    # Check that X and Y are valid distance matrices.
    if (
        spatial.distance.is_valid_dm(np.nan_to_num(X)) == False
        and spatial.distance.is_valid_y(X) == False
    ):
        raise ValueError("X is not a valid condensed or redundant distance matrix")
    if (
        spatial.distance.is_valid_dm(np.nan_to_num(Y)) == False
        and spatial.distance.is_valid_y(Y) == False
    ):
        raise ValueError("Y is not a valid condensed or redundant distance matrix")

    # If X or Y is a redundant distance matrix, reduce it to a condensed distance matrix.
    if len(X.shape) == 2:
        X = spatial.distance.squareform(X, force="tovector", checks=False)
    if len(Y.shape) == 2:
        Y = spatial.distance.squareform(Y, force="tovector", checks=False)

    # Check for size equality.
    if len(X) != len(Y):
        raise ValueError("X and Y are not of equal size")

    # Check for minimum size.
    if len(X) < 3:
        raise ValueError("X and Y should represent at least 3 objects")

    # Check finiteness of X and Y
    if not np.isfinite(X).all():
        raise ValueError(
            "X cannot contain NaNs (but Y may contain NaNs, so consider reordering X and Y)"
        )
    finite_Y = np.isfinite(Y)
    if not ignore_nans and not finite_Y.all():
        raise ValueError('Y may contain NaNs, but "ignore_nans" must be set to True')
    if ignore_nans and finite_Y.all():
        ignore_nans = False  # ignore_nans is True but Y contains no nans

    # If Spearman correlation is requested, convert X and Y to ranks.
    method = method.lower()
    if method == "spearman":
        X, Y = stats.rankdata(X), stats.rankdata(Y)
        Y[~finite_Y] = np.nan  # retain any nans, so that these can be ignored later

    # Check for valid method parameter.
    elif method != "pearson":
        raise ValueError('The method should be set to "pearson" or "spearman"')

    # Now we're ready to start the computation of correlation
    # coefficients of permuted matrices using a number of
    # optimizations:
    #
    # 1. Rather than compute correlation coefficients, we'll just
    #    compute the covariances. This works because the denominator
    #    in the equation for the correlation coefficient will yield
    #    the same result however the objects are permuted, making it
    #    redundant. Removing the denominator leaves us with the
    #    covariance.
    #
    # 2. Rather than permute the Y distances and derive the residuals
    #    to calculate the covariance with the X distances, we'll
    #    represent the Y residuals in the matrix and shuffle those
    #    directly.
    #
    # 3. If the number of possible permutations is less than the
    #    number of permutations that were requested, we'll run a
    #    deterministic test where we try all possible permutations
    #    rather than sample the permutation space. This gives a
    #    faster, deterministic result.

    # Calculate the X and Y residuals, which will be used to compute the
    # covariance under each permutation.
    X_residuals = X - np.mean(X[finite_Y])
    Y_residuals = Y - np.mean(Y[finite_Y])

    # Expand the Y residuals to a redundant matrix.
    Y_residuals_as_matrix = spatial.distance.squareform(
        Y_residuals, force="tomatrix", checks=False
    )

    m = len(Y_residuals_as_matrix)  # number of objects
    n = np.math.factorial(m)  # number of possible matrix permutations

    # If the number of requested permutations is greater than the number of
    # possible permutations (m!) or the perms parameter is set to 0, then run a
    # deterministic Mantel test
    if perms >= n or perms == 0:
        if ignore_nans:
            correlations = deterministic_test_with_nans(m, n, X, Y_residuals_as_matrix)
        else:
            correlations = deterministic_test(m, n, X_residuals, Y_residuals_as_matrix)
        # correlations[0] is the veridical correlation

    else:
        if ignore_nans:
            correlations = stochastic_test_with_nans(m, perms, X, Y_residuals_as_matrix)
        else:
            correlations = stochastic_test(m, perms, X_residuals, Y_residuals_as_matrix)
        correlations[0] = sum(X_residuals[finite_Y] * Y_residuals[finite_Y]) / np.sqrt(
            sum(X_residuals[finite_Y] ** 2) * sum(Y_residuals[finite_Y] ** 2)
        )  # compute veridical correlation and place in positon 0

    return correlations


def mantel_test_from_correlations(correlations, tail="two-tail"):
    """
    Takes correlations previously computed (see compute_correlations()
    function) and return the Mantel test values. The Mantel test is a
    significance test of the correlation between two distance
    matrices.

    Parameters
    ----------
    correlations : array of floats
            The correlations computed by the compute_correlations() function.
    tail : str, optional
            Which tail to test in the calculation of the empirical p-value; either
            'upper', 'lower', or 'two-tail' (default: 'two-tail').

    Returns
    -------
    r : float
            Veridical correlation
    p : float
            Empirical p-value
    m : float
            Arithmetic mean of correlation coefficients computed for permuted matrices
    s : float
            Standard deviation of correlation coefficients computed for permuted matrices

    """
    # Check for valid tail parameter.
    tail = tail.lower()
    if tail not in ["upper", "lower", "two-tail"]:
        raise ValueError('The tail should be set to "upper", "lower", or "two-tail"')

    r = correlations[0]

    if tail == "upper":
        p = sum(correlations >= r) / len(correlations)
    elif tail == "lower":
        p = sum(correlations <= r) / len(correlations)
    elif tail == "two-tail":
        p = sum(abs(correlations) >= abs(r)) / len(correlations)

    m = np.mean(correlations)
    s = np.std(correlations)

    return r, p, m, s


def mantel_test(X, Y, perms=10000, method="pearson", tail="two-tail", ignore_nans=False):
    """
    Takes two distance matrices (either redundant matrices or condensed vectors)
    and performs a Mantel test. The Mantel test is a significance test of the
    correlation between two distance matrices.

    Parameters
    ----------
    X : array_like
            First distance matrix (condensed or redundant).
    Y : array_like
            Second distance matrix (condensed or redundant), where the order of
            elements corresponds to the order of elements in the first matrix.
    perms : int, optional
            The number of permutations to perform (default: 10000). A larger
            number gives more reliable results but takes longer to run. If the
            number of possible permutations is smaller, all permutations will
            be tested. This can be forced by setting perms to 0.
    method : str, optional
            Type of correlation coefficient to use; either 'pearson' or 'spearman'
            (default: 'pearson').
    tail : str, optional
            Which tail to test in the calculation of the empirical p-value; either
            'upper', 'lower', or 'two-tail' (default: 'two-tail').
    ignore_nans : bool, optional
            Ignore NaN values in the Y matrix (default: False). This can be
            useful if you have missing values in one of the matrices.

    Returns
    -------
    r : float
            Veridical correlation
    p : float
            Empirical p-value
    m : float
            Arithmetic mean of correlation coefficients computed for permuted matrices
    s : float
            Standard deviation of correlation coefficients computed for permuted matrices
    """

    correlations = compute_correlations(X, Y, perms=perms, method=method, ignore_nans=ignore_nans)

    return mantel_test_from_correlations(correlations, tail=tail)


def test(X, Y, perms=10000, method="pearson", tail="two-tail", ignore_nans=False):
    """
    Takes two distance matrices (either redundant matrices or condensed vectors)
    and performs a Mantel test. The Mantel test is a significance test of the
    correlation between two distance matrices.

    Parameters
    ----------
    X : array_like
            First distance matrix (condensed or redundant).
    Y : array_like
            Second distance matrix (condensed or redundant), where the order of
            elements corresponds to the order of elements in the first matrix.
    perms : int, optional
            The number of permutations to perform (default: 10000). A larger
            number gives more reliable results but takes longer to run. If the
            number of possible permutations is smaller, all permutations will
            be tested. This can be forced by setting perms to 0.
    method : str, optional
            Type of correlation coefficient to use; either 'pearson' or 'spearman'
            (default: 'pearson').
    tail : str, optional
            Which tail to test in the calculation of the empirical p-value; either
            'upper', 'lower', or 'two-tail' (default: 'two-tail').
    ignore_nans : bool, optional
            Ignore NaN values in the Y matrix (default: False). This can be
            useful if you have missing values in one of the matrices.

    Returns
    -------
    r : float
            Veridical correlation
    p : float
            Empirical p-value
    z : float
            Standard score (z-score)
    """

    r, p, m, s = mantel_test(X, Y, perms=perms, method=method, ignore_nans=ignore_nans)

    z = (r - m) / s

    return r, p, z


def plot(correlations, plot, tail="two-tail",
         significance_level=0.05,
         gaussian_background_color='blue', gaussian_background_alpha=0.1,
         gaussian_color='blue', gaussian_alpha=0.3,
         gaussian_curve_color='blue', gaussian_curve_alpha=0.3,
         hist_fill_color='orange', hist_edge_color='green', hist_alpha=0.7,
         acceptance_color='green', rejection_color='red'):
    """
    Plot the correlations previously computed (see
    compute_correlations() function) and on the given matplotlib axes,
    as well as the theoretical normal distribution with its confidence
    interval highlighted.

    Parameters
    ----------
    correlations : array of floats
            The correlations computed by the compute_correlations() function.
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
    Nothing, but draw on the given plot object.
    Since this function requires the matplotlib library, an exception
    is raised if this package isn't available.
    """
    if not _has_matplotlib:
        raise Exception("In order to produce histograms, you need to install the 'matplotlib' library first.")

    r, p, m, s = mantel_test_from_correlations(correlations, tail=tail)

    tail = tail.lower()
    if tail == "upper":
        z0 = stats.norm.ppf(1e-5)
        z1 = stats.norm.ppf(1 - significance_level)
    elif tail == "lower":
        z0 = stats.norm.ppf(significance_level / 2)
        z1 = stats.norm.ppf(1-1e-5)
    elif tail == "two-tail":
        z0 = stats.norm.ppf(significance_level / 2)
        z1 = stats.norm.ppf(1 - significance_level / 2)

    left  = z0 * s + m
    right = z1 * s + m
    x = np.linspace(left, right, 100)
    y = stats.norm.pdf(x, m, s)

    lower = -5 * s + m
    upper = 5 * s + m
    x_all = np.linspace(lower, upper, 100)
    y_all = stats.norm.pdf(x_all, m, s)

    plot.fill_between(x_all, y_all, 0, color=gaussian_background_color, alpha=gaussian_background_alpha)
    plot.fill_between(x, y, 0, color=gaussian_color, alpha=gaussian_alpha)
    plot.plot(x_all, y_all, color=gaussian_curve_color, alpha=gaussian_curve_alpha)
    plot.vlines(x=[left, right], ymin=0, ymax=[stats.norm.pdf(left, m, s), stats.norm.pdf(right, m, s)], linestyle='-', color=gaussian_curve_color, alpha=gaussian_curve_alpha)

    plot.hist(correlations, bins=20, range=(lower, upper), density=True, color=hist_fill_color, edgecolor=hist_edge_color, alpha=hist_alpha)

    plot.set_xlim(left=lower, right=upper)

    threshold_color = acceptance_color if r <= right else rejection_color
    plot.axvline(x=r, linestyle=':', color=threshold_color)
    plot.annotate('{:.2f}'.format(r), xy=(r, 0.9), color=threshold_color)



def deterministic_test(m, n, X_residuals, Y_residuals_as_matrix):
    Y_residuals_permuted = np.zeros((m ** 2 - m) // 2)
    covariances = np.zeros(n)
    for i, order in enumerate(permutations(range(m))):
        Y_residuals_as_matrix_permuted = Y_residuals_as_matrix[order, :][:, order]
        spatial.distance._distance_wrap.to_vector_from_squareform_wrap(
            Y_residuals_as_matrix_permuted, Y_residuals_permuted
        )
        covariances[i] = (X_residuals * Y_residuals_permuted).sum()
    denominator = np.sqrt(sum(X_residuals ** 2) * sum(Y_residuals_permuted ** 2))
    return covariances / denominator


def deterministic_test_with_nans(m, n, X, Y_residuals_as_matrix):
    Y_residuals_permuted = np.zeros((m ** 2 - m) // 2)
    correlations = np.zeros(n)
    for i, order in enumerate(permutations(range(m))):
        Y_residuals_as_matrix_permuted = Y_residuals_as_matrix[order, :][:, order]
        spatial.distance._distance_wrap.to_vector_from_squareform_wrap(
            Y_residuals_as_matrix_permuted, Y_residuals_permuted
        )
        # Since on each permutation we will be ignoring different values in X,
        # the X_residuals need to be recomputed each time depending on which
        # values in permuted Y are finite.
        finite_Y_permuted = np.isfinite(Y_residuals_permuted)
        reduced_X = X[finite_Y_permuted]
        reduced_X_residuals = reduced_X - reduced_X.mean()
        reduced_Y_residuals = Y_residuals_permuted[finite_Y_permuted]
        covariance = (reduced_X_residuals * reduced_Y_residuals).sum()
        # The denominator will be different on each permutation
        denominator = np.sqrt(
            sum(reduced_X_residuals ** 2) * sum(reduced_Y_residuals ** 2)
        )
        correlations[i] = covariance / denominator
    return correlations


def stochastic_test(m, n, X_residuals, Y_residuals_as_matrix):
    Y_residuals_permuted = np.zeros((m ** 2 - m) // 2)
    covariances = np.zeros(n)
    order = np.arange(m)
    for i in range(1, n):
        np.random.shuffle(order)
        Y_residuals_as_matrix_permuted = Y_residuals_as_matrix[order, :][:, order]
        spatial.distance._distance_wrap.to_vector_from_squareform_wrap(
            Y_residuals_as_matrix_permuted, Y_residuals_permuted
        )
        covariances[i] = (X_residuals * Y_residuals_permuted).sum()
    denominator = np.sqrt(sum(X_residuals ** 2) * sum(Y_residuals_permuted ** 2))
    return covariances / denominator


def stochastic_test_with_nans(m, n, X, Y_residuals_as_matrix):
    Y_residuals_permuted = np.zeros((m ** 2 - m) // 2)
    correlations = np.zeros(n)
    order = np.arange(m)
    for i in range(1, n):
        np.random.shuffle(order)
        Y_residuals_as_matrix_permuted = Y_residuals_as_matrix[order, :][:, order]
        spatial.distance._distance_wrap.to_vector_from_squareform_wrap(
            Y_residuals_as_matrix_permuted, Y_residuals_permuted
        )
        # Since on each permutation we will be ignoring different values in X,
        # the X_residuals need to be recomputed each time depending on which
        # values in permuted Y are finite.
        finite_Y_permuted = np.isfinite(Y_residuals_permuted)
        reduced_X = X[finite_Y_permuted]
        reduced_X_residuals = reduced_X - reduced_X.mean()
        reduced_Y_residuals = Y_residuals_permuted[finite_Y_permuted]
        covariance = (reduced_X_residuals * reduced_Y_residuals).sum()
        # The denominator will be different on each permutation
        denominator = np.sqrt(
            sum(reduced_X_residuals ** 2) * sum(reduced_Y_residuals ** 2)
        )
        correlations[i] = covariance / denominator
    return correlations
