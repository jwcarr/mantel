import itertools
import numpy as np
import scipy.spatial as spatial
import scipy.stats as stats
from typing import Iterable
from numpy.typing import ArrayLike

######################################
# Available correlation/tail methods #
######################################

AVAILABLE_CORRELATION_METHODS = ('pearson', 'spearman')

AVAILABLE_TAIL_METHODS = ('lower', 'upper', 'two-tail')


######################################################
# Functions to test the correctness of method values #
######################################################

def _value_in_collection(value: str, collection: Iterable, message: str):
    """
    Ensure that the given string value is in the given collection
    (case matters)
    """
    if value not in collection:
        raise ValueError(message + '\n Available values are: "'
                         + '", "'.join(collection)
                         + '"')


def check_correlation_method(method: str):
    _value_in_collection(method.lower(),
                        AVAILABLE_CORRELATION_METHODS,
                        'Correlation method ("{}") is not available'.format(method))


def check_tail_method(method: str):
    _value_in_collection(method.lower(),
                        AVAILABLE_TAIL_METHODS,
                        'Tail method ("{}") is not available'.format(method))


######################################################
# Functions to compute some statistical informations #
######################################################

def probability_from_sample(value: float, sample: Iterable, tail: str) -> float:
    """
    Return the empirical p-value of the given value according to the
    given sample and the given tail method.

    Parameters
    ----------
    value : float
            The value for which the empirical probability is computed
    sample : array of floats
            The observed values.
    tail : str, optional
            Which tail to test in the calculation of the empirical p-value; either
            'upper', 'lower', or 'two-tail'.

    Returns
    -------
    p : float
            Empirical p-value
    """

    # Check for valid tail parameter.
    tail = tail.lower()
    check_tail_method(tail)

    if tail == "upper":
        p = sum(sample >= value) / len(sample)
    elif tail == "lower":
        p = sum(sample <= value) / len(sample)
    elif tail == "two-tail":
        p = sum(abs(sample) >= abs(value)) / len(sample)

    return p


def confidence_interval(mean: float = 0, std: float = 1,
                        significance_level: float = 0.05,
                        tail: str = 'two-tail') -> (float, float):
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

    # Check for valid significance level
    if not (0 < significance_level < 1):
        raise ValueError("The significance level must be in the range ]0, 1[")
    # Check for valid tail parameter.
    tail = tail.lower()
    check_tail_method(tail)

    if tail == "upper":
        lower_bound = stats.norm.ppf(1e-10, mean, std)
        upper_bound = stats.norm.ppf(1 - significance_level, mean, std)
    elif tail == "lower":
        lower_bound = stats.norm.ppf(significance_level, mean, std)
        upper_bound = stats.norm.ppf(1-1e-10, mean, std)
    elif tail == "two-tail":
        lower_bound = stats.norm.ppf(significance_level / 2, mean, std)
        upper_bound = stats.norm.ppf(1 - significance_level / 2, mean, std)
    return (lower_bound, upper_bound)


###############################################################################
# Correlation computation for distance matrices according to their properties #
###############################################################################

def deterministic_test(m, n, X_residuals, Y_residuals_as_matrix):
    """
    Performs a deterministic Mantel test for input condensed matrices
    of length m * m such that they have non NaN values. This leads to
    the computation of n (= m!) correlations coefficients.
    """
    Y_residuals_permuted = np.zeros((m ** 2 - m) // 2)
    covariances = np.zeros(n)
    for i, order in enumerate(itertools.permutations(range(m))):
        Y_residuals_as_matrix_permuted = Y_residuals_as_matrix[order, :][:, order]
        spatial.distance._distance_wrap.to_vector_from_squareform_wrap(
            Y_residuals_as_matrix_permuted, Y_residuals_permuted
        )
        covariances[i] = (X_residuals * Y_residuals_permuted).sum()
    denominator = np.sqrt(sum(X_residuals ** 2) * sum(Y_residuals_permuted ** 2))
    return covariances / denominator


def deterministic_test_with_nans(m, n, X, Y_residuals_as_matrix):
    """
    Performs a deterministic Mantel test for input condensed matrices
    of length m * m such that they have NaN values. This leads to the
    computation of n (= m!) correlations coefficients.
    """
    Y_residuals_permuted = np.zeros((m ** 2 - m) // 2)
    correlations = np.zeros(n)
    for i, order in enumerate(itertools.permutations(range(m))):
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
    """
    Performs a stochastic Mantel test for input condensed matrices of
    length m * m such that they have no NaN values. This leads to the
    computation of n (< m!) correlations coefficients.
    """
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
    """
    Performs a stochastic Mantel test for input condensed matrices of
    length m * m such that they have NaN values. This leads to the
    computation of n (< m!) correlations coefficients.
    """
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
