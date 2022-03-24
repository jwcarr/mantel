from ._plot import plot_correlations as _plot_correlations
from ._utils import \
    check_correlation_method as _check_correlation_method, \
    check_tail_method as _check_tail_method, \
    probability_from_sample as _probability_from_sample, \
    confidence_interval as _confidence_interval, \
    deterministic_test as _deterministic_test, \
    deterministic_test_with_nans as _deterministic_test_with_nans, \
    stochastic_test as _stochastic_test, \
    stochastic_test_with_nans as _stochastic_test_with_nans
import numpy as np
import scipy.spatial as spatial
import scipy.stats as stats
from typing import Iterable, Tuple
try:
    from numpy.typing import ArrayLike
except:
    from typing import Type
    ArrayLike = Type[typing.TypeVar("array-like")]


class Mantel:
    """
    Constructs an object that performs a Mantel test.

    The Mantel test is a significance test of the correlation between
    two distance matrices, where the null hypothesis states that the
    two matrices aren't correlated. If the veridical correlation
    between the two matrices are outside from the confidence interval,
    then the null hypothesis should be rejected in flavor of the
    alternative hypothesis (the two matrices are correlated).
    """

    def __init__(self, X: ArrayLike, Y: ArrayLike,
                 permutations: int = 10000,
                 method: str = "pearson",
                 ignore_nans: bool = False,
                 tail : str = 'two-tail',
                 significance_level: float = 0.05):
        """
        Build a Mantel instance that allows to compare two distance
        matrices (either redundant matrices or condensed vectors) and
        computes the correlation between symmetric permutations two
        distance matrices.

        Parameters
        ----------
        X : array_like
            First distance matrix (condensed or redundant).
        Y : array_like
            Second distance matrix (condensed or redundant), where the
            order of elements corresponds to the order of elements in
            the first matrix.
        permutations : int, optional
            The number of permutations to perform (default: 10000). A
            larger number gives more reliable results but takes longer
            to run. If the number of possible permutations is smaller,
            all permutations will be tested. This can be forced by
            setting permutations to 0.
        method : str, optional
            Type of correlation coefficient to use; either 'pearson'
            or 'spearman' (default: 'pearson').
        ignore_nans : bool, optional
            Ignore NaN values in the Y matrix (default: False). This
            can be useful if you have missing values in one of the
            matrices.
        tail : str, optional
            Which tail to test in the calculation of the empirical p-value; either
            'upper', 'lower', or 'two-tail' (default: 'two-tail').
        significance_level : float
            Significance level of the null hypothesis of the Mantel
            test (default: 5%)
        """
        self._X = np.asarray(X)
        self._Y = np.asarray(Y)
        self._permutations = permutations
        _check_correlation_method(method)
        self._method = method.lower()
        self._ignore_nans = ignore_nans
        self.tail_method = tail
        self.significance_level = significance_level

        # Check that X and Y are valid distance matrices.
        if (
                spatial.distance.is_valid_dm(np.nan_to_num(self._X)) == False
                and spatial.distance.is_valid_y(self._X) == False
        ):
            raise ValueError("X is not a valid condensed or redundant distance matrix")
        if (
                spatial.distance.is_valid_dm(np.nan_to_num(self._Y)) == False
                and spatial.distance.is_valid_y(self._Y) == False
        ):
            raise ValueError("Y is not a valid condensed or redundant distance matrix")

        # If X or Y is a redundant distance matrix, reduce it to a condensed distance matrix.
        if len(self._X.shape) == 2:
            self._X = spatial.distance.squareform(self._X, force="tovector", checks=False)
        if len(self._Y.shape) == 2:
            self._Y = spatial.distance.squareform(self._Y, force="tovector", checks=False)

        # Check for size equality.
        if len(self._X) != len(self._Y):
            raise ValueError("X and Y are not of equal size")

        # Check for minimum size.
        if len(self._X) < 3:
            raise ValueError("X and Y should represent at least 3 objects")

        # Check finiteness of X and Y
        if not np.isfinite(self._X).all():
            raise ValueError(
                "X cannot contain NaNs (but Y may contain NaNs, so consider reordering X and Y)"
            )

        self._finite_Y = np.isfinite(self._Y)
        if not self._ignore_nans and not self._finite_Y.all():
            raise ValueError('Y may contain NaNs, but "ignore_nans" must be set to True')
        if self._ignore_nans and self._finite_Y.all():
            self._ignore_nans = False  # ignore_nans is True but Y contains no nans

        # If Spearman correlation is requested, convert X and Y to ranks.
        if self._method == "spearman":
            self._X, self._Y = stats.rankdata(self._X), stats.rankdata(self._Y)
            self._Y[~self._finite_Y] = np.nan  # retain any nans, so that these can be ignored later

        self._compute_correlations()


    @property
    def X(self) -> ArrayLike:
        """
        The first distance matrix in condensed form. If the
        correlation in use is the correlation of Spearman (see
        correlation_method) then it corresponds to the ranks of the
        distances.
        """
        return self._X


    @property
    def Y(self) -> ArrayLike:
        """
        The second distance matrix in condensed form. If the
        correlation in use is the correlation of Spearman (see
        correlation_method) then it corresponds to the ranks of the
        distances.
        """
        return self._Y


    @property
    def permutations(self) -> int:
        """
        The number of permutations used. If it was greater than the
        number possible permutations or if the permutations parameter
        was set to 0, then this attribute is the number of possible
        permutations.
        """
        return self._permutations


    @property
    def correlation_method(self) -> str:
        """
        The correlation method used to compare the input
        matrices. Currently, only 'pearson' and 'spearman' values are
        available.
        """
        return self._method


    @property
    def ignore_nans(self) -> bool:
        """
        If the input matrices contains NaN values and if the Mantel
        test was created setting this parameter to True, then it is
        True, otherwise it is False (even if current instance was
        built with this parameter set to True).
        """
        return self._ignore_nans


    @property
    def tail_method(self) -> str:
        """
        The way the null hypothesis of the Mantel test is tested. The
        veridical correlation coefficient can be tested against the
        'lower' tail of the distribution, the 'upper' tail of the
        distribution or the 'two-tail' of the distribution.
        """
        return self._tail


    @tail_method.setter
    def tail_method(self, tail: str) -> None:
        _check_tail_method(tail)
        self._tail = tail


    @property
    def significance_level(self) -> float:
        return self._significance_level


    @significance_level.setter
    def significance_level(self, significance_level: float) -> None:
        if not (0 < significance_level < 1):
            raise ValueError("The significance level must be in the range ]0, 1[")
        self._significance_level = significance_level


    @property
    def correlations(self) -> Iterable[float]:
        """
        The array containing the correlation coefficients computed by
        the Mantel test. The first element (at index 0) is the
        correlation coefficient of the input distance matrices (see
        veridical_correlation attribute).
        """
        return self._correlations


    @property
    def veridical_correlation(self) -> float:
        """
        The correlation coefficient of the input distance matrices.
        """
        return self._correlations[0]


    @property
    def mean(self) -> float:
        """
        The mean of the correlations computed by the current Mantel
        test.
        """
        return self._mean


    @property
    def std(self) -> float:
        """
        The standard deviation of the correlations computed by the
        current Mantel test.
        """
        return self._std


    @property
    def p_value(self) -> float:
        """
        The current Mantel test empirical p-value.
        """
        return self._p_value


    @property
    def z_score(self) -> float:
        """
        The current Mantel test z-score.
        """
        return (self.veridical_correlation - self.mean) / self.std


    @property
    def confidence_interval(self) -> Tuple[float, float]:
        """
        The confidence interval of the current Mantel test.
        """
        return _confidence_interval(mean = self.mean,
                                    std = self.std,
                                    significance_level = self.significance_level,
                                    tail = self.tail_method)


    @property
    def null_hypothesis(self) -> bool:
        """
        True if the current Mantel test result is to accept the null
        hypothesis according to the given significance level (and
        false otherwise).
        """
        l, u = _confidence_interval(significance_level = self.significance_level,
                                    tail = self.tail_method)
        return l <= self.z_score <= u


    def _compute_correlations(self) -> None:
        """
        Computes the correlations (and set the corresponding
        attribute) between symmetric permutations of the two distance
        matrices.

        This method set the mean, std and p_value (and thus the
        z_score) properties of the current instance.
        """

        # We're ready to start the computation of correlation
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
        X_residuals = self._X - np.mean(self._X[self._finite_Y])
        Y_residuals = self._Y - np.mean(self._Y[self._finite_Y])

        # Expand the Y residuals to a redundant matrix.
        Y_residuals_as_matrix = spatial.distance.squareform(
            Y_residuals, force="tomatrix", checks=False
        )

        m = len(Y_residuals_as_matrix)  # number of objects
        n = np.math.factorial(m)  # number of possible matrix permutations

        # If the number of requested permutations is greater than the number of
        # possible permutations (m!) or the permutations parameter is set to 0, then run a
        # deterministic Mantel test
        if self._permutations >= n or self._permutations == 0:
            self._permutations = n
            if self._ignore_nans:
                correlations = _deterministic_test_with_nans(m, n, self._X, Y_residuals_as_matrix)
            else:
                correlations = _deterministic_test(m, n, X_residuals, Y_residuals_as_matrix)
            # correlations[0] is the veridical correlation
        else:
            if self._ignore_nans:
                correlations = _stochastic_test_with_nans(m, self._permutations, self._X, Y_residuals_as_matrix)
            else:
                correlations = _stochastic_test(m, self._permutations, X_residuals, Y_residuals_as_matrix)
            correlations[0] = sum(X_residuals[self._finite_Y] * Y_residuals[self._finite_Y]) / np.sqrt(
                sum(X_residuals[self._finite_Y] ** 2) * sum(Y_residuals[self._finite_Y] ** 2)
            )  # compute veridical correlation and place in position 0

        self._correlations = correlations
        self._mean = np.mean(correlations)
        self._std = np.std(correlations)
        self._p_value = _probability_from_sample(value = self.veridical_correlation,
                                                 sample = self.correlations,
                                                 tail = self.tail_method)

    def plot_correlations(self, plot,
                          gaussian_background_color='blue', gaussian_background_alpha=0.1,
                          gaussian_color='blue', gaussian_alpha=0.3,
                          gaussian_curve_color='blue', gaussian_curve_alpha=0.3,
                          hist_fill_color='orange', hist_edge_color='green', hist_alpha=0.7,
                          acceptance_color='green', rejection_color='red') -> None:
        """
        Plot the correlations previously computed (see
        compute_correlations() function) and on the given matplotlib
        axes, as well as the theoretical normal distribution with its
        confidence interval highlighted.

        Parameters
        ----------
        plot : matplotlib.axes.Axes
            The matplotlib figure to draw on.
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
        _plot_correlations(self.correlations, plot,
                           significance_level = self.significance_level,
                           gaussian_background_color = gaussian_background_color,
                           gaussian_background_alpha = gaussian_background_alpha,
                           gaussian_color = gaussian_color,
                           gaussian_alpha = gaussian_alpha,
                           gaussian_curve_color = gaussian_curve_color,
                           gaussian_curve_alpha = gaussian_curve_alpha,
                           hist_fill_color = hist_fill_color,
                           hist_edge_color = hist_edge_color,
                           hist_alpha = hist_alpha,
                           acceptance_color = acceptance_color,
                           rejection_color = rejection_color)
