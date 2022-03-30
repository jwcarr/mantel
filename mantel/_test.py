from ._mantel import Mantel

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

    mantel_test = Mantel(X, Y, perms, method, ignore_nans)
    return mantel_test.correlations


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

    mantel_test = Mantel(X, Y, permutations=perms, method=method, ignore_nans=ignore_nans, tail=tail)
    r = mantel_test.veridical_correlation
    p = mantel_test.p_value
    z = mantel_test.z_score
    return r, p, z
