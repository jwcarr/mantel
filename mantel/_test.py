from itertools import permutations
import numpy as np
from scipy import spatial, stats


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

    # Check for valid tail parameter.
    tail = tail.lower()
    if tail not in ["upper", "lower", "two-tail"]:
        raise ValueError('The tail should be set to "upper", "lower", or "two-tail"')

    # Now we're ready to start the Mantel test using a number of optimizations:
    #
    # 1. We don't need to recalculate the pairwise distances between the objects
    #    on every permutation. They've already been calculated, so we can use a
    #    simple matrix shuffling technique to avoid recomputing them. This works
    #    like memoization.
    #
    # 2. Rather than compute correlation coefficients, we'll just compute the
    #    covariances. This works because the denominator in the equation for the
    #    correlation coefficient will yield the same result however the objects
    #    are permuted, making it redundant. Removing the denominator leaves us
    #    with the covariance.
    #
    # 3. Rather than permute the Y distances and derive the residuals to calculate
    #    the covariance with the X distances, we'll represent the Y residuals in
    #    the matrix and shuffle those directly.
    #
    # 4. If the number of possible permutations is less than the number of
    #    permutations that were requested, we'll run a deterministic test where
    #    we try all possible permutations rather than sample the permutation
    #    space. This gives a faster, deterministic result.

    # Calculate the X and Y residuals, which will be used to compute the
    # covariance under each permutation.
    X_residuals = X - np.mean(X[finite_Y])
    Y_residuals = Y - np.mean(Y[finite_Y])

    # Expand the Y residuals to a redundant matrix.
    Y_residuals_as_matrix = spatial.distance.squareform(
        Y_residuals, force="tomatrix", checks=False
    )

    # Get the number of objects.
    m = len(Y_residuals_as_matrix)

    # Calculate the number of possible matrix permutations.
    n = np.math.factorial(m)

    # Initialize an empty array to store temporary permutations of Y_residuals.
    Y_residuals_permuted = np.zeros(len(Y_residuals))

    # If the number of requested permutations is greater than the number of
    # possible permutations (m!) or the perms parameter is set to 0, then run a
    # deterministic Mantel test ...
    if perms >= n or perms == 0:

        # Initialize an empty array to store the covariances.
        covariances = np.zeros(n)

        # Iterate over all permutations of row/column orders.
        for i, order in enumerate(permutations(range(m))):

            # Take a permutation of the matrix.
            Y_residuals_as_matrix_permuted = Y_residuals_as_matrix[order, :][:, order]

            # Condense the permuted version of the matrix. Rather than use
            # distance.squareform(), we call directly into the C wrapper for speed.
            spatial.distance._distance_wrap.to_vector_from_squareform_wrap(
                Y_residuals_as_matrix_permuted, Y_residuals_permuted
            )

            # Compute and store the covariance.
            if ignore_nans:
                finite_Y_permuted = np.isfinite(Y_residuals_permuted)
                # since we are now ignoring different values in X, the X
                # residuals need to be recomputed because mean is different
                reduced_X = X[finite_Y_permuted]
                reduced_X_residuals = reduced_X - reduced_X.mean()
                reduced_Y_residuals = Y_residuals_permuted[finite_Y_permuted]
                covariances[i] = (reduced_X_residuals * reduced_Y_residuals).sum()
            else:
                covariances[i] = (X_residuals * Y_residuals_permuted).sum()

    # ... otherwise run a stochastic Mantel test.
    else:

        # Initialize an empty array to store the covariances.
        covariances = np.zeros(perms)

        # Initialize an array to store the permutation order.
        order = np.arange(m)

        # Store the veridical covariance in 0th position...
        covariances[0] = sum(X_residuals[finite_Y] * Y_residuals[finite_Y])

        # ...and then run the random permutations.
        for i in range(1, perms):

            # Choose a random order in which to permute the rows and columns.
            np.random.shuffle(order)

            # Take a permutation of the matrix.
            Y_residuals_as_matrix_permuted = Y_residuals_as_matrix[order, :][:, order]

            # Condense the permuted version of the matrix. Rather than use
            # distance.squareform(), we call directly into the C wrapper for speed.
            spatial.distance._distance_wrap.to_vector_from_squareform_wrap(
                Y_residuals_as_matrix_permuted, Y_residuals_permuted
            )

            # Compute and store the covariance.
            if ignore_nans:
                finite_Y_permuted = np.isfinite(Y_residuals_permuted)
                # since we are now ignoring different values in X, the X
                # residuals need to be recomputed because mean is different
                reduced_X = X[finite_Y_permuted]
                reduced_X_residuals = reduced_X - reduced_X.mean()
                reduced_Y_residuals = Y_residuals_permuted[finite_Y_permuted]
                covariances[i] = (reduced_X_residuals * reduced_Y_residuals).sum()
            else:
                covariances[i] = (X_residuals * Y_residuals_permuted).sum()

    # Calculate the veridical correlation coefficient from the veridical covariance.
    r = covariances[0] / np.sqrt(
        sum(X_residuals[finite_Y] ** 2) * sum(Y_residuals[finite_Y] ** 2)
    )

    # Calculate the empirical p-value for the upper or lower tail.
    if tail == "upper":
        p = sum(covariances >= covariances[0]) / len(covariances)
    elif tail == "lower":
        p = sum(covariances <= covariances[0]) / len(covariances)
    elif tail == "two-tail":
        p = sum(abs(covariances) >= abs(covariances[0])) / len(covariances)

    # Calculate the standard score.
    z = (covariances[0] - np.mean(covariances)) / np.std(covariances)

    return r, p, z
