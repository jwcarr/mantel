# MantelTest v1.2.5
# http://jwcarr.github.io/MantelTest/
#
# Copyright (c) 2014-2015 Jon W. Carr
# Licensed under the terms of the MIT License

from itertools import permutations
from math import factorial
from numpy import arange, asarray, random, zeros
from scipy.spatial import distance
from scipy.stats import pearsonr, spearmanr, kendalltau

def Test(X, Y, perms=10000, method='pearson', tail='upper'):
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
      The number of permutations to perform (default: 10000). A larger number
      gives more reliable results but takes longer to run. If the actual number
      of possible permutations is smaller, the program will enumerate all
      permutations. Enumeration can be forced by setting this argument to 0.
  method : str, optional
      Type of correlation coefficient to use; either 'pearson', 'spearman', or
      'kendall' (default: 'pearson'). N.B. the time complexity of Kendall's tau
      scales exponentially with matrix size, so it is slow for large matrices.
  tail : str, optional
      Which tail to test in the calculation of the empirical p-value; either
      'upper' or 'lower' (default: 'upper').

  Returns
  -------
  r : float
      Veridical correlation
  p : float
      Empirical p-value
  z : float
      Standard score (z-score)
  """

  # Ensure that X and Y are formatted as Numpy arrays.

  X = asarray(X, dtype=float)
  Y = asarray(Y, dtype=float)

  # Check that X and Y are valid distance matrices/vectors.

  if distance.is_valid_dm(X) == False and distance.is_valid_y(X) == False:
    raise ValueError('X is not a valid (condensed) distance matrix')

  if distance.is_valid_dm(Y) == False and distance.is_valid_y(Y) == False:
    raise ValueError('Y is not a valid (condensed) distance matrix')

  # Figure out whether X and Y are matrices or vectors and convert both to
  # vectors and one to a matrix (as needed).

  # X is vector and Y is vector
  if len(X.shape) == 1 and len(Y.shape) == 1:
    Y_as_matrix = distance.squareform(Y, force='tomatrix', checks=False)

  # X is vector and Y is matrix
  elif len(X.shape) == 1 and len(Y.shape) == 2:
    Y_as_matrix = Y
    Y = distance.squareform(Y, force='tovector', checks=False)

  # X is matrix and Y is vector
  elif len(X.shape) == 2 and len(Y.shape) == 1:
    Y_as_matrix = X
    X, Y = Y, distance.squareform(X, force='tovector', checks=False)

  # X is matrix and Y is matrix
  elif len(X.shape) == 2 and len(Y.shape) == 2:
    Y_as_matrix = Y
    X = distance.squareform(X, force='tovector', checks=False)
    Y = distance.squareform(Y, force='tovector', checks=False)

  # Check for size equality.

  if X.shape[0] != Y.shape[0]:
    raise ValueError('X and Y are not of equal size')

  # Check for minimum size.

  if X.shape[0] < 3:
    raise ValueError('X and Y should represent at least 3 objects')

  # Assign the relevant correlation function to the variable 'correlate'.

  if method == 'pearson':
    correlate = pearsonr

  elif method == 'spearman':
    correlate = spearmanr

  elif method == 'kendall':
    correlate = kendalltau

  else:
    raise ValueError('The method should be set to "pearson", "spearman", or "kendall"')

  # Check for valid tail parameter.

  if tail != 'upper' and tail != 'lower':
    raise ValueError('The tail should be set to "upper" or "lower"')

  m = Y_as_matrix.shape[0] # Number of objects
  n = factorial(m) # Number of matrix permutations

  # Initialize an empty array to store temporary vector permutations of Y.
  Y_permuted = zeros(Y.shape[0], dtype=float)

  # If the number of requested permutations is greater than the number of
  # possible permutations (m!) or the perms parameter is set to 0, then run a
  # deterministic Mantel test ...

  if perms >= n or perms == 0:

    # Initialize an empty array to store the correlations.
    corrs = zeros(n, dtype=float)

    # Enumerate all permutations of row/column orders.
    orders = permutations(range(m))

    perms = 0

    for order in orders:

      # Take a permutation of the matrix.
      Y_as_matrix_permuted = Y_as_matrix[order, :][:, order]

      # Condense the permuted version of the matrix. Rather than use
      # distance.squareform(), we call directly into the C wrapper for speed.
      distance._distance_wrap.to_vector_from_squareform_wrap(Y_as_matrix_permuted, Y_permuted)

      # Compute the correlation coefficient and store it to corrs.
      corrs[perms] = correlate(X, Y_permuted)[0]

      perms += 1

  # ... otherwise run a stochastic Mantel test.

  else:

    # Initialize an empty array to store the correlations.
    corrs = zeros(perms, dtype=float)

    # Initialize an array to store the permutation order.
    order = arange(m)

    # Store the veridical correlation coefficient in first position.
    corrs[0] = correlate(X, Y)[0]

    for i in range(1, perms):

      # Choose a random order in which to permute the rows and columns.
      random.shuffle(order)

      # Take a permutation of the matrix.
      Y_as_matrix_permuted = Y_as_matrix[order, :][:, order]

      # Condense the permuted version of the matrix. Rather than use
      # distance.squareform(), we call directly into the C wrapper for speed.
      distance._distance_wrap.to_vector_from_squareform_wrap(Y_as_matrix_permuted, Y_permuted)

      # Compute the correlation coefficient and store it to corrs.
      corrs[i] = correlate(X, Y_permuted)[0]

  # Assign veridical correlation to r.
  r = corrs[0]

  # Calculate the empirical p-value for the upper or lower tail.

  if tail == 'upper':
    p = (corrs >= r).sum() / float(perms)

  elif tail == 'lower':
    p = (corrs <= r).sum() / float(perms)

  # Calculate the standard score.

  m = corrs.mean()
  sd = corrs.std()
  z = (r - m) / sd

  return r, p, z
