# MantelTest v1.1.2
# https://github.com/jwcarr/MantelTest
#
# Copyright (c) 2014-2015 Jon W. Carr
#
# The MIT License (MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from numpy import asarray, random, zeros
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
      gives a more reliable Z-score but takes longer to run.
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

  # Ensure X and Y are arrays.

  X = asarray(X, dtype=float)
  Y = asarray(Y, dtype=float)

  # Check that X and Y are valid distance matrices.

  if distance.is_valid_dm(X) == False and distance.is_valid_y(X) == False:
    raise ValueError('X is not a valid distance matrix')

  if distance.is_valid_dm(Y) == False and distance.is_valid_y(Y) == False:
    raise ValueError('Y is not a valid distance matrix')

  # Figure out whether X and Y are matrices or vectors and convert both to
  # vectors and one to a matrix (as needed).

  # X is vector and Y is vector
  if len(X.shape) == 1 and len(Y.shape) == 1:
    Y_as_matrix = distance.squareform(Y, 'tomatrix', False)

  # X is vector and Y is matrix
  elif len(X.shape) == 1 and len(Y.shape) == 2:
    Y_as_matrix = Y
    Y = distance.squareform(Y, 'tovector', False)

  # X is matrix and Y is vector
  elif len(X.shape) == 2 and len(Y.shape) == 1:
    Y_as_matrix = X
    X, Y = Y, distance.squareform(X, 'tovector', False)

  # X is matrix and Y is matrix
  elif len(X.shape) == 2 and len(Y.shape) == 2:
    Y_as_matrix = Y
    X = distance.squareform(X, 'tovector', False)
    Y = distance.squareform(Y, 'tovector', False)

  # Check for size equality.

  if X.shape[0] != Y.shape[0]:
    raise ValueError('X and Y are not of equal size')

  # Assign the relevant correlation function to the variable 'correlate'.

  if method == 'pearson':
    correlate = pearsonr

  elif method == 'spearman':
    correlate = spearmanr

  elif method == 'kendall':
    correlate = kendalltau

  else:
    raise ValueError('The method should be set to "pearson", "spearman", or "kendall"')

  # Determine the size of the matrix (i.e. number of rows/columns).
  n = Y_as_matrix.shape[0]

  # Initialize an empty array to store temporary vector permutations of Y.
  Y_permuted = zeros(Y.shape[0], dtype=float)

  # Initialize an empty array to store the Monte Carlo sample correlations.
  MC_corrs = zeros(perms, dtype=float)

  # Monte Carlo loop.

  for i in range(perms-1):

    # Choose a random order in which to permute the rows/columns of the matrix.
    order = random.permutation(n)

    # Take a permutation of the matrix.
    Y_as_matrix_permuted = Y_as_matrix[order, :][:, order]

    # Condense the permuted version of the matrix into a vector. Rather than use
    # distance.squareform(), we call directly into the C wrapper for speed.
    distance._distance_wrap.to_vector_from_squareform_wrap(Y_as_matrix_permuted, Y_permuted)

    # Compute the correlation coefficient and store it to MC_corrs.
    MC_corrs[i] = correlate(X, Y_permuted)[0]

  # Compute the veridical correlation coefficient.
  r = correlate(X, Y)[0]

  # Include the veridical correlation among the Monte Carlo correlations to
  # prevent the p-value from being 0.
  MC_corrs[perms-1] = r

  # Calculate the empirical p-value for the upper or lower tail.

  if tail == 'upper':
    p = (MC_corrs >= r).sum() / float(perms)

  elif tail == 'lower':
    p = (MC_corrs <= r).sum() / float(perms)

  else:
    raise ValueError('The tail should be set to "upper" or "lower"')

  # Calculate the standard score.

  m = MC_corrs.mean()
  sd = MC_corrs.std()
  z = (r - m) / sd

  return r, p, z
