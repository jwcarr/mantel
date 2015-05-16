# MantelTest v1.2.1
# http://jwcarr.github.io/MantelTest/
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

from itertools import permutations
from math import factorial
from numpy import asarray, random, sqrt, zeros
from scipy.spatial import distance
from scipy.stats import rankdata

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
      of possilbe permutations is smaller, the program will enumerate all
      permutations. Enumeration can be forced by setting this argument to 0.
  method : str, optional
      Type of correlation coefficient to use; either 'pearson' or 'spearman'
      (default: 'pearson').
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

  # Check that X and Y are valid distance matrices/vectors.

  if distance.is_valid_dm(X) == False and distance.is_valid_y(X) == False:
    raise ValueError('X is not a valid distance matrix')

  if distance.is_valid_dm(Y) == False and distance.is_valid_y(Y) == False:
    raise ValueError('Y is not a valid distance matrix')

  # If X or Y is a matrix, condense it to a vector.

  if len(X.shape) == 2:
    X = distance.squareform(X, force='tovector', checks=False)

  if len(Y.shape) == 2:
    Y = distance.squareform(Y, force='tovector', checks=False)

  # Check for size equality.

  if X.shape[0] != Y.shape[0]:
    raise ValueError('X and Y are not of equal size')

  # Check for minimum size.

  if X.shape[0] < 3:
    raise ValueError('X and Y should represent at least 3 objects')

  # If Spearman correlation is requested, convert X and Y to ranks.

  if method == 'spearman':
    X = rankdata(X)
    Y = rankdata(Y)

  elif method != 'pearson':
    raise ValueError('The method should be set to "pearson" or "spearman"')

  # Most parts of the correlation coefficient will be the same for every
  # permutation and can therefore be computed outside the loop.

  X_res = X - X.mean() # X residuals
  Y_res = Y - Y.mean() # Y residuals
  X_ss = (X_res * X_res).sum() # X sum-of-squares
  Y_ss = (Y_res * Y_res).sum() # Y sum-of-squares
  denominator = sqrt(X_ss * Y_ss) # Denominator of the correlation coefficient

  # Although Y_res will be the same set of numbers on every permutation, the
  # order will be different each time. Therefore, we reformat Y_res as a matrix
  # so that we can take matrix permutations of the Y residuals.
  Y_res_as_matrix = distance.squareform(Y_res, force='tomatrix', checks=False)

  # Determine the size of the matrix (i.e. number of rows/columns).
  n = Y_res_as_matrix.shape[0]

  # Initialize an empty array to store temporary vector permutations of Y_res.
  Y_res_permuted = zeros(Y_res.shape[0], dtype=float)

  # Either enumerate all permutations ...

  if perms >= factorial(n) or perms == 0:

    # Initialize an empty array to store the correlations.
    corrs = zeros(factorial(n), dtype=float)

    # Enumerate all permutations of row/column orders.
    orders = permutations(range(n))

    perms = 0

    for order in orders:

      # Take a permutation of the matrix.
      Y_res_as_matrix_permuted = Y_res_as_matrix[order, :][:, order]

      # Condense the permuted version of the matrix. Rather than use
      # distance.squareform(), we call directly into the C wrapper for speed.
      distance._distance_wrap.to_vector_from_squareform_wrap(Y_res_as_matrix_permuted, Y_res_permuted)

      # Compute the correlation coefficient and store it to corrs.
      corrs[perms] = (X_res * Y_res_permuted).sum() / denominator

      perms += 1

  # ... or randomly sample from the space of permutations.

  else:

    # Initialize an empty array to store the correlations.
    corrs = zeros(perms, dtype=float)

    # Store the veridical correlation coefficient first.
    corrs[0] = (X_res * Y_res).sum() / denominator

    for i in range(1, perms):

      # Choose a random order in which to permute the rows and columns.
      order = random.permutation(n)

      # Take a permutation of the matrix.
      Y_res_as_matrix_permuted = Y_res_as_matrix[order, :][:, order]

      # Condense the permuted version of the matrix. Rather than use
      # distance.squareform(), we call directly into the C wrapper for speed.
      distance._distance_wrap.to_vector_from_squareform_wrap(Y_res_as_matrix_permuted, Y_res_permuted)

      # Compute the correlation coefficient and store it to corrs.
      corrs[i] = (X_res * Y_res_permuted).sum() / denominator

  # Assign veridical correlation to r.
  r = corrs[0]

  # Calculate the empirical p-value for the upper or lower tail.

  if tail == 'upper':
    p = (corrs >= r).sum() / float(perms)

  elif tail == 'lower':
    p = (corrs <= r).sum() / float(perms)

  else:
    raise ValueError('The tail should be set to "upper" or "lower"')

  # Calculate the standard score.

  m = corrs.mean()
  sd = corrs.std()
  z = (r - m) / sd

  return r, p, z
