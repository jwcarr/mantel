# MantelTest v1.2.3
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
      of possible permutations is smaller, the program will enumerate all
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

  # Ensure that X and Y are formatted as Numpy arrays.

  X = asarray(X, dtype=float)
  Y = asarray(Y, dtype=float)

  # Check that X and Y are valid distance matrices/vectors.

  if distance.is_valid_dm(X) == False and distance.is_valid_y(X) == False:
    raise ValueError('X is not a valid (condensed) distance matrix')

  if distance.is_valid_dm(Y) == False and distance.is_valid_y(Y) == False:
    raise ValueError('Y is not a valid (condensed) distance matrix')

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

  # Check for valid method parameter.

  elif method != 'pearson':
    raise ValueError('The method should be set to "pearson" or "spearman"')

  # Check for valid tail parameter.

  if tail != 'upper' and tail != 'lower':
    raise ValueError('The tail should be set to "upper" or "lower"')

  # Calculate the X and Y residuals, which will be used to compute the
  # covarience under each permutation.

  X_residuals = X - X.mean()
  Y_residuals = Y - Y.mean()

  # Let's assume that we're going to permute the Y objects. Although the
  # Y_residuals will be the same set of numbers on every permutation, their
  # order will be different each time. Here we reformat Y_residuals as a matix
  # in order to take matrix permutations. Matrix permuting the residuals is a
  # shortcut that avoids the computation of the entire correlation coefficient
  # on every permutation.

  Y_res_as_matrix = distance.squareform(Y_residuals, force='tomatrix', checks=False)

  # Determine the size of the matrix (i.e. number of rows/columns/objects).
  n = Y_res_as_matrix.shape[0]

  # Initialize an empty array to store temporary permutations of Y_residuals.
  Y_res_permuted = zeros(Y_residuals.shape[0], dtype=float)

  # If the number of requested permutations is greater than the number of
  # possible permutations (n!) or the perms parameter is set to 0, then run a
  # deterministic Mantel test ...

  if perms >= factorial(n) or perms == 0:

    # Initialize an empty array to store the covariences.
    covariences = zeros(factorial(n), dtype=float)

    # Enumerate all permutations of row/column orders.
    orders = permutations(range(n))

    perms = 0

    for order in orders:

      # Take a permutation of the matrix.
      Y_res_as_matrix_permuted = Y_res_as_matrix[order, :][:, order]

      # Condense the permuted version of the matrix. Rather than use
      # distance.squareform(), we call directly into the C wrapper for speed.
      distance._distance_wrap.to_vector_from_squareform_wrap(Y_res_as_matrix_permuted, Y_res_permuted)

      # Compute and store the covarience.
      covariences[perms] = (X_residuals * Y_res_permuted).sum()

      perms += 1

  # ... otherwise run a stochastic Mantel test.

  else:

    # Initialize an empty array to store the covariences.
    covariences = zeros(perms, dtype=float)

    # Store the veridical covarience in first position.
    covariences[0] = (X_residuals * Y_residuals).sum()

    for i in range(1, perms):

      # Choose a random order in which to permute the rows and columns.
      order = random.permutation(n)

      # Take a permutation of the matrix.
      Y_res_as_matrix_permuted = Y_res_as_matrix[order, :][:, order]

      # Condense the permuted version of the matrix. Rather than use
      # distance.squareform(), we call directly into the C wrapper for speed.
      distance._distance_wrap.to_vector_from_squareform_wrap(Y_res_as_matrix_permuted, Y_res_permuted)

      # Compute and store the covarience.
      covariences[i] = (X_residuals * Y_res_permuted).sum()

  # Calculate the veridical correlation coefficient.
  r = covariences[0] / sqrt((X_residuals * X_residuals).sum() * (Y_residuals * Y_residuals).sum())

  # Calculate the empirical p-value for the upper or lower tail.

  if tail == 'upper':
    p = (covariences >= covariences[0]).sum() / float(perms)

  elif tail == 'lower':
    p = (covariences <= covariences[0]).sum() / float(perms)

  # Calculate the standard score.
  z = (covariences[0] - covariences.mean()) / covariences.std()

  return r, p, z
