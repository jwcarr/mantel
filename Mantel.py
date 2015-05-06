#!/usr/bin/env python

from scipy import asarray, random, spatial, stats, zeros

def Test(X, Y, perms=100000, method='pearson'):
  """
  Takes two distance matrices (either redundant matrices or condensed
  vectors) and performs a Mantel test.

  Parameters
  ----------
  X : array_like
      First distance matrix (condensed or redundant)
  Y : array_like
      Second distance matrix (condensed or redundant), where the order
      of elements corresponds to the order of elements in the first matrix
  perms : int, optional
      The number of permutations to perform
  method : str, optional
      Type of correlation coefficient to use; either 'pearson', 'spearman',
      or 'kendall'

  Returns
  -------
  z : float
      A standard score (z-score)
  r : float
      Veridical correlation
  m : float
      Mean of Monte Carlo sample correlations
  sd : float
      Standard deviation of Monte Carlo sample correlations
  """

  X = asarray(X, dtype=float)
  Y = asarray(Y, dtype=float)

  if spatial.distance.is_valid_dm(X) == False and spatial.distance.is_valid_y(X) == False:
    raise ValueError('X is not a valid distance matrix')

  if spatial.distance.is_valid_dm(Y) == False and spatial.distance.is_valid_y(Y) == False:
    raise ValueError('Y is not a valid distance matrix')

  # X is vector and Y is vector
  if len(X.shape) == 1 and len(Y.shape) == 1:
    Y_as_matrix = spatial.distance.squareform(Y, 'tomatrix', False)

  # X is vector and Y is matrix
  elif len(X.shape) == 1 and len(Y.shape) == 2:
    Y_as_matrix = Y
    Y = spatial.distance.squareform(Y, 'tovector', False)

  # X is matrix and Y is vector
  elif len(X.shape) == 2 and len(Y.shape) == 1:
    Y_as_matrix = X
    X, Y = Y, spatial.distance.squareform(X, 'tovector', False)

  # X is matrix and Y is matrix
  elif len(X.shape) == 2 and len(Y.shape) == 2:
    Y_as_matrix = Y
    X = spatial.distance.squareform(X, 'tovector', False)
    Y = spatial.distance.squareform(Y, 'tovector', False)

  else:
    raise ValueError('X and Y should have 1 or 2 dimensions')

  if X.shape[0] != Y.shape[0]:
    raise ValueError('X and Y are not of equal size')

  if method == 'pearson':
    correlate = stats.pearsonr

  elif method == 'spearman':
    correlate = stats.spearmanr

  elif method == 'kendall':
    correlate = stats.kendalltau

  else:
    raise ValueError('The correlation method should be set to "pearson", "spearman", or "kendall"')

  r = correlate(X, Y)[0] # Veridical correlation
  n = Y_as_matrix.shape[0] # Matrix size (N x N)
  MC_corrs = zeros(perms, dtype=float) # Empty array for Monte Carlo correlations

  for i in xrange(perms):
    permutation = random.permutation(n) # Random order in which to permute the matrix
    Y_as_matrix_permuted = Y_as_matrix[permutation, :][:, permutation] # Permute the matrix
    Y_permuted = spatial.distance.squareform(Y_as_matrix_permuted, 'tovector', False) # Convert back to vector
    MC_corrs[i] = correlate(X, Y_permuted)[0] # Store the correlation between X and permuted Y

  m = MC_corrs.mean() # Mean of Monte Carlo correlations
  sd = MC_corrs.std() # Standard deviation of Monte Carlo correlations
  z = (r - m) / sd # Z-score

  return z, r, m, sd
