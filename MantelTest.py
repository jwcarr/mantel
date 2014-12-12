#!/usr/bin/env python

from scipy import arange, array, mean, random, spatial, stats, std, zeros



# MantelTest()
#   Takes two lists of pairwise distances and performs a Mantel test. Returns
#   the veridical correlation (r), the mean (m) and standard deviation (sd)
#   of the Monte Carlo sample correlations, a Z-score (z) quantifying the
#   significance of the veridical correlation, and a p-value for a normality
#   test on the distribution of sample correlations (norm).

def MantelTest(distances1, distances2, randomizations=10000, correlation_measure='pearson'):
  ValidateInput(distances1, distances2, randomizations, correlation_measure)
  vector1 = array(distances1, dtype=float)
  vector2 = array(distances2, dtype=float)
  r, p = Correlate(vector1, vector2, correlation_measure)
  m, sd, norm = MonteCarlo(vector1, vector2, randomizations, correlation_measure)
  z = (r-m)/sd
  return r, p, m, sd, z, norm



# Correlate()
#   Correlates two vectors using a given type of correlation coefficient

def Correlate(vector1, vector2, correlation_measure):
  if correlation_measure == 'pearson':
    return stats.pearsonr(vector1, vector2)
  elif correlation_measure == 'spearman':
    return stats.spearmanr(vector1, vector2)
  elif correlation_measure == 'kendall':
    return stats.kendalltau(vector1, vector2)
  else:
    raise ValueError('The correlation_measure should be set to "pearson", "spearman", or "kendall"')



# MonteCarlo()
#   Takes two vectors. Measures the correlation between vector 1 and vector 2
#   many times, shuffling vector 2 on each iteration. Returns the mean and
#   standard deviation of the correlations, and a p-value for a normality test
#   of the distribution of correlations.

def MonteCarlo(vector1, vector2, randomizations, correlation_measure):
  correlations = zeros(randomizations, dtype=float)
  for i in xrange(0, randomizations):
    correlations[i] = Correlate(vector1, MatrixShuffle(vector2), correlation_measure)[0]
  return mean(correlations), std(correlations), stats.normaltest(correlations)[1]



# MatrixShuffle()
#   Takes a vector, converts it to a distance matrix, shuffles the matrix, and
#   returns the upper triangle as a vector.

def MatrixShuffle(vector):
  matrix = spatial.distance.squareform(vector, 'tomatrix')
  shuffled_vector = zeros(vector.shape[0], dtype=float)
  n = matrix.shape[0]
  shuffle_order = arange(n, dtype=int)
  random.shuffle(shuffle_order)
  c = 0
  for i in xrange(0, n-1):
    for j in xrange(i+1, n):
      shuffled_vector[c] = matrix[shuffle_order[i], shuffle_order[j]]
      c += 1
  return shuffled_vector



# ValidateInput()
#   Validates input arguments and raises an error if a problem is identified.

def ValidateInput(distances1, distances2, randomizations, correlation_measure):
  if type(randomizations) != int:
    raise ValueError('The number of randomizations should be an integer')
  if type(correlation_measure) != str or correlation_measure not in ['pearson', 'spearman', 'kendall']:
    raise ValueError('The correlation_measure should be set to "pearson", "spearman", or "kendall"')
  if type(distances1) != list or type(distances2) != list:
    raise ValueError('The sets of pairise distances should be Python lists')
  if len(distances1) != len(distances2):
    raise ValueError('The two sets of pairwise distances should be of the same length')
  if spatial.distance.is_valid_y(array(distances1, dtype=float)) == False:
    raise ValueError('The first set of pairwise distances is invalid')
  if spatial.distance.is_valid_y(array(distances2, dtype=float)) == False:
    raise ValueError('The second set of pairwise distances is invalid')
