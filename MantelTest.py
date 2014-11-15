#!/usr/bin/env python

from scipy import arange, array, corrcoef, mean, random, spatial, stats, std, zeros



# MantelTest()
#   Takes two lists of pairwise distances and performs a Mantel test. Returns
#   the veridical correlation (r), the mean (m) and standard deviation (sd)
#   of the Monte Carlo sample correlations, a Z-score (z) quantifying the
#   significance of the veridical correlation, and a p-value for a normality
#   test on the distribution of sample correlations (norm).

def MantelTest(distances1, distances2, randomizations=10000):
  if ValidateInput(distances1, distances2, randomizations) == False:
    return None
  vector1 = array(distances1, dtype=float)
  vector2 = array(distances2, dtype=float)
  r = corrcoef([vector1, vector2])[0, 1]
  m, sd, norm = MonteCarlo(vector1, vector2, randomizations)
  z = (r-m)/sd
  return r, m, sd, z, norm



# MonteCarlo()
#   Takes two vectors. Measures the correlation between vector 1 and vector 2
#   many times, shuffling vector 2 on each iteration. Returns the mean and
#   standard deviation of the correlations, and a p-value for a normality test
#   of the distribution of correlations.

def MonteCarlo(vector1, vector2, randomizations):
  correlations = zeros(randomizations, dtype=float)
  for i in xrange(0, randomizations):
    correlations[i] = corrcoef([vector1, MatrixShuffle(vector2)])[0, 1]
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
#   Validates input arguments and returns an error message if a problem is
#   identified. Returns True otherwise.

def ValidateInput(distances1, distances2, randomizations):
  if type(randomizations) == int:
    if type(distances1) == list and type(distances2) == list:
      if len(distances1) == len(distances2):
        if spatial.distance.is_valid_y(array(distances1, dtype=float)) == True:
          if spatial.distance.is_valid_y(array(distances2, dtype=float)) == True:
            return True
        print('Error: the sets of pairwise distances are invalid')
        return False
      print('Error: the sets of distances should be of the same length')
      return False
    print('Error: the sets of distances should be lists')
    return False
  print('Error: the number of randomizations should be an integer')
  return False
