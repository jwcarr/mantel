MantelTest v1.1.2
=================

Efficient Python implementation of the Mantel test (Mantel, 1967). The Mantel test is a significance test of the correlation between two distance matrices.


Description
-----------

This implementation of the Mantel test takes two distance matrices (either redundant matrices or condensed vectors) and returns: the veridical correlation, the empirical p-value, and a standard score (z-score).

Optionally, you can specify the number of randomizations to perform (a larger number gives a more reliable p-value and z-score but takes longer to run), which type of correlation coefficient to use (Pearson’s *r*, Spearman’s *ρ*, or Kendall’s *τ*), and which tail to test in the calculation of the empirical p-value.

There are currently two versions of the code: ```Mantel.py``` and ```Mantel_with_Kendall.py```. ```Mantel.py``` is significantly more efficient but does not support Kendall’s *τ*.


Requirements
------------

SciPy: http://scipy.org


Usage example
-------------

First, let’s import the module:

```python
import Mantel
```

Let’s say we have a set of four objects and we want to correlate X (the distances between the four objects using one measure) with Y (the corresponding distances between the four objects using another measure). For example, your “objects” might be species of animal, and your two measures might be genetic distance and geographical distance (the hypothesis being that species that live far away from each other will tend to be more genetically different).

For four objects, there are six pairwise distances. First you should compute the pairwise distances for each measure and store the distances in two lists or arrays (i.e. condensed distance vectors). Alternatively, you can compute the full redundant distance matrices; this program will accept either format. No distance functions are included in this module, since the metrics you use will be specific to your particular data.

Let’s say our data looks like this:

```python
#         E.g. species A through D
#         A~B  A~C  A~D  B~C  B~D  C~D
dists1 = [0.2, 0.4, 0.3, 0.6, 0.9, 0.4] # E.g. genetic distances
dists2 = [0.3, 0.3, 0.2, 0.7, 0.8, 0.3] # E.g. geographical distances
```

We pass the data to the ```Test()``` function of the ```Mantel``` module and optionally specify the number of randomizations to perform, a correlation method to use (either ‘pearson’, ‘spearman’, or ‘kendall’), and which tail to test (either ‘upper’ or ‘lower’). In this case, we’ll use the Pearson correlation and test the upper tail, since we’re expecting to find a positive correlation.

```python
Mantel.Test(dists1, dists2, 10000, 'pearson', 'upper')
```

This will measure the veridical Pearson correlation between the two sets of pairwise distances. It then repeatedly measures the correlation again and again under Monte Carlo permutations of one of the distance matrices to produce a distribution of correlations under the null hypothesis. Finally, it computes the empirical p-value (the proportion of sample correlations that were greater than or equal to the veridical correlation) and compares the veridical correlation with the mean and standard deviation of the sample correlations to generate a z-score.

In this example, the program will return something like the following:

```python
# r                    p                     z
(0.91489361702127669, 0.041621999999999999, 2.0409128395130831)
```

Since the p-value is less than 0.05 and the z-score is greater than 1.96, we can conclude that there is a significant correlation between these two sets of distances. This suggests that the species that live closer together tend to be more genetically related, while those that live further apart tend to be less genetically related.


Additional notes
----------------

The Mantel test uses a Monte Carlo method to sample the space of possible permutations of one of the distance matrices. This is most useful when the size of your distance matrix is sufficiently large that it is intractable to compute all possible permutations. In practice, this method is best suited to cases where you have more than 9 objects. For 9 objects or fewer, it is possible to try all possible permutations in a reasonable amount of time (< 1 minute). For 13 objects, you’d be looking at multiple days of computation to try all possible permutations, for 15 objects you’d be looking at multiple years, and 23 objects would take longer than the current age of the universe!


References and links
--------------------

Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. *Cancer Research*, *27*(2), 209–220.

*Mantel Test* on Wikipedia: https://en.wikipedia.org/wiki/Mantel_test

A guide to the Mantel test for linguists: http://www.jonwcarr.net/blog/2014/9/19/a-guide-to-the-mantel-test-for-linguists
