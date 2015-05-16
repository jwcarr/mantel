MantelTest v1.2.0
=================

Efficient Python implementation of the Mantel test (Mantel, 1967). The Mantel test is a significance test of the correlation between two distance matrices.


Description
-----------

This implementation of the Mantel test takes two distance matrices (either redundant matrices or condensed vectors) and returns: the veridical correlation, the empirical p-value, and a standard score (z-score).

Optionally, you can specify the number of permutations to test against (a larger number gives a more reliable p-value and z-score but takes longer to run), which type of correlation coefficient to use (Pearson’s *r*, Spearman’s *ρ*, or Kendall’s *τ*), and which tail to test in the calculation of the empirical p-value.

There are currently two versions of the code: ```Mantel.py``` and ```Mantel_with_Kendall.py```. ```Mantel.py``` is significantly faster (especially when using the Spearman correlation), but it does not support Kendall’s *τ*. ```Mantel_with_Kendall.py``` supports all three correlation methods.


Requirements
------------

- Python 2 or 3 (tested in 2.6, 2.7, and 3.4)
- SciPy/NumPy (any version since 2008 should be fine)


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

We pass the data to the ```Test()``` function of the ```Mantel``` module and optionally specify the number of permutations to test against, a correlation method to use (either ‘pearson’, ‘spearman’, or ‘kendall’), and which tail to test (either ‘upper’ or ‘lower’). In this case, we’ll use the Pearson correlation and test the upper tail, since we’re expecting to find a positive correlation.

```python
Mantel.Test(dists1, dists2, perms=10000, method='pearson', tail='upper')
```

This will measure the veridical Pearson correlation between the two sets of pairwise distances. It then repeatedly measures the correlation again and again under permutations of one of the distance matrices to produce a distribution of correlations under the null hypothesis. Finally, it computes the empirical p-value (the proportion of correlations that were greater than or equal to the veridical correlation) and compares the veridical correlation with the mean and standard deviation of the correlations to generate a z-score.

In this example, the program will return the following:

```python
# r                    p                     z
(0.91489361702127669, 0.041666666666666664, 2.0404024922610229)
```

Since the p-value is less than 0.05 and the z-score is greater than 1.96, we can conclude that there is a significant correlation between these two sets of distances. This suggests that the species that live closer together tend to be more genetically related, while those that live further apart tend to be less genetically related.


Additional notes
----------------

In the example above, we requested 10,000 permutations (the default). However, for four objects there are only 4! = 24 possible permutations. If the number of requested permutations is greater than the actual number of possible permutations for a given matrix size, then the program ignores your request and just tests the veridical against the possible permutations. This gives a deterministic result and can be forced by setting the ```perms``` argument to ```0```. If, however, the number of possible permutations is greater than your request, the program randomly samples the space of possible permutations the requested number of times. This is useful because, in the case of large matrices, it may be intractable to compute all possible permutations. For example, for 13 objects, you’d be looking at multiple days of computation, for 15 objects you’d be looking at multiple years, and 23 objects would take longer than the current age of the universe!


References and links
--------------------

Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. *Cancer Research*, *27*(2), 209–220.

*Mantel Test* on Wikipedia: https://en.wikipedia.org/wiki/Mantel_test

A guide to the Mantel test for linguists: http://www.jonwcarr.net/blog/2014/9/19/a-guide-to-the-mantel-test-for-linguists
