MantelTest
==========

Efficient Python code for performing a Mantel test (Mantel, 1967). The Mantel test is a significance test of the correlation between two distance matrices.


Description
-----------

This implementation of the Mantel test takes two distance matrices (either redundant matrices or condensed vectors) and returns: the veridical correlation, the mean and standard deviation of the Monte Carlo sample correlations, and a Z-score quantifying the significance of the veridical correlation.

Optionally, you can specify the number of randomizations to perform (a larger number gives a more precise z-score but takes longer to run) and which type of correlation coefficient to use (Pearson’s *r*, Spearman’s *ρ*, or Kendall’s *τ*).

A Z-score greater that 1.96 (or less than -1.96) indicates a significant correlation at *α* = 0.05.

There are currently two versions: ```Mantel.py``` and ```Mantel_with_Kendall.py```. ```Mantel.py``` is significantly more efficient but does not support Kendall’s *τ*.


Requirements
------------

SciPy: http://scipy.org


Usage
-----

First, let’s import the module:

```python
import Mantel
```

Let’s say we have a set of four items and we want to correlate (i) the distances between the four items using one measure with (ii) the corresponding distances between the four items using another measure. For example, your “items” might be species of animal, and your two measures might be genetic distance and geographical distance (the hypothesis being that species that live far away from each other will tend to be more genetically different).

For four items, there are six pairwise distances. First you should compute the pairwise distances for each measure and store the distances in two lists or arrays (i.e. condensed distance vectors). Alternatively, you can compute the full redundant distance matrices; this program will accept either format. No distance functions are included in this module, since the metrics you use will be specific to your particular data.

Let’s say our data looks like this:

```python
#         E.g. species A through D
#         A~B  A~C  A~D  B~C  B~D  C~D
dists1 = [0.2, 0.4, 0.3, 0.6, 0.9, 0.4] # E.g. genetic distances
dists2 = [0.5, 0.3, 0.3, 0.7, 0.3, 0.6] # E.g. geographical distances
```

We input the data into the Mantel test and optionally specify the number of randomizations to perform and a correlation method to use (either ‘pearson’, ‘spearman’, or ‘kendall’):

```python
Mantel.Test(dists1, dists2, 10000, 'pearson')
```

This will measure the veridical Pearson correlation between the two lists of pairwise distances. It then proceeds to repeatedly measure the correlation again and again under Monte Carlo permutations of one of the distance matrices. Finally, it compares the veridical correlation with the mean and standard deviation of the Monte Carlo sample correlations to generate a Z-score.

In this example, the program would return the following 4-tuple:

```python
# z          r         m         sd
(-0.205..., -0.091..., 0.001..., 0.446...)
```

Since the Z-score is not greater than 1.96 (nor less than -1.96), we cannot conclude that there is a significant correlation between the two sets of distances.


Additional notes
----------------

The Mantel test uses a Monte Carlo method to sample the space of possible permutations of one of the distance matrices. This is most useful when the size of your distance matrix is sufficiently large that it is intractable to compute all possible permutations. In practice, this method is best suited to cases where you have more than 9 items. For 9 items or fewer, it is possible to try all possible permutations in a reasonable amount of time (< 1 minute). For 13 items, you’d be looking at multiple days of computation to try all possible permutations, for 15 items you’d be looking at multiple years, and 23 items would take longer than the current age of the universe!


References and links
--------------------

Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. *Cancer Research*, *27*(2), 209–220.

*Mantel Test* on Wikipedia: https://en.wikipedia.org/wiki/Mantel_test

A guide to the Mantel test for linguists: http://www.jonwcarr.net/blog/2014/9/19/a-guide-to-the-mantel-test-for-linguists
