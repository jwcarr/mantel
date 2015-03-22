MantelTest
==========

Python code for performing a Mantel test (Mantel, 1967). The Mantel test is a significance test for the correlation between two distance matrices.


Description
-----------

The MantelTest() function takes two lists of pairwise distances and returns: the veridical correlation and p-value, the mean and standard deviation of the Monte Carlo sample correlations, a Z-score quantifying the significance of the veridical correlation, and a p-value for a normality test of the distribution of sample correlations.

Optionally, you can specify the number of randomizations to perform (a larger number gives a more precise z-score but takes longer to run) and which type of correlation coefficient to use (Pearson's r, Spearman's r, or Kendall's τ).

A Z-score greater that 1.96 (or less than -1.96) indicates a significant correlation at α = 0.05.


Requirements
------------

SciPy: http://scipy.org


Usage
-----

First, let’s import the module:

```python
import MantelTest
```

Let's say we have a set of four items and we want to correlate the distances between the four items using one measure with the corresponding distances between the four items using another measure. For example, your “items” might be species of animal, and your two measures might be genetic distance and geographical distance (the hypothesis being that species that live far away from each other will tend to be more genetically different).

For four items, there are six pariwise distances. First we derive the pairwise distances for each measure, which can be represented in two vectors (i.e. condensed distance matrices). In Python, these can simply be lists. For example:

```python
dists1 = [0.2, 0.4, 0.3, 0.6, 0.9, 0.4]
dists2 = [0.5, 0.3, 0.3, 0.7, 0.3, 0.6]
```

We plug these two lists of pairwise distances into the MantelTest and optionally specify the number of randomizations to perform and a correlation method:

```python
MantelTest.MantelTest(dists1, dists2, 10000, 'pearson')
```

This will measure the veridical correlation between the two lists of pairwise distances. Then we proceed to repeatedly measure the correlation again and again under random permutations of one of the distance matrices. Finally, we compare our veridical correlation with the mean and standard deviation of the Monte Carlo sample correlations to generate a Z-score.

In this example, the program would return the following 6-tuple:

```python
(-0.090752..., 0.864244..., 0.000122..., 0.444419..., -0.204481..., 0.0)
```

Since the fifth number (the Z-score) is not greater than 1.96 (nor less than -1.96), we cannot say that there is a significant correlation between the two sets of distances.


Additional notes
----------------

The Mantel test uses a Monte Carlo method to sample the space of possible permutations of one of the distance matrices. This is most useful when the size of your matrix is sufficiently large that it is intractable to compute all possible permutations. In practice, this method is best suited to matrices larger than 9×9 (i.e. you have more than 9 items). Smaller matrices could be computed deterministically.

For the use of a Z-score to be valid, the sample correlations should be normally distributed. This program runs a D’Agostino normality test and returns its p-value (the final return value). If this is not less than 0.05, then the Z-score should not be relied on to make inferences about the correlation between the distance matrices.


References and links
--------------------

Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. *Cancer Research*, *27*(2), 209–220.

*Mantel Test* on Wikipedia: https://en.wikipedia.org/wiki/Mantel_test

A guide to the Mantel test for linguists: http://www.jonwcarr.net/blog/2014/9/19/a-guide-to-the-mantel-test-for-linguists
