MantelTest
==========

Python code for performing a Mantel test (Mantel, 1967). The Mantel test is a statistical test of the correlation between two distance matrices.


Usage
-----

The MantelTest() function takes two sets of pairwise distances and returns: the veridical correlation (r), the mean (m) and standard deviation (sd) of the Monte Carlo sample correlations, and the Z-score (z).

Optionally, you can specify the number of randomizations to perform. A larger number gives a more precise score but takes longer to run.

A Z-score greater that 1.96 (or less than -1.96) indicates a significant correlation at α = 0.05.


Dependencies
------------

This module requires SciPy: http://scipy.org


Example
-------

Let's say we have two sets of four items that we want to correlate. For four items, there are six pariwise distances. For example:

> dists1 = [0.2, 0.4, 0.3, 0.6, 0.9, 0.4]

> dists2 = [0.5, 0.3, 0.3, 0.7, 0.3, 0.6]

We then plug these two sets of pairwise distances into the MantelTest and optionally specify the number of randomizations:

> MantelTest(dists1, dists2, 10000)

This returns the following 4-tuple:

> (-0.090752970081348555, 0.00012251650960981139, 0.44441974239796322, -0.20448120981444243)

In this case, the final number (the Z-score) is not greater than 1.96 or less than -1.96, so we cannot say that there is a significant correlation between the two sets.


Additional notes
----------------

The Mantel test uses a Monte Carlo method to sample the space of possible permutations of one of the distance matrices. This is most useful when the size of your matrix is sufficiently large that it becomes intractable to compute all possible permutations. In practice, this method is best suited to matrices larger than 9×9. Smaller matrices could be computed deterministically.


References and links
--------------------

Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. Cancer Research, 27(2), 209–220.

Wikipedia: https://en.wikipedia.org/wiki/Mantel_test

A guide to the Mantel test for linguists: http://www.jonwcarr.net/blog/2014/9/19/a-guide-to-the-mantel-test-for-linguists
