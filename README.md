MantelTest
==========

Python code for performing a Mantel test (Mantel, 1967). The Mantel test is a statistical test of the correlation between two distance matrices.

Usage
-----

The MantelTest() function takes two distance matrices and returns: the veridical correlation (x), the mean and standard deviation of the Monte Carlo sample (m and s), and the Z-score (z).

Two additional parameters can also be specified:

- *kind* = "matrix" or "vector". Here you specify which format you're supplying: a full redundant distance matrix with 0s on the diagonal, or a condensed distance matrix (i.e. a vector representation of the upper triangle).

- *simulations* = an integer that determines the number of Monte Carlo permutations that will be run. A larger number gives a more precise score. In practice, 1000 should be sufficient.

In general, a Z-score greater that 1.96 indicates a significant correlation between the two distance matrices.

Requirements
------------

SciPy: http://scipy.org

Example 1: Correlation of two redundant distance matrices
---------------------------------------------------------

> dists1 = [[0,2,3],
            [2,0,1],
            [3,1,0]]

> dists2 = [[0,4,3],
            [4,0,2],
            [3,2,0]]

> MantelTest(dists1, dists2, "matrix", 1000)

> (0.5, -0.00165, 0.70411098379442771, 0.71245870543962686)

Example 2: Correlation of two condensed distance matrices
---------------------------------------------------------

> dists1 = [2,3,1]

> dists2 = [4,3,2]

> MantelTest(dists1, dists2, "vector", 1000)

> (0.5, -0.0073000000000000001, 0.7124933052316994, 0.71200669013307916)

Example 1 and Example 2 above are identical datasets which is why they both approximate the same Z-score. Each example is simply a different way of representing the same distance matrices. In this case, the Z-score is low, so we cannot say that there's a significant correlation between the two matrices.

Additional notes
----------------

This implementation uses a Monte Carlo method to sample the space of possible permutations of one of the distance matrices. This is most useful when the size of your matrix is sufficiently large that it becomes intractable to compute all possible permutations. In practice, this method is best suited to matrices larger than 9×9. Smaller matrices could be computed deterministically.

References
----------

Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. Cancer Research, 27(2), 209–220.

Wikipedia: https://en.wikipedia.org/wiki/Mantel_test