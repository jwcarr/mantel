from math import isclose
import numpy as np
import mantel


np.random.seed(117)


def test_stochasic_hundred_thousand():
    n_objects = 27
    n_distances = (n_objects**2 - n_objects) // 2
    dists1 = np.random.random(n_distances)
    dists2 = np.random.random(n_distances)
    result = mantel.test(dists1, dists2, perms=100_000)
    assert isclose(result.r, 0, abs_tol=0.1)
    assert len(result.correlations) == 100_000
    assert result.p > 0.05
