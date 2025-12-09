from math import isclose
import numpy as np
import mantel


condensed_X = [0.2, 0.4, 0.3, 0.6, 0.9, 0.4]
condensed_Y = [0.3, 0.3, 0.2, 0.7, 0.8, 0.3]
condensed_Y_with_nan = [0.3, 0.3, 0.2, 0.7, 0.8, np.nan]

redundant_X = [
    [0.0, 0.2, 0.4, 0.3],
    [0.2, 0.0, 0.6, 0.9],
    [0.4, 0.6, 0.0, 0.4],
    [0.3, 0.9, 0.4, 0.0],
]
redundant_Y = [
    [0.0, 0.3, 0.3, 0.2],
    [0.3, 0.0, 0.7, 0.8],
    [0.3, 0.7, 0.0, 0.3],
    [0.2, 0.8, 0.3, 0.0],
]


def test_condensed_pearson_upper():
    result = mantel.test(condensed_X, condensed_Y, method="pearson", tail="upper")
    assert isclose(result.r, 0.9148936170212766)
    assert isclose(result.p, 0.041666666666666664)
    result = mantel.test(condensed_Y, condensed_X, method="pearson", tail="upper")
    assert isclose(result.r, 0.9148936170212766)
    assert isclose(result.p, 0.041666666666666664)


def test_redundant_pearson_upper():
    result = mantel.test(redundant_X, redundant_Y, method="pearson", tail="upper")
    assert isclose(result.r, 0.9148936170212766)
    assert isclose(result.p, 0.041666666666666664)
    result = mantel.test(redundant_Y, redundant_X, method="pearson", tail="upper")
    assert isclose(result.r, 0.9148936170212766)
    assert isclose(result.p, 0.041666666666666664)


def test_condensed_spearman_upper():
    result = mantel.test(condensed_X, condensed_Y, method="spearman", tail="upper")
    assert isclose(result.r, 0.8316554899054915)
    assert isclose(result.p, 0.041666666666666664)
    result = mantel.test(condensed_Y, condensed_X, method="spearman", tail="upper")
    assert isclose(result.r, 0.8316554899054915)
    assert isclose(result.p, 0.041666666666666664)


def test_redundant_spearman_upper():
    result = mantel.test(redundant_X, redundant_Y, method="spearman", tail="upper")
    assert isclose(result.r, 0.8316554899054915)
    assert isclose(result.p, 0.041666666666666664)
    result = mantel.test(redundant_Y, redundant_X, method="spearman", tail="upper")
    assert isclose(result.r, 0.8316554899054915)
    assert isclose(result.p, 0.041666666666666664)


def test_mixed_pearson_upper():
    result = mantel.test(condensed_X, redundant_Y, method="pearson", tail="upper")
    assert isclose(result.r, 0.9148936170212766)
    assert isclose(result.p, 0.041666666666666664)
    result = mantel.test(condensed_Y, redundant_X, method="pearson", tail="upper")
    assert isclose(result.r, 0.9148936170212766)
    assert isclose(result.p, 0.041666666666666664)


def test_mixed_spearman_upper():
    result = mantel.test(condensed_X, redundant_Y, method="spearman", tail="upper")
    assert isclose(result.r, 0.8316554899054915)
    assert isclose(result.p, 0.041666666666666664)
    result = mantel.test(condensed_Y, redundant_X, method="spearman", tail="upper")
    assert isclose(result.r, 0.8316554899054915)
    assert isclose(result.p, 0.041666666666666664)


def test_condensed_pearson_upper_with_nan():
    result = mantel.test(
        condensed_X,
        condensed_Y_with_nan,
        method="pearson",
        tail="upper",
        ignore_nans=True,
    )
    assert isclose(result.r, 0.920327285673818)
    assert isclose(result.p, 0.041666666666666664)


def test_condensed_spearman_upper_with_nan():
    result = mantel.test(
        condensed_X,
        condensed_Y_with_nan,
        method="spearman",
        tail="upper",
        ignore_nans=True,
    )
    assert isclose(result.r, 0.8459061860169812)
    assert isclose(result.p, 0.125)
