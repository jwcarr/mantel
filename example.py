import mantel

# Example with condendensed distance matrices

dists1 = [0.2, 0.4, 0.3, 0.6, 0.9, 0.4]
dists2 = [0.3, 0.3, 0.2, 0.7, 0.8, 0.3]

result = mantel.test(dists1, dists2, method="pearson", tail="upper")
print(result)

# Example with redundant distance matrices

dists1 = [[0.0,0.2,0.4,0.3],
          [0.2,0.0,0.6,0.9],
          [0.4,0.6,0.0,0.4],
          [0.3,0.9,0.4,0.0]]
dists2 = [[0.0,0.3,0.3,0.2],
          [0.3,0.0,0.7,0.8],
          [0.3,0.7,0.0,0.3],
          [0.2,0.8,0.3,0.0]]

result = mantel.test(dists1, dists2, method="pearson", tail="upper")
print(result)

try:
    import matplotlib.pyplot as plt
    import numpy as np

    tails = ['lower', 'upper', 'two-tail']
    fig, axs = plt.subplots(2, len(tails))

    plt.suptitle(r"Mantel tests for $\alpha$=5%")
    correlations = mantel.compute_correlations(dists1, dists2, method="pearson")
    result = mantel.mantel_test_from_correlations(correlations)
    print("{} => {}".format(result, (result[0] - result[2]) /result[3]))
    for i, tail in enumerate(tails):
        axs[0, i].set_title("Small example with tail set to {}".format(tail))
        mantel.plot(correlations, axs[0, i], tail=tail)

    N = 50
    dists1 = np.random.rand(N, N)
    dists2 = np.random.rand(N, N)
    for i in range(N):
        dists1[i][i] = dists2[i][i] = 0
        for j in range(N):
            dists1[i][j] = dists1[j][i]
            dists2[i][j] = dists2[j][i]
    correlations = mantel.compute_correlations(dists1, dists2, method="pearson")
    for i, tail in enumerate(tails):
        mantel.plot(correlations, axs[1, i], tail=tail)
        mantel.plot(correlations, axs[1, i], tail=tail)
        axs[1, i].set_title("Random distance matrices example with tail set to {}".format(tail))

    plt.tight_layout()
    plt.show()

except ImportError:
    print("In order to produce histograms, we recommand you installing the matplotlib package")
