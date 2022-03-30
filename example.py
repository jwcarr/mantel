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

# Same example of using the object oriented programming

t = mantel.Mantel(dists1, dists2, method="pearson", tail="upper")
r = t.veridical_correlation
p = t.p_value
z = t.z_score
print((r, p, z))

if t.null_hypothesis:
    print("The Mantel test result is to accept the null hypothesis")
    print("=> The distance matrices aren't correlated.")
else:
    print("The Mantel test result is to reject the null hypothesis")
    print("=> The distance matrices are correlated.")
print("   [The probablilty to take the wrong decision is less than {}%.]".format(t.significance_level*100))



# Graphical plot of correlations (requires matplotlib to be installed)

try:
    import matplotlib.pyplot as plt
    import numpy as np
    from mantel._utils import AVAILABLE_TAIL_METHODS

    fig, axs = plt.subplots(2, len(AVAILABLE_TAIL_METHODS))

    plt.suptitle(r"Mantel tests for $\alpha$=5%")
    correlations = mantel.compute_correlations(dists1, dists2, method="pearson")
    for i, tail in enumerate(AVAILABLE_TAIL_METHODS):
        axs[0, i].set_title("Small example with tail set to {}".format(tail))
        mantel.plot_correlations(correlations, axs[0, i], tail=tail)

    # Example using the object oriented programming with random
    # matrices. The Mantel test should lead to accept the null
    # hypothesis most of the time

    # Building the random matrices
    N = 50
    dists1 = np.random.rand(N, N)
    dists2 = np.random.rand(N, N)
    for i in range(N):
        dists1[i][i] = dists2[i][i] = 0
        for j in range(N):
            dists1[i][j] = dists1[j][i]
            dists2[i][j] = dists2[j][i]

    # Building a new instance of the Mantel test for the random matrices.
    t = mantel.Mantel(dists1, dists2, method="pearson")
    for i, tail in enumerate(AVAILABLE_TAIL_METHODS):
        t.plot_correlations(axs[1, i])
        axs[1, i].set_title("Random distance matrices example with tail set to {}".format(tail))

    plt.tight_layout()
    plt.show()

except ImportError:
    print("In order to produce histograms, you need to install the matplotlib library")
