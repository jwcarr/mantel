import numpy as np
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


# Example with random data (implying no correlation)

n_objects = 27
n_distances = (n_objects**2 - n_objects) // 2

dists1 = np.random.random(n_distances)
dists2 = np.random.random(n_distances)

result = mantel.test(dists1, dists2)

print(result.r)
print(result.p)
print(result.z)

print(result.p < 0.05)

print(result.correlations)
print(result.mean)
print(result.std)


# Plotting example (requires matplotlib)

fig, axis = mantel.plot(result)
fig.savefig('example.svg')
