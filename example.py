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
