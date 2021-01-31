import mantel

dists1 = [0.2, 0.4, 0.3, 0.6, 0.9, 0.4] # E.g. genetic distances
dists2 = [0.3, 0.3, 0.2, 0.7, 0.8, 0.3] # E.g. geographical distances

result = mantel.test(dists1, dists2, perms=10000, method='pearson', tail='upper')

print(result)
