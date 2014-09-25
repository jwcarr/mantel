from scipy import mean, random, spatial, stats, std

#   Takes two distance matrices and performs a Mantel test. Returns a
#   Z-score quantifying the strength of the correlation between the two
#   distance matrices.

def MantelTest(distances1, distances2, kind="matrix", simulations=1000):
    if kind == "matrix" or kind == "m":
        vector1 = spatial.distance.squareform(distances1, "tovector")
        vector2 = spatial.distance.squareform(distances2, "tovector")
        matrix2 = distances2
    elif kind == "vector" or kind == "v":
        vector1 = distances1
        vector2 = distances2
        matrix2 = spatial.distance.squareform(distances2, "tomatrix")
    else:
        print "Error: The parameter 'kind' should be set to 'matrix' or 'vector'."
        return None
    x = stats.pearsonr(vector1, vector2)[0]
    m, sd = MonteCarlo(vector1, matrix2, simulations)
    return (x-m)/sd

#   Takes a vector and matrix. Shuffles the matrix some number of times and
#   measures the correlation with the vector for each shuffle. Returns the
#   mean and standard deviation of the correlations.

def MonteCarlo(vector, matrix, simulations):
    correlations = []
    for i in xrange(0, simulations):
        vector2 = spatial.distance.squareform(ShuffleMatrix(matrix), "tovector")
        correlations.append(stats.pearsonr(vector, vector2)[0])
    return mean(correlations), std(correlations)

#   Shuffles the rows and columns of a matrix, maintaining the order of
#   elements along the columns and down the rows.

def ShuffleMatrix(matrix):
    n = len(matrix)
    shuffled_matrix = []    
    order = range(0, n)
    random.shuffle(order)
    for i in xrange(0, n):
        row = []
        for j in xrange(0, n):
            row.append(matrix[order[i]][order[j]])
        shuffled_matrix.append(row)
    return shuffled_matrix
