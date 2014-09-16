from scipy import stats, std, mean
from random import shuffle

def MantelTest(distance_matrix_1, distance_matrix_2, simulations=1000):
    vector_1 = VectorizeMatrix(distance_matrix_1)
    vector_2 = VectorizeMatrix(distance_matrix_2)
    x = stats.pearsonr(vector_1, vector_2)[0]
    m, sd = MonteCarlo(vector_1, distance_matrix_2, simulations)
    z = (x-m)/sd
    return x, m, sd, z

def MonteCarlo(vector_1, distance_matrix_2, simulations):
    correlations = []
    for i in xrange(0, simulations):
        distance_matrix_2_prime = ShuffleMatrix(distance_matrix_2)
        vector_2_prime = VectorizeMatrix(distance_matrix_2_prime)
        correlations.append(stats.pearsonr(vector_1, vector_2_prime)[0])
    return mean(correlations), std(correlations)

def VectorizeMatrix(matrix):
    vector = []
    n = len(matrix)
    for i in range(0, n):
        for j in range(i+1, n):
            vector.append(matrix[i][j])
    return vector

def ShuffleMatrix(matrix):
    n = len(matrix)
    shuffled_matrix = [[] for i in range(0, n)]    
    new_order = range(0, n)
    shuffle(new_order)
    for i in range(0, n):
        for j in range(0, n):
            shuffled_matrix[i].append(matrix[new_order[i]][new_order[j]])
    return shuffled_matrix
