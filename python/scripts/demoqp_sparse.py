import numpy as np
import scipy.sparse as sp
import qpSWIFT_sparse_bindings as qpSWIFT

# Solver Options
# For information about Solver options please refer to qpSWIFT documentation
opts = {'MAXITER': 30, 'VERBOSE': 1, 'OUTPUT': 2}

# Cost Function
P_dense = np.array([[5.0, 1.0, 0.0],
                    [1.0, 2.0, 1.0],
                    [0.0, 1.0, 4.0]])
P = sp.csc_matrix(P_dense)

c = np.array([1.0, 2.0, 1.0])

# Inequality Constraints
G_dense = np.array([[-4.0, -4.0, 0.0],
                    [0.0, 0.0, -1.0]])
G = sp.csc_matrix(G_dense)

h = np.array([-1.0, -1.0])

# Equality Constraints
A_dense = np.array([[1.0, -2.0, 1.0]])
A = sp.csc_matrix(A_dense)
b = np.array([3.0])

# Equality Constrained QP
reseq = qpSWIFT.run_sparse(c, h, P, G, A, b, opts)

# Inequality Constrained QP
res = qpSWIFT.run_sparse(c, h, P, G, opts=opts)

# Solution
print(res['sol'])
print(reseq['sol'])
