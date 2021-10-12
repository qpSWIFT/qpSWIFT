from scipy import sparse
import numpy as np

filePtr = open("some.txt",'w')

n = 4
m = 5
p = 3

P = sparse.random(n, n, density=0.25,format='csc')
A = sparse.random(p, n, density=0.25,format='csc')
G = sparse.random(m, n, density=0.25,format='csc')

s = np.random.rand(1,m)
delta_s = np.random.rand(1,m)
z = np.random.rand(1,m)
delta_z = np.random.rand(1,m)

stat = np.array([n,m,p,P.nnz,A.nnz,G.nnz])

np.savetxt(filePtr,[stat],delimiter=',',fmt='%d')

## P Matrix

np.savetxt(filePtr,[P.data],delimiter=',')
np.savetxt(filePtr,[P.indices],delimiter=',',fmt='%d')
np.savetxt(filePtr,[P.indptr],delimiter=',',fmt='%d')

## A Matrix

np.savetxt(filePtr,[A.data],delimiter=',')
np.savetxt(filePtr,[A.indices],delimiter=',',fmt='%d')
np.savetxt(filePtr,[A.indptr],delimiter=',',fmt='%d')

## G Matrix

np.savetxt(filePtr,[G.data],delimiter=',')
np.savetxt(filePtr,[G.indices],delimiter=',',fmt='%d')
np.savetxt(filePtr,[G.indptr],delimiter=',',fmt='%d')

## axu variables
np.savetxt(filePtr,s,delimiter=',')
np.savetxt(filePtr,delta_s,delimiter=',')
np.savetxt(filePtr,z,delimiter=',')
np.savetxt(filePtr,delta_z,delimiter=',')






filePtr.close()