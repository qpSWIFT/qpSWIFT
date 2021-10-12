import numpy as np 
import qpSWIFT

import tracemalloc


from pympler.tracker import SummaryTracker
tracker = SummaryTracker()

###  Solver Options
### For information about Solver options please refer to qpSWIFT
### documentation

tracemalloc.start()

opts = {"MAXITER":30,"VERBOSE":0}

### Cost Function

P = np.array([[5.0,1.0,0.0],
         [1.0, 2.0, 1.0],
         [0.0, 1.0, 4.0]])
    
c = np.array([1.0,2.0,1.0])

### Inequality Constraints
G = np.array([[-4.0,-4.0,0.0],
        [0.0,0.0,-1.0]])

h = np.array([-1.0,-1.0])

### Equality Constraints
A = np.array([[1.0, -2.0, 1.0]])

b = np.array([3.0])

### Equality Constrained QP
# snapshot1 = tracemalloc.take_snapshot()
reseq = qpSWIFT.run(c,h,P,G,A,b,opts)
# snapshot2 = tracemalloc.take_snapshot()

### Inequality Constrained QP
res = qpSWIFT.run(c,h,P,G,opts=opts)

tracker.print_diff()
### Solution
# print(res['sol'])

# print(reseq['sol'])`

# snapshot = tracemalloc.take_snapshot()
# top_stats = snapshot.statistics('lineno')

# top_stats = snapshot2.compare_to(snapshot1, 'lineno')
# print("[ Top 10 ]")
# for stat in top_stats[:10]:
#     print(stat)
