#import numpy as np
#import scipy.optimize as opt

#matW = np.array([[-1,2],[4,-3]]).T
#vecT = np.array([20,8])

#rowNrmW = np.sqrt(np.sum(matW**2,axis=1))

#matA = np.vstack((rowNrmW,matW.T)).T
#vecB = vecT

#vecC = np.array([0,-1,0])

#res = opt.linprog(vecC,A_ub=matA,b_ub = vecB)
#print(res)
import numpy as np
import polytope as pc
import pypoman.duality as poly
A = np.array([[2, 1],
              [2, -3],
              [-1, 4],
              [-3, -2],
              [-1,1]])

b = np.array([22, 8, 20, -6,0])

p = pc.Polytope(A, b)
#p = pc.box2poly([[0, 2], [0, 1]])

p.dim # number of dimensions of ambient Euclidean space
p.volume # measure in ambient space
print(p.chebR) # Chebyshev ball radius
print(p.chebXc) # Chebyshev ball center
p.cheby
p.bounding_box

x = poly.compute_polytope_halfspaces([[7.556,6.889],[6.667,6.667],[1.2,1.2],[2.615,-0.923],[9.25,3.5]])
print(x)
