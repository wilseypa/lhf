import numpy as np
import scipy.optimize as opt

matW = np.array([[2,2,-3,1,-1],[1,-3,-2,-1,4]]).T
vecT = np.array([22,8,-6,-2,20])

rowNrmW = np.sqrt(np.sum(matW**2,axis=1))

matA = np.vstack((rowNrmW,matW.T)).T
vecB = vecT

vecC = np.array([-1,0,0])

res = opt.linprog(vecC,A_ub=matA,b_ub = vecB)
print(res)
