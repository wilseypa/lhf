import numpy as np
import matplotlib.pyplot as plt

import pypoman.duality as poly
import polytope as pc

from sklearn.decomposition import PCA

fig, ax = plt.subplots()
# creating the dataset
x=np.array([[1,3,0],[2,2,0],[4,2,0],[5,9,0],[2.5,7,0]])
plt.scatter(x[:,0], x[:,1])
x43=np.array([[1,3,0],[2,2,0],[4,2,0],[5,9,0],[2.5,7,0],[1,3,0]])

print(x)
pca = PCA(n_components=3).fit(x)
print(pca.components_)
pca = PCA(n_components=2).fit(x)
pca_x = pca.transform(x)
print(pca.components_)
plt.scatter(pca_x[:,0], pca_x[:,1], alpha=0.3)
#cxfh

x23 = poly.compute_polytope_halfspaces(pca_x)
p23 = pc.Polytope(x23[0], x23[1])
print(p23.chebR) # Chebyshev ball radius
print(p23.chebXc) # Chebyshev ball center
#dxfh
print(type(pca_x))
print(type(p23.chebXc))
pca_x = np.append(pca_x,[p23.chebXc], axis=0)
print(pca_x)
x2 = pca.inverse_transform(pca_x)
plt.scatter(x[:,0], x[:,1], alpha=0.3)
plt.plot(x43[:,0],x43[:,1])
plt.plot(pca_x[:,0],pca_x[:,1])

print(x2)
plt.scatter(x2[:,0], x2[:,1])
circle1 = plt.Circle((2.73301,3.46085), 1.46085185, color='r')
circle2 = plt.Circle((-1.13685,0.18198), 1.46085185, color='r')
ax.add_patch(circle1)
ax.add_patch(circle2)
plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.axis('equal')
plt.show()
print(pca_x)
