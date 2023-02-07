import numpy as np
from pypoman import compute_polytope_vertices
from pypoman import compute_polytope_halfspaces
import cdd
import pandas as pd
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from scipy.spatial import HalfspaceIntersection
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3
'''
class Polyhedron:
    def __init__(self, halfspaces):
        self.halfspaces = halfspaces

    def get_generators(self):
        generators = []
        for halfspace1 in self.halfspaces:
            is_generator = True
            for halfspace2 in self.halfspaces:
                if all(halfspace1[i] <= halfspace2[i] for i in range(len(halfspace1))):
                    continue
                else:
                    is_generator = False
                    break
            if is_generator:
                generator = [-halfspace1[-1] / coord for coord in halfspace1[:-1]]
                generators.append(generator)
        return generators

        
halfspaces = [[1, 1, 2, -5], [-2, 0, 3, -7], [0, -1, -1, 4]]
p = Polyhedron(halfspaces)
generators = p.get_generators()
print(generators)




vertices = map(
    np.array,
    [[-1, -1], [0, -1.57735], [1,-1], [1,9],[-1,9]],
)


V = np.vstack(vertices)
print(V)
t = np.ones((V.shape[0], 1))  # first column is 1 for vertices
print(t)
tV = np.hstack([t, V])
print(tV)
mat = cdd.Matrix(tV, number_type='float')
print(mat)
mat.rep_type = cdd.RepType.GENERATOR
print(mat)
P = cdd.Polyhedron(mat)
print("Poly " ,P)
bA = np.array(P.get_inequalities())
print(bA)
if bA.shape == (0,):  # bA == []
	print(A)
	print(b)
# the polyhedron is given by b + A x >= 0 where bA = [b|A]
else:
	b, A = np.array(bA[:, 0]), -np.array(bA[:, 1:])
	print(A)
	print(b)
	
	
vertices = map(
    np.array,
    [[-1, -1], [0, -1.57735], [1,-1], [1,9],[-1,9]],
)

A, b = compute_polytope_halfspaces(generators)
print(A)
print(b)


vertices = compute_polytope_vertices(A, b)
print(vertices)
print(A)
print(b)


***************************
'''
import matplotlib.pyplot as plt
hyperplanes = np.loadtxt("input3out.out",delimiter=",")
points = np.loadtxt("inputfile3.txt",delimiter=",")

hyperplanedf = pd.DataFrame(hyperplanes,columns=['X','Y','B'])

A = hyperplanedf[['X','Y']]
B = hyperplanedf[['B']]
for i in range(0,len(B),6):
	a = A.loc[i:(i+5)].to_numpy()
	b = B.loc[i:(i+5)].to_numpy()
	#print(a,b)
	try:
		vertices = compute_polytope_vertices(a, b)
	except Exception:
		pass
	if(len(vertices)>0):
		#print(vertices)
		vertices = np.concatenate(vertices, axis=0 )
		l = int(len(vertices)/2)
		vertices = vertices.reshape(l,2)
		#print(vertices)
		hull = ConvexHull(vertices)
		for simplex in hull.simplices:
			plt.plot(vertices[simplex, 0], vertices[simplex, 1],lw=1)
		#plt.plot(vertices[hull.vertices,0], vertices[hull.vertices,1], 'r--', lw=2)
		#plt.plot(vertices[hull.vertices[0],0], vertices[hull.vertices[0],1], 'ro')
	#plt.show()
	#input()
plt.scatter(points[:,0], points[:,1],s=3)
plt.xlim([-1.5, 1.5])
plt.ylim([-1.5, 1.5])
plt.savefig("2dPartiotiontreeSphere.pdf", bbox_inches = 'tight',pad_inches = 0)
plt.show()

	#input()


'''
class Faces():
    def __init__(self,tri, sig_dig=12, method="convexhull"):
        self.method=method
        self.tri = np.around(np.array(tri), sig_dig)
        self.grpinx = list(range(len(tri)))
        norms = np.around([self.norm(s) for s in self.tri], sig_dig)
        _, self.inv = np.unique(norms,return_inverse=True, axis=0)

    def norm(self,sq):
        cr = np.cross(sq[2]-sq[0],sq[1]-sq[0])
        return np.abs(cr/np.linalg.norm(cr))

    def isneighbor(self, tr1,tr2):
        a = np.concatenate((tr1,tr2), axis=0)
        return len(a) == len(np.unique(a, axis=0))+2

    def order(self, v):
        if len(v) <= 3:
            return v
        v = np.unique(v, axis=0)
        n = self.norm(v[:3])
        y = np.cross(n,v[1]-v[0])
        y = y/np.linalg.norm(y)
        c = np.dot(v, np.c_[v[1]-v[0],y])
        if self.method == "convexhull":
            h = ConvexHull(c)
            return v[h.vertices]
        else:
            mean = np.mean(c,axis=0)
            d = c-mean
            s = np.arctan2(d[:,0], d[:,1])
            return v[np.argsort(s)]

    def simplify(self):
        for i, tri1 in enumerate(self.tri):
            for j,tri2 in enumerate(self.tri):
                if j > i: 
                    if self.isneighbor(tri1,tri2) and \
                       self.inv[i]==self.inv[j]:
                        self.grpinx[j] = self.grpinx[i]
        groups = []
        for i in np.unique(self.grpinx):
            u = self.tri[self.grpinx == i]
            u = np.concatenate([d for d in u])
            u = self.order(u)
            groups.append(u)
        return groups


hyperplanes = np.loadtxt("input3out.out",delimiter=",")
points = np.loadtxt("inputfile3.txt",delimiter=",")

hyperplanedf = pd.DataFrame(hyperplanes,columns=['X','Y','Z','B'])
globalg = []
A = hyperplanedf[['X','Y','Z']]
B = hyperplanedf[['B']]
for i in range(0,len(B),12):
	a = A.loc[i:(i+11)].to_numpy()
	b = B.loc[i:(i+11)].to_numpy()
	#print(a,b)
	try:
		vertices = compute_polytope_vertices(a, b)
	except Exception:
		pass
	if(len(vertices)>0):
		vertices = np.concatenate(vertices, axis=0 )
		l = int(len(vertices)/3)
		vertices = vertices.reshape(l,3)
		hull = ConvexHull(vertices,qhull_options='QJ')
		simplices = hull.simplices
		org_triangles = [vertices[s] for s in simplices]	
		f = Faces(org_triangles)
		g = f.simplify()
		#print(g)
		#print("rohit")
		globalg = globalg + g
		#print(len(globalg))
		#print(globalg)
		#input()

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
colors = list(map("C{}".format, range(len(globalg))))

pc = a3.art3d.Poly3DCollection(globalg,  facecolor=colors, 
                                   edgecolor="k",lw=0.1, alpha=0.1)
ax.add_collection3d(pc)

ax.dist=10
ax.azim=30
ax.elev=10
ax.set_xlim([0,0.5])
ax.set_ylim([0,0.5])
ax.set_zlim([0,0.5])
ax.scatter3D(points[:,0], points[:,1],points[:,2],s=3)
plt.savefig("3dPartiotiontree.pdf", bbox_inches = 'tight',pad_inches = 0)
plt.show()
'''
