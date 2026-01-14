import matplotlib.pyplot as plt
import csv
from mpl_toolkits import mplot3d
from collections import defaultdict
from matplotlib import cm

'''This code is the plotting code for visualizing the strands of the betas mesh as part of constructing polytopes.
Need the vertices.csv, faces.csv, simplices.csv, and the point cloud being used in this directory.
'''

#parses the original point cloud data
def read_pc():
    verts = {}
    with open('torus.csv', newline='') as pc_file:
        reader = csv.reader(pc_file)
        for idx, row in enumerate(reader):
            verts[idx] = (float(row[0]), float(row[1]), float(row[2]))
        return verts

#parses the face data (edges in our current cases)
def read_edges():
    edges = []
    with open('faces.csv', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            v1, v2 = map(int, row['vertex_ids'].split())
            edges.append((v1, v2))
        return edges
'''
def read_simplices():
    simplices = [] #list of faces
    with open ('simplices.csv', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
'''
            
#plot point cloud
ax = plt.axes(projection="3d")
verts = read_pc()
edges = read_edges()
x = [v[0] for v in verts.values()]
y = [v[1] for v in verts.values()]
z = [v[2] for v in verts.values()]
ax.scatter(x, y, z)
#plot faces (edges in 3D case)
'''
for v1, v2 in edges:
    xline = [verts[v1][0], verts[v2][0]]
    yline = [verts[v1][1], verts[v2][1]]
    zline = [verts[v1][2], verts[v2][2]]

    ax.plot(xline, yline, zline, linewidth=0.8)
'''


plt.tight_layout()
plt.show()

