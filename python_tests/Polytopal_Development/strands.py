import matplotlib.pyplot as plt
import csv
from mpl_toolkits import mplot3d
from collections import defaultdict
from matplotlib import cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

'''This code is the plotting code for visualizing the strands of the betas mesh as part of constructing polytopes.
Need the vertices.csv, faces.csv, simplices.csv, and the point cloud being used in this directory.
'''

#parses the original point cloud data
def read_pc():
    verts = {}
    with open('swiss_roll.csv', newline='') as pc_file:
        reader = csv.reader(pc_file)
        for idx, row in enumerate(reader):
            verts[idx] = (float(row[0]), float(row[1]), float(row[2]))
        return verts

#parses the face data (edges in our current cases)
def read_edges():
    edges = []
    face_to_verts = {}
    with open('faces.csv', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            face_id = int(row['face_id'])
            v1, v2 = map(int, row['vertex_ids'].split())
            edges.append((v1, v2))
            face_to_verts[face_id] = (v1, v2)
        return edges, face_to_verts

#collect strands and simplices data first, than assign color to each strands.
def read_strands():
    simplices = []
    with open ('simplices.csv', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sid = int(row['simplex_id'])
            f1, f2, f3 = map(int, row['face_ids'].split())
            strand = int(row['strand_id'])
            simplices.append((sid, (f1, f2, f3), strand))
        return simplices
            
verts = read_pc()
edges, face_to_verts = read_edges()
simplices = read_strands()
#building strand adjacency graph for coloring problem
face_to_strands = defaultdict(set)

for _, faces, strand in simplices:
    for f in faces:
        face_to_strands[f].add(strand)

strand_adj = defaultdict(set)
for strands in face_to_strands.values():
    for s in strands:
        strand_adj[s].update(strands - {s})

def generate_colors(n):
    cmap = plt.get_cmap("hsv")
    return [cmap(i / n) for i in range(n)]

palette = generate_colors(200) 

strand_color = {}

for strand in sorted(strand_adj.keys()):
    used = {strand_color[n] for n in strand_adj[strand] if n in strand_color}

    for color in palette:
        if color not in used:
            strand_color[strand] = color
            break
    else:
        raise RuntimeError("Not enough colors to color strands without conflict")

#plot point cloud
ax = plt.axes(projection="3d")

x = [v[0] for v in verts.values()]
y = [v[1] for v in verts.values()]
z = [v[2] for v in verts.values()]
ax.scatter(x, y, z)
#plot faces (edges in 3D case)

for v1, v2 in edges:
    xline = [verts[v1][0], verts[v2][0]]
    yline = [verts[v1][1], verts[v2][1]]
    zline = [verts[v1][2], verts[v2][2]]

    ax.plot(xline, yline, zline, linewidth=0.8)

triangles = []
colors = []

for _, (f1, f2, f3), strand in simplices:
    # Collect all vertex IDs from the three faces
    v1 = face_to_verts[f1]
    v2 = face_to_verts[f2]
    v3 = face_to_verts[f3]

    # Combine and remove duplicates while preserving deterministic order
    vertex_ids = []
    for vid in (*v1, *v2, *v3):
        if vid not in vertex_ids:
            vertex_ids.append(vid)

    if len(vertex_ids) != 3:
        raise ValueError("Invalid simplex, not triangular")

    triangles.append([verts[vid] for vid in vertex_ids])
    colors.append(strand_color[strand])


surf = Poly3DCollection(triangles, alpha=0.75)
surf.set_facecolor(colors)
surf.set_edgecolor("k")
ax.add_collection3d(surf)

plt.tight_layout()
plt.show()

