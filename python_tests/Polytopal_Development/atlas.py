import matplotlib.pyplot as plt
import csv
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from collections import defaultdict
import random

'''
Visualizes charts of the beta mesh.
Requires:
    [point_cloud].csv        (point cloud)
    atlas_simplices.csv   (chart output from C++)
'''

# -----------------------------
# Read point cloud
# -----------------------------
def read_pc():
    verts = {}
    with open('torus.csv', newline='') as pc_file:
        reader = csv.reader(pc_file)
        for idx, row in enumerate(reader):
            verts[idx] = (float(row[0]), float(row[1]), float(row[2]))
    return verts


# -----------------------------
# Read atlas simplices
# -----------------------------
def read_charts():
    charts = defaultdict(list)

    with open('atlas.csv', newline='') as f:
        reader = csv.DictReader(f)

        for row in reader:
            chart_id = int(row['chart_id'])
            vertex_ids = list(map(int, row['vertex_ids'].split()))

            charts[chart_id].append(vertex_ids)

    return charts


# -----------------------------
# Generate colors for charts
# -----------------------------
import colorsys

def generate_colors(n):
    colors = []
    golden_ratio = 0.61803398875
    h = 0.0

    for i in range(n):
        h = (h + golden_ratio) % 1
        s = 0.6 + 0.4 * ((i % 3) / 2)   # vary saturation slightly
        v = 0.85                        # keep colors bright
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        colors.append((r, g, b, 1.0))   # RGBA

    return colors


# -----------------------------
# Load data
# -----------------------------
verts = read_pc()
charts = read_charts()

chart_ids = sorted(charts.keys())
colors = generate_colors(len(chart_ids))

chart_color = {cid: colors[i] for i, cid in enumerate(chart_ids)}


# -----------------------------
# Setup 3D plot
# -----------------------------
fig = plt.figure()
ax = plt.axes(projection="3d")

# Plot point cloud
x = [v[0] for v in verts.values()]
y = [v[1] for v in verts.values()]
z = [v[2] for v in verts.values()]

ax.scatter(x, y, z, s=5, alpha=0.5)


# -----------------------------
# Build triangles
# -----------------------------
triangles = []
triangle_colors = []

for cid, simplices in charts.items():

    for simplex in simplices:

        if len(simplex) != 3:
            raise ValueError("Non-triangular simplex detected")

        triangles.append([verts[v] for v in simplex])
        triangle_colors.append(chart_color[cid])


# -----------------------------
# Plot surfaces
# -----------------------------
surf = Poly3DCollection(triangles, alpha=0.75)
surf.set_facecolor(triangle_colors)
surf.set_edgecolor("k")

ax.add_collection3d(surf)

plt.tight_layout()
plt.show()