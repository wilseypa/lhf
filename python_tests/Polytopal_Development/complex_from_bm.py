import matplotlib.pyplot as plt
import csv
import argparse
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# -----------------------------
# Argument parsing
# -----------------------------
parser = argparse.ArgumentParser(description="Visualize 3D Polytopal Complex")
parser.add_argument('-vh', '--hide-vertices', action='store_true',
                    help="Hide vertex scatter plot")
args = parser.parse_args()


# -----------------------------
# Read vertices (supports ND, we use first 3 dims)
# -----------------------------
def read_vertices(filename):
    verts = {}
    with open(filename, newline='') as f:
        reader = csv.DictReader(f)

        coord_keys = [k for k in reader.fieldnames if k.startswith('x')]

        for row in reader:
            vid = int(row['id'])
            coords = [float(row[k]) for k in coord_keys]

            # Ensure at least 3D (pad if needed)
            while len(coords) < 3:
                coords.append(0.0)

            verts[vid] = coords[:3]

    return verts


# -----------------------------
# Read cells
# -----------------------------
def read_cells(filename):
    cells = []
    with open(filename, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            vertex_ids = list(map(int, row['vertex_indices'].split()))
            cells.append(vertex_ids)
    return cells


# -----------------------------
# Load data
# -----------------------------
verts = read_vertices("pc_vertices.csv")
cells = read_cells("pc_cells.csv")


# -----------------------------
# Plot setup (3D)
# -----------------------------
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


# -----------------------------
# Plot vertices (optional)
# -----------------------------
if not args.hide_vertices:
    x = [v[0] for v in verts.values()]
    y = [v[1] for v in verts.values()]
    z = [v[2] for v in verts.values()]

    ax.scatter(x, y, z, s=20)


# -----------------------------
# Draw cells
# -----------------------------
faces_to_draw = []

for cell in cells:
    pts = [verts[vid] for vid in cell]

    if len(cell) == 2:
        # Edge
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        zs = [p[2] for p in pts]
        ax.plot(xs, ys, zs)

    elif len(cell) >= 3:
        # Face (store for batch drawing)
        faces_to_draw.append(pts)


# -----------------------------
# Draw faces
# -----------------------------
if faces_to_draw:
    poly_collection = Poly3DCollection(faces_to_draw, alpha=0.2)
    ax.add_collection3d(poly_collection)


# -----------------------------
# Formatting
# -----------------------------
ax.set_box_aspect([1, 1, 1])  # equal scaling

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

plt.title("3D Polytopal Complex")
plt.tight_layout()
plt.show()