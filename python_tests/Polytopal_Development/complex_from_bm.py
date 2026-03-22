import matplotlib.pyplot as plt
import csv

# -----------------------------
# Read vertices
# -----------------------------
def read_vertices(filename):
    verts = {}
    with open(filename, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            vid = int(row['id'])
            x = float(row['x0'])
            y = float(row['x1'])
            verts[vid] = (x, y)
    return verts


# -----------------------------
# Read cells (edges)
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
# Plot setup
# -----------------------------
fig, ax = plt.subplots()

# Plot vertices
x = [v[0] for v in verts.values()]
y = [v[1] for v in verts.values()]
ax.scatter(x, y, s=20)

# -----------------------------
# Draw edges (cells)
# -----------------------------
for cell in cells:
    if len(cell) != 2:
        continue  # skip non-edges for now

    v0, v1 = cell
    x_vals = [verts[v0][0], verts[v1][0]]
    y_vals = [verts[v0][1], verts[v1][1]]

    ax.plot(x_vals, y_vals)


# -----------------------------
# Formatting
# -----------------------------
ax.set_aspect('equal')
plt.title("Polytopal Complex (Edges)")
plt.tight_layout()
plt.show()