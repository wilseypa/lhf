import numpy as np

'''
This file generates some of the topological structures specifically in 3-dimensions. Being developed for testing and plotting
of the strand identification. Right now just does the klein bottle, I may add more.
'''
def pc_klein_bottle(n, s=1.0, noise=0.0):
    # Sample parameters uniformly
    u = np.random.uniform(0.0, 2.0 * np.pi, n)
    v = np.random.uniform(0.0, 2.0 * np.pi, n)

    # Standard immersed Klein bottle in R^3
    x = np.cos(u) * (np.cos(u / 2.0) * (np.sqrt(2.0) + np.cos(v))
                     + np.sin(u / 2.0) * np.sin(v))
    y = np.sin(u) * (np.cos(u / 2.0) * (np.sqrt(2.0) + np.cos(v))
                     + np.sin(u / 2.0) * np.sin(v))
    z = (-np.sin(u / 2.0) * (np.sqrt(2.0) + np.cos(v))
         + np.cos(u / 2.0) * np.sin(v))

    # Stack and scale
    points = s * np.column_stack((x, y, z))

    # Add noise if requested
    if noise > 0.0:
        points += np.random.normal(scale=noise, size=points.shape)

    # Write CSV (no header)
    np.savetxt("klein_bottle.csv", points, delimiter=",")

def pc_torus(n, s=1.0, noise=0.0, h=1):
    # Sample parameters
    u = np.random.uniform(0.0, 2.0 * np.pi, n)
    v = np.random.uniform(0.0, 2.0 * np.pi, n)

    # Major/minor radii
    R = 2.0
    r = 0.6

    # Frequency-modulated winding to induce h holes
    # (immersed connected sum of tori)
    Ru = R + r * np.cos(v)
    angle = h * u

    x = Ru * np.cos(angle)
    y = Ru * np.sin(angle)
    z = r * np.sin(v)

    # Stack and scale
    points = s * np.column_stack((x, y, z))

    # Add noise if requested
    if noise > 0.0:
        points += np.random.normal(scale=noise, size=points.shape)

    # Write CSV (no header)
    np.savetxt("torus.csv", points, delimiter=",")

pc_torus(10, 1, 0.5)