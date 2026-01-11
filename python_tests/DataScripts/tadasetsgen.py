import collections
import tadasets
import numpy as np
import pandas as pd
import csv

torus = tadasets.torus(n=50, c=2, a=1, ambient=4, noise=0.2)
swiss_roll = tadasets.swiss_roll(n=50, r=4, ambient=10, noise=1.2)
dsphere = tadasets.dsphere(n=50, d=12, r=3.14, ambient=14, noise=0.14)
inf_sign = tadasets.infty_sign(n=50, noise=0.1)
eyeglasses = tadasets.eyeglasses(n=50, r1=1, r2=2, neck_size=.5, noise=0.1, ambient=4)
theta = np.linspace(0, 8*np.pi, 50)
x = np.cos(theta)
y = np.sin(theta)
z = theta / np.pi
points = np.stack([x, y, z], axis=1)

np.savetxt('torus.csv', (torus), delimiter=',')
np.savetxt('swiss_roll.csv', (swiss_roll), delimiter=',')
np.savetxt('dsphere.csv', (dsphere), delimiter=',')
np.savetxt('inf_sign.csv', (inf_sign), delimiter=',')
np.savetxt('eyeglasses.csv', (eyeglasses), delimiter=',')
pd.DataFrame(points, columns=['x', 'y', 'z']).to_csv('helix.csv', index=False)
