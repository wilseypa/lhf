import collections
import tadasets
import numpy as np
import csv

torus = tadasets.torus(n=2000, c=2, a=1, ambient=200, noise=0.2)
swiss_roll = tadasets.swiss_roll(n=2000, r=4, ambient=10, noise=1.2)
dsphere = tadasets.dsphere(n=1000, d=12, r=3.14, ambient=14, noise=0.14)
inf_sign = tadasets.infty_sign(n=3000, noise=0.1)

np.savetxt('torus.csv', (torus), delimiter=',')
np.savetxt('swiss_roll.csv', (swiss_roll), delimiter=',')
np.savetxt('dsphere.csv', (dsphere), delimiter=',')
np.savetxt('inf_sign.csv', (inf_sign), delimiter=',')
