
import numpy as np
from . import dataGen as dg

def fish(n, d, a = 1):
    fish = '(_x[0]**2 + _x[1]**2)**2 - ((' + str(a) + ') * _x[0] * (_x[0]**2 - _x[1]**2))'
    return np.array(dg.gibbsSampling(n,d,fish,2,0.1,0.5,2.5))
    
    
def infinity(n, d, a = 1):
    infty = '(_x[0]**2 + _x[1]**2)**2 - ((' + str(a) + ') * (_x[0]**2 - _x[1]**2))'
    return np.array(dg.gibbsSampling(n,d,infty,2,0.1,0.5,2.5))

def salmon(n, d, a=0.5, b=0.5):
    salmon = '(_x[0]**2 - ' + str(a) + '**2)**2 + (_x[1]**2 - ' + str(a) + '**2)**2 - ' + str(b) + '**4'
    return np.array(dg.gibbsSampling(n,d,salmon,2,0.1,0.5,2.5))
    
def spiral(n, d, r1=0.3, r2=0.5, a=0.2):
    if(d <= 2):
        print("Dimension must be greater than 2 for spiral object")
        exit(1)
    
    spiral = '(_x[0]-'+str(r2)+'*math.cos(_x[2]/'+str(a)+'))**2 + (_x[1]-'+str(r2) + '*math.sin(_x[2]/'+str(a)+'))**2 - ' + str(r1) + '**2'
    return np.array(dg.gibbsSampling(n,d,spiral,2,0.1,0.5,2.5))

def torus(n, d, r1=0.15, r2=0.7):
    if(d <= 2):
        print("Dimension must be greater than 2 for torus object")
        exit(1)
    
    torus = '((_x[0]**2 + _x[1]**2)**0.5 - ' + str(r2) + ')**2+_x[2]**2-' + str(r1) + '**2'
    return np.array(dg.gibbsSampling(n,d,torus,2,0.1,0.5,2.5))

def ellipsoid(n, d, r=1, a=0.2, b=0.5, c=0.8):
    if(d <= 2):
        print("Dimension must be greater than 2 for ellipsoid object")
        exit(1)
        
    ellipsoid = '(_x[0]**2)/'+str(a)+'**2+(_x[1]**2)/'+str(b)+'**2+(_x[2]**2)/'+str(c)+'**2-'+str(r)+'**2'
    return np.array(dg.gibbsSampling(n,d,ellipsoid,2,0.1,0.5,2.5))

#Commenting out for now, need to confirm further
#def klein(n, d):
#    if(d <= 2):
#        print("Dimension must be greater than 2 for klein object")
#        exit(1)
#    
#    klein = '(_x[0]**2+_x[1]**2+_x[2]**2+2*_x[1]-1)*((_x[0]**2+_x[1]**2+_x[2]**2-2*_x[1]-1)**2-8*_x[2]**2)+16*_x[0]*_x[2]*((_x[0]**2+_x[1]**2+_x[2]**2-2*_x[1]-1))'
#    return np.array(dg.gibbsSampling(n,d,klein,2,0.1,0.5,2.5))
    
