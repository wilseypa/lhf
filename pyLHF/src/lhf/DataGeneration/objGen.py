
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
    
def spiral2d(n, d, a=1):
    spiral = '_x[1]/_x[0] - math.tan(' + str(a) + '/(math.sqrt(_x[0]**2 + _x[1]**2)))'
    return np.array(dg.gibbsSampling(n,d,spiral,2,0.1,0.5,2.5))
    
def spiral(n, d, r1=0.3, r2=0.5, a=0.2):
    if(d <= 2):
        print("Dimension must be greater than 2 for spiral object; calling spiral2d")
        return spiral2d(n,d,a)
    
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
    
def coolTrig(n, d):
    trig = 'math.sin(math.cos(math.tan(_x[0]*_x[1]))) - math.sin(math.cos(math.tan(_x[0]))) + math.sin(math.cos(math.tan(_x[1])))'
    return np.array(dg.gibbsSampling(n,d,trig,2,0.1,0.5,2.5))
    
    
def circularWave(n, d):
    cirWave = 'math.sin(math.sqrt(_x[0]**2 + _x[1]**2))'
    return np.array(dg.gibbsSampling(n,d,cirWave,2,0.1,0.5,2.5))
    
def heart(n, d):
    if(d <= 2):
        print("Dimension must be greater than 2 for heart object")
        exit(1)
    
    heart = '(_x[0]**2 + (9*_x[1]**2)/4+_x[2]**2 - 1)**3 - _x[0]**2*_x[2]**3 - (9*_x[1]**2 * _x[2]**3)/80'
    return np.array(dg.gibbsSampling(n,d,heart,2,0.1,0.5,2.5))
    
def heartCurve(n, d):
    heart = '(_x[0]**2 + _x[1]**2 - 1)**3 - _x[0]**2*_x[1]**3'
    return np.array(dg.gibbsSampling(n,d,heart,2,0.1,0.5,2.5))
    
def roseCurve(n, d):
    rose = '(_x[0]**2 + _x[1]**2)**3 - 4*_x[0]**2*_x[1]**2'
    return np.array(dg.gibbsSampling(n,d,rose,2,0.1,0.5,2.5))
    
def butterflyCurve(n, d):
    butterfly = '_x[0]**6 + _x[1]**6 - _x[0]**2'
    return np.array(dg.gibbsSampling(n,d,butterfly,2,0.1,0.5,2.5))
    
def sphere(n, d, r=0.7):
    sphere = '_x[0]**2 + _x[1]**2 + _x[2]**2 - '+str(r)+'**2'
    return np.array(dg.gibbsSampling(n,d,sphere,2,0.1,0.5,2.5))
    
def dsphere(n, d, r=1.0):
    dsphere = ']**2+'.join(['_x['+str(i) for i in range(d)]) + ']**2 -'+str(r)+'**2'
    return np.array(dg.gibbsSampling(n,d,dsphere,2,0.1,0.5,2.5))
    
def hyperboloid(n, d, r=1.0, a=0.2, b=0.5, c=0.8):
    if(d <= 2):
        print("Dimension must be greater than 2 for hyperboloid object")
        exit(1)
    hyperboloid = '(_x[0]**2)/'+str(a)+'**2 + (_x[1]**2)/'+str(b)+'**2 - (_x[2]**2)/'+str(c)+'**2-'+str(r)+'**2'
    return np.array(dg.gibbsSampling(n,d,hyperboloid,2,0.1,0.5,2.5))
    
def hyperbolicParaboloid(n, d, r=1.0, a=0.2, b=0.5, c=0.8):
    if(d <= 2):
        print("Dimension must be greater than 2 for hyperbolic paraboloid object")
        exit(1)
    paraboloid = '-(_x[0]**2)/'+str(a)+'**2 + (_x[1]**2)/'+str(b)+'**2-_x[2]-'+str(r)+'**2'
    return np.array(dg.gibbsSampling(n,d,paraboloid,2,0.1,0.5,2.5))
    
def ellipticCone(n, d, r=0.1, a=0.2, b=0.5, c=0.8):
    if(d <= 2):
        print("Dimension must be greater than 2 for elliptic cone object")
        exit(1)
        
    cone = '(_x[0]**2)/'+str(a)+'**2+(_x[1]**2)/'+str(b)+'**2-(_x[2]**2)/'+str(c)+'**2-'+str(r)+'**2'
    return np.array(dg.gibbsSampling(n,d,cone,2,0.1,0.5,2.5))
    
def circle(n, d, r=0.7):
    circle = '(_x[0]**2) + (_x[1]**2) - '+str(r)+'**2'
    return np.array(dg.gibbsSampling(n,d,circle,2,0.1,0.5,2.5))
    
def sinCurve(n, d, r=1.0):
    sin = '_x[1]-'+str(r)+'*math.sin(_x[0]*10)'
    return np.array(dg.gibbsSampling(n,d,sin,2,0.1,0.5,2.5))

#Commenting out remainder for now, need to confirm further
#def klein(n, d):
#    if(d <= 2):
#        print("Dimension must be greater than 2 for klein object")
#        exit(1)
#    
#    klein = '(_x[0]**2+_x[1]**2+_x[2]**2+2*_x[1]-1)*((_x[0]**2+_x[1]**2+_x[2]**2-2*_x[1]-1)**2-8*_x[2]**2)+16*_x[0]*_x[2]*((_x[0]**2+_x[1]**2+_x[2]**2-2*_x[1]-1))'
#    return np.array(dg.gibbsSampling(n,d,klein,2,0.1,0.5,2.5))
    
#def bulg3d(n, d):
#    if(d <= 2):
#        print("Dimension must be greater than 2 for bulg3d object")
#        exit(1)
        
#    bulg3d = 'math.sqrt(_x[0]**2 + _x[1]**2) * math.sin(math.sqrt(_x[0]**2 + _x[1]**2))/math.sqrt(_x[0]**2 + _x[1]**2)'
#    return np.array(dg.gibbsSampling(n,d,bulg3d,2,0.1,0.5,2.5))
