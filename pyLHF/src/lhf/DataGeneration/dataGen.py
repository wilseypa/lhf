
import numpy as np
from datetime import datetime
import math
import random
from scipy.stats import qmc
from scipy.stats import ortho_group

def genFilledCube (dim=3, k=12, range=(-1.0, 1.0)) :
    """
    Generate a nearly uniform distribution of points in a cube using the sobol sequence generator in scipy.  

    Parameters
    ----------
    dim : int 
        The number of dimensions for the desired cube (default 3).
    k : int 
        A power of 2 that defines the number of points to generate (default 11, for 2**12 total points).  The maximum value for this is 31.
    range : (float, float) 
        The range of values to generate for each coordinate axis (default (-1.0, 1.0)).

    Returns
    -------
    numpy.array(2**k, 3)
        numpy matrix of the generated data.
    """
    return np.array(qmc.scale(qmc.Sobol(d=dim).random_base2(m=k), [range[0]] * dim, [range[1]] * dim))


def buildObj (objDfn='False', dataIn=[], center=None) :
    """
    Select points from a matrix of data (dataIn) that satisfy an evaluation of the objDfn.  objDfn is a boolean function (or
    string that evaluates to a boolean) that is
    evaluated against each row of dataIn.  If objDfn returns true, the row is appended to the output data.  This function assumes
    that the evaluation string will reference the elements of the row by the variable name "_x[i]" where i is the index.  

    In general object construction (and the objDfn equation) is assumed to be centered at the origin.  The center parameter
    permits the user to define a different location for the object.

    Parameters
    ----------
    objDfn : boolean function/string
        A python boolean function or expression (presented as a string)  (default False).
    dataIn : numpy.array(*, *)
        A matix of points to be considered for inclusion within the expression defining the object.
    center : (*)
        The cartesian coordinates of the center of the object (default [0]*dimension)

    Returns
    -------
    numpy.array(2**k, 3)
        numpy matrix of the input data that satisfies the objDfn boolean.
    """
    outData = []
    if center is None :
        center = np.zeros(len(dataIn[0]))
    for d in dataIn :
        _x = d - center
        if eval(objDfn) : outData.append(d)
    return np.array(outData)


def embedding(data, targetDim=3, center=None) :
    """
    Embed data in an abient space of dimension, optionally with the data recentered.  The input data is assumed to be centered at
    the orgin of object space.  If otherwise unspecified (by the center parameter), the additional dimensional coordinates are set
    to zero. 

    Parameters
    ----------
    data : numpy.array(n, m)
        m dimensional data to embed/recenter.
    targetDim : int 
        The ambient dimension, must be >= object dimension (default 3).
    center : [float] * m
        The coordinates in the ambient space for the recentering.

    Returns
    -------
    numpy.array(2**k, 3)
        numpy matrix of the embedded data
    """

    n, d = data.shape
    if targetDim < d :
        print('Aborting: Can only embed a point cloud into the same or higher dimensional space.  Source Dim: {}, Target Dim: {}'.format(d,targetDim))
        sys.exit(-1)
    if center is None :
        if targetDim == d :
            return data
    
    # Embed the point cloud in the target dimension by setting all additional axes to be 0
    newPtCld = np.zeros((n, targetDim))
    newPtCld[:,:d] = data
    
    # if requested, recenter the data
    if center is None :
        return newPtCld
    else :
        return newPtCld + center

    
def rotate(data, rotation='random') :
    """
    Perform a rotation of the data.  Currently only supporting random linear rotations.

    Parameters
    ----------
    data : numpy.array(n, m)
        The m-dimensional data to rotate.
    rotation : string 
        Type of rotation to perform; options are: 'random' (default 'random').

    Returns
    -------
    numpy.array(n, m)
        rotated numpy matrix of the input data.
    """
    if rotation == 'random' :
        # Generate a random orthogonal marix to rotate the point cloud and then
        # Transform the point cloud using matrix multiplication
        return (np.dot(data, ortho_group.rvs(dim = len(data[0]))))
    else :
        errStr = 'Unknown rotation request: {}'.format(rotation)
        sys.exit(errStr)



def sparseEnough(points, current, sparseSphereR):
    """
    Check if the new point is sparse enough - i.e. if it's outside of other points with a sphere around them of radius sparseSphereR

    Parameters
    ----------
    points : numpy.array(l, m)
        The input points to test for sparse
    current : numpy.array(n, m)
        The currently generated data
    sparseSpherR : float
        dispersion factor for moving around the mainfold surface
    
    Returns
    -------
    bool
        True if the point is valid (sparse), False if the point is too close to another point (dist < sparseSphereR)
    """
    
    for x in points:
        if(math.dist(x, current) < sparseSphereR):
            return False
    return True

def normal(current, eps, seed = datetime.now().timestamp()):
    """
    ...
    
    Parameters
    ----------
    current : numpy.array(n, m)
        The currently generated data
    eps : float
        Dilation parameter - smaller values for closer to surface, larger values for rough surface
    seed : random.seed()
        Seed parameter for controlling the function
    
    Returns
    -------
    bool
        True if the point is valid (sparse), False if the point is too close to another point (dist < sparseSphereR)
    """
    random.seed(seed)
    val = np.random.normal(loc=current, scale=eps, size=None)
    
    if(val<5 and val>-5):
        return val
    else:
        return current
        
def actFunction(x, disp):
    return 1-math.sqrt(1-x**disp)


def funcEval(a, x):
    _x = x    
    return eval(a)
    

def validation1(a, current, eps, disp, seed=datetime.now().timestamp()):
    """
    ...
    
    Parameters
    ----------
    a : ...
        TBD...
    current : numpy.array(n, m)
        The currently generated data
    eps : float
        Dilation parameter - smaller values for closer to surface, larger values for rough surface
    disp : float
        dispersion factor for moving around the mainfold surface
    
    Returns
    -------
    bool
        True if the point is valid, False otherwise
    """
    x = np.array(current)
    dilationVector = (x * eps) / np.linalg.norm(x)
    
    x1 = x + dilationVector
    x2 = x - dilationVector
    
    mag1 = funcEval(a, x + dilationVector)
    mag2 = funcEval(a, x - dilationVector)
    
    if(mag1 * mag2 > 0):
        return False
    else:
        if mag1 < 0: mag1 *=-1
        if mag2 < 0: mag2 *=-1
        if mag2 > mag1:
            x = mag1 / mag2
        else:
            x = mag2 / mag1
            
        acceptProb = actFunction(x, disp)
        
        random.seed(seed)
        if(np.random.uniform(0,1) < acceptProb):
            return True
        else:
            return False
            
def validation2(a, current, eps, disp, seed=datetime.now().timestamp()):
    """
    ...
    
    Parameters
    ----------
    a : ...
        TBD...
    current : numpy.array(n, m)
        The currently generated data
    eps : float
        Dilation parameter - smaller values for closer to surface, larger values for rough surface
    disp : float
        dispersion factor for moving around the mainfold surface
    
    Returns
    -------
    bool
        True if the point is valid, False otherwise
    """
    x = np.array(current)
    res = funcEval(a, x)
    
    if(abs(res) > eps):
        return False
    else:
        
        x = 0
        if res > 0:
            x = 1-(res/eps)
        else:
            x = 1+(res/eps)
        
        acceptProb = actFunction(x, disp)
        
        random.seed(seed)
        if(np.random.uniform(0,1) < acceptProb):
            return True
        else:
            return False
            
        

def gibbsSampling(n, d, a, disp, eps, step, maxStep = 2.5, maxHeat = 1000, sparseSphereR = 0.01):
    """
    Perform gibbs sampling from a manifold:
        .Sample data from a parametric equation wth dilation (noise) epsilon
        .Works well for higher dimensional spaces
        .Crawl the surface of the manifold in epsilon-vicinity
        .Minimize potential to crawl away from parametric equation

    Parameters
    ----------
    n : int
        The number of points to generate
    d : int
        The dimension of data to generate
    a : ...
        TBD...
    disp : float
        dispersion factor for moving around the mainfold surface
    eps : float
        Dilation parameter - smaller values for closer to surface, larger values for rough surface
    neigh: float
        Size of jump in neighborhood - smaller values result in localized exploring, larger may lead to more rejections
    maxJump: float
        Maximum size of a jump when walking the manifold
    
    
    Returns
    -------
    numpy.array(n, m)
        numpy matrix of generated data.
    ?? float
        n/i (number of points over iterations)
    """
    
    outData = []
    
    #Initialize to a random point from R^{dim}
    current = [0 for i in range(d)]
    
    onManifold = False
    heat = 0
    stepS = step
    indx = 0
    
    while(True):
        newCurrent = current.copy()
        newCurrent[indx] = normal(newCurrent[indx], stepS)
        a1 = validation1(a, newCurrent, eps, disp)
        a2 = validation2(a, newCurrent, eps, disp)
        
        if((a1 or a2) and sparseEnough(outData, newCurrent, sparseSphereR)):
            outData.append(newCurrent)
            current = newCurrent
            onManifold = True
            heat = 0
            stepS = step
        
        else:
            heat += 1
            
            if heat > maxHeat and stepS < maxStep:
                stepS = stepS * 2
                heat = 0
            else:
                pass
                
            if not onManifold:
                current = newCurrent
            else:
                pass
        
        if(len(outData) == n):
            break;
        
        indx = (indx+1) % d
        
    return outData
        
    
    
    
    
    


