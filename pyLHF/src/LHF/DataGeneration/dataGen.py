
import numpy as np
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
