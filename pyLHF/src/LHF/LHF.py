import ctypes
import numpy as np
import os

        
        
class bettiTable(ctypes.Structure):
    
    """
    This structure mirrors the definition of retBettiTable in LHF.hpp
    
    """       
    _fields_ = [
        ("dim", ctypes.c_int),                      # int dim
        ("birth", ctypes.c_double),                 # double birth
        ("death", ctypes.c_double),                 # double death
        ("boundarySize", ctypes.c_int),          # int boundarySize
        ("boundaryEntries", ctypes.POINTER(ctypes.c_uint))    # unsigned* boundaryEntries
        ]
    
    
        
        
        
class pipePacket(ctypes.Structure):
    """
    This structure mirrors the definition of retPipePacket in LHF.hpp
    
    
    """        
    _fields_ = [  
        ("size_betti", ctypes.c_int),                       # int size_betti
        ("LHF_size", ctypes.c_int),                         # int LHF_size
        ("LHF_dim", ctypes.c_int),                          # int LHF_dim
        ("workData_size", ctypes.c_int),                    # int workData_size
        ("bettiTable", ctypes.POINTER(bettiTable)),         # BRET* BettiTable
        ("inputData", ctypes.POINTER(ctypes.c_double)),     # double* inputData
        ("distMatrix", ctypes.POINTER(ctypes.c_double)),    # double* distMatrix
        ("workData", ctypes.POINTER(ctypes.c_double)),      # double* workData
        ("centroidLabels", ctypes.POINTER(ctypes.c_uint)),  # unsigned* centroidLabels
        ("stats", ctypes.c_char_p),                         # char* stats
        ("runLog", ctypes.c_char_p),                        # char* runLog
        ("ident", ctypes.c_char_p)                          # char* ident
        ]

class LHF:
    """
    This class represents the LHF library.

    Attributes:
        lib (ctypes.CDLL): The LHF library.
        args (dict): A dictionary containing the arguments to pass to the LHF function.
        default (dict): A dictionary containing the default arguments for the LHF function.
        data (list): A list containing the data to be used by the LHF function.

    """
    # Use RTLD_LAZY mode due to undefined symbols
    script_dir = os.path.dirname(__file__)
    filename = os.path.join(script_dir, "libLHFlib.so")
    # print(package_dir)
    # lib = ctypes.CDLL("./env/lib/python3.8/site-packages/libLHF/libLHFlib.so", mode=1)
    lib = ctypes.CDLL(filename, mode=1)
    args = {}
    default = {"threads": "30", "mpi": "0", "dimensions": "2", "outputFile": "output", "epsilon": "5", "debug": "0", "complexType": "simplexArrayList", "preprocessor": "", "upscale": "false"}
    data = []

    # Some notes here:
    ##
    # These are returned by C++ in a C function; need to be wrapped into a C structure in LHF
    # This indicates we can return only necessary (C-format) entries - TBD
    ##
    # Some entries have known sizes (based on LHF class):
    # bettiBoundaryTable = (? x 4) - betti boundary struct
    # workData = ((? <= LHF.size) x LHF.dim) - centroids
    # centroidLabels = (LHF.size x 1) - centroid labels
    # inputData = (LHF.size x LHF.dim) - original data
    # distMatrix = (LHF.size x LHF.size) - distance matrix
    ##
    ##

    def __init__(self):
        """
        This class represents MyClass.

        Attributes:
            lib (ctypes.CDLL): The library used by MyClass.
        """

        self.lib.free_pipeWrap.argtypes = [ctypes.c_void_p]
        self.lib.free_pipeWrap.restype = None
        
        self.lib.pyLHFWrapper.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.POINTER(ctypes.c_double)]
        self.lib.pyLHFWrapper.restype = ctypes.c_void_p
        

    def mergeArgs(self):
        """
        This method merges the default arguments and the user-supplied arguments.

        """
        for arg in self.default:
            if arg not in self.args.keys():
                self.args[arg] = self.default[arg]

    def args2string(self, inList):
        """
        Converts a dictionary of arguments into a string.

        Args:
            inList (dict): A dictionary of arguments.

        Returns:
            bytes: A byte string containing the dictionary arguments.
        """
        ret = ""
        for a in inList:
            ret += a + " " + str(inList[a])+" "

        return ret.encode('utf-8')
        
    def runPH(self, data):
        """
        Runs the persistent homology computation using the LHF shared object library and returns the pipepacket result.

        Args:
            data (numpy.ndarray): A 2-dimensional array of data points.

        Returns:
            numpy.ndarray: An array of computed persistent homology data.
        """

        #Get data sizes to pass to C
        self.args["datasize"] = len(data)
        self.args["datadim"] = len(data[0])
        
        self.mergeArgs()
        temp = self.args2string(self.args)
        
        na = ctypes.c_char_p(temp)
        
        self.data = data.flatten().tolist()
        self.data = (ctypes.c_double * len(self.data))(*self.data)
        
        retAddr = self.lib.pyLHFWrapper(len(temp), na, self.data)
        
        retPH = pipePacket.from_address(retAddr)
        
        bettiTable = self.decodeReturn(retPH.size_betti, retPH.bettiTable)
        self.lib.free_pipeWrap(retAddr)
        
        return bettiTable, retPH
        
        
    def decodeReturn(self,s_betti, bettis):
        """
        Decodes the computed persistent homology data returned by the PyBind library.

        Args:
            retPH (pipePacketAtt): A data structure containing the computed persistent homology data.

        Returns:
            numpy.ndarray: An array of decoded persistent homology data.
        """
        piResults = []
        for idx in range(s_betti):
            temp = [bettis[idx].dim, bettis[idx].birth, bettis[idx].death, [bettis[idx].boundaryEntries[i] for i in range(bettis[idx].boundarySize)]]
            
            piResults.append(temp)
        
        return piResults
        

