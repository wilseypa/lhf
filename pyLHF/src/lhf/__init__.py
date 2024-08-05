from . import DataGeneration
from . import OutputAnalysis
from . import utilities
from . import triangulation
from . import curves


import ctypes
from time import process_time
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

class createPipeline:
    """
    This class represents building a pipeline using the LHF library.

    Attributes:
        lib (ctypes.CDLL): The LHF library.
        args (dict): A dictionary containing the arguments to pass to the LHF function.
        default (dict): A dictionary containing the default arguments for the LHF function.
        data (list): A list containing the data to be used by the LHF function.

    """
    # Use RTLD_LAZY mode due to undefined symbols
    script_dir = os.path.dirname(__file__)
    filename = os.path.join(script_dir, "libLHFlib.so")
    lib = ctypes.CDLL(filename, mode=1)
    config = {}
    default = {"threads": "30", "mpi": "0", "dimensions": "2", "outputFile": "output", "epsilon": "5", "debug": "0", "complexType": "simplexArrayList", "preprocessor": "", "upscale": "false"}
    data = []



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
        

    def mergeConfig(self):
        """
        This method merges the default configuration arguments and the user-supplied configuration.

        """
        for cfg in self.default:
            if cfg not in self.config.keys():
                self.config[cfg] = self.default[cfg]

    def config2string(self, inList):
        """
        Converts a dictionary of configuration arguments into a string.

        Args:
            inList (dict): A dictionary of configuration arguments.

        Returns:
            bytes: A byte string containing the dictionary configuration.
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
        self.config["datasize"] = len(data)
        self.config["datadim"] = len(data[0])
        
        self.mergeConfig()
        temp = self.config2string(self.config)
        
        na = ctypes.c_char_p(temp)
        
        self.data = data.flatten().tolist()
        self.data = (ctypes.c_double * len(self.data))(*self.data)
        
        start = process_time()
        retAddr = self.lib.pyLHFWrapper(len(temp), na, self.data)
        elapsed = (process_time() - start)
        
        retPH = pipePacket.from_address(retAddr)
        
        bettiTable = self.decodeReturn(retPH.size_betti, retPH.bettiTable)
        self.lib.free_pipeWrap(retAddr)
        
        return bettiTable, retPH, elapsed
        
        
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
        

