import ctypes
import numpy as np
import os



class pipePacket:
    """
    PipePacket

    Class that contains information about the pipeline's output data

    Attributes:
    ----------
    betti : np.array
        The betti numbers.
    inputData : np.array
        The input data.
    distMatrix : np.array
        The distance matrix.
    workData : np.array
        The centroids.
    centroidLabels : np.array
        The labels of the centroids.
    stats : str
        The statistics of the computation.
    runLog : str
        The log of the computation.
    ident : str
        The identity of the pipeline.
    """
    def __init__(self,arr1,inputData,distMatrix,workData,centroidLabels,stats,runLog,ident):
        self.betti = arr1
        self.inputData = inputData
        self.distMatrix = distMatrix
        self.workData = workData
        self.centroidLabels = centroidLabels
        self.stats = stats
        self.runLog = runLog
        self.ident =  ident
        #######coming soon
        #bettidim 
        #???

class BettiTable:
    """
    Class that stores betti numbers

    Attributes:
    ----------
    dim : np.array
        The dimension of the betti numbers.
    birth : np.array
        The birth of the betti numbers.
    death : np.array
        The death of the betti numbers.
    """
    def __init__(self,arr1,arr2,arr3):
        self.dim = arr1
        self.birth = arr2
        self.death = arr3

#("dim", ctypes.POINTER(ctypes.c_int)),
class bettiBoundaryTableEntry(ctypes.Structure):
    """
    Class that represents an entry in a betti boundary table.

    Attributes:
    ----------
    dim : int
        The dimension of the betti number.
    birth : double
        The birth of the betti number.
    death : double
        The death of the betti number.
    """
    _fields_ = [  # ("dim", ctypes.POINTER(ctypes.c_int)),
        ("dim", ctypes.c_int),
        ("birth", ctypes.c_double),
        ("death", ctypes.c_double)]


class bettiBoundaryTable(ctypes.Structure):
    """
    Class that represents a betti boundary table.

    Attributes:
    ----------
    size : int
        The size of the betti boundary table.
    bettis : ctypes.c_void_p
        The betti numbers in the betti boundary table.
    """
    _fields_ = [("size", ctypes.c_int),
                ("bettis", ctypes.c_void_p)]


class pipePacketAtt(ctypes.Structure):
    """
    Class that represents the attributes of a pipe packet.

    Attributes:
    ----------
    size : int
        The size of the pipe packet.
    LHF_size : int
        The size of the LHF.
    LHF_dim : int
        The dimension of the LHF.
    workData_size : int
        The size of the work data.
    bettiTable : ctypes.c_void_p
        The betti table.
    inputData : ctypes.c_void_p
        The input data.
    distMatrix : ctypes.c_void_p
        The distance matrix.
    workData : ctypes.c_void_p
        The work data.
    centroidLabels : ctypes.c_void_p
        The centroid labels.
    stats : str
        The statistics of the computation.
    runLog : str
        The log of the computation.
    ident : str
        The identity of the pipeline.
    """
    _fields_ = [("size", ctypes.c_int),
                ("LHF_size", ctypes.c_int),
                ("LHF_dim", ctypes.c_int),
                ("workData_size", ctypes.c_int),
                ("bettiTable", ctypes.c_void_p),
                ("inputData", ctypes.c_void_p),
                ("distMatrix", ctypes.c_void_p),
                ("workData", ctypes.c_void_p),
                ("centroidLabels", ctypes.c_void_p),
                ("stats", ctypes.c_char_p),
                ("runLog", ctypes.c_char_p),
                ("ident", ctypes.c_char_p)]


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

        self.lib.testFunc.argtypes = [ctypes.c_int, ctypes.c_char_p]
        self.lib.testFunc.restype = None

        self.lib.freeWrapper.argtypes = [ctypes.c_void_p]
        self.lib.freeWrapper.restype = None
        
        self.lib.pyRunWrapper.argtypes = [
            ctypes.c_int, ctypes.c_char_p, ctypes.POINTER(ctypes.c_double)]
        self.lib.pyRunWrapper.restype = None

        self.lib.pyRunWrapper2.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.POINTER(ctypes.c_double)]
        self.lib.pyRunWrapper2.restype = ctypes.c_void_p

    def mergeArgs(self):
        """
        This method merges the default arguments and the user-supplied arguments.

        """
        for arg in self.default:
            if arg not in self.args.keys():
                self.args[arg] = self.default[arg]
                
    def allocation(size):
        """
        This method allocates memory for C structures.

        Args:
            size (int): The size of the memory block to allocate.

        Returns:
            None
        """
        class pybettiBoundaryTableEntry(ctypes.Structure):
            _fields_ = [("dim", ctypes.c_int),
                        #("Bettidim", ctypes.POINTER(ctypes.c_uint)),
                        ("Bettidim", ctypes.c_double * size),
                        ("birth", ctypes.c_double * size),
                        ("death", types.c_double * size)]

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

    def runPH(self):
        """
        Runs the persistent homology computation using the PyBind library and returns the computed boundary data.

        Returns:
            numpy.ndarray: An array of computed persistent homology data.
        """
        # Create char* for passing to C++
        self.mergeArgs()
        temp = self.args2string(self.args)

        return self.lib.pyRunWrapper(len(temp), ctypes.c_char_p(temp))

    def testFunc(self, num, st):
        """
        A test function to check the PyBind functionality.

        Args:
            num (int): An integer argument.
            st (str): A string argument.

        Returns:
            bytes: A byte string containing the computed result.
        """
        return self.lib.testFunc(num, ctypes.c_char_p(st.encode('utf-8')))
        
    def runPH(self, data):
        """
        Runs the persistent homology computation using the PyBind library and returns the computed boundary data.

        Args:
            data (numpy.ndarray): A 2-dimensional array of data points.

        Returns:
            numpy.ndarray: An array of computed persistent homology data.
        """

        print("Calling C++ LHF shared library...")
        #Get data sizes to pass to C
        self.args["datasize"] = len(data)
        self.args["datadim"] = len(data[0])
        
        self.mergeArgs()
        
        temp = self.args2string(self.args)
        
        na = ctypes.c_char_p(temp)
        
        self.data = data.flatten().tolist()
        #self.data = ctypes.cast(self.data, ctypes.POINTER(ctypes.c_double))
        
        self.data = (ctypes.c_double * len(self.data))(*self.data)
        
        retAddr = self.lib.pyRunWrapper2(len(temp), na, self.data)
        retPH = pipePacketAtt.from_address(retAddr)
        
        a = self.decodeReturn(retPH)
        self.lib.freeWrapper(retAddr)
        
        return a
        
        
    def decodeReturn(self,retPH):
        """
        Decodes the computed persistent homology data returned by the PyBind library.

        Args:
            retPH (pipePacketAtt): A data structure containing the computed persistent homology data.

        Returns:
            numpy.ndarray: An array of decoded persistent homology data.
        """
        print("Total Boundaries", retPH.size)
        bettiBoundaryTableEntries = type("array", (ctypes.Structure, ), {
            # data members
            "_fields_": [("arr", bettiBoundaryTableEntry * retPH.size)]
        })

        retBounds = bettiBoundaryTableEntries.from_address(retPH.bettiTable)
        
        piResults = []
        for i in range(retPH.size):
            piResults.append([retBounds.arr[i].dim,retBounds.arr[i].birth,retBounds.arr[i].death])
        
        return np.array(piResults)
        
        

    def runPH3(self, cmd_input):
        """
        Runs the persistent homology computation using the PyBind library and returns the computed boundary data.

        Args:
            cmd_input (list): A list of command-line arguments.

        Returns:
            tuple: A tuple of computed persistent homology data.
        """
        
        # Create char* for passing to C++
        # att = {"inputFile" : str(cmd_input[1])}
        # print(self.args)
        temp = self.args2string(self.args)

        argv = cmd_input
        argc = len(argv)

        print("RunPH3", temp)
        
        p = ((ctypes.POINTER(ctypes.c_char))*len(argv))()
        for i, arg in enumerate(argv):  # not sys.argv, but argv!!!
            # print(i)

            enc_arg = arg.encode('utf-8')
            p[i] = ctypes.create_string_buffer(enc_arg)

        #na = ctypes.cast(temp, ctypes.POINTER(ctypes.POINTER(ctypes.c_char)))
        na = ctypes.c_char_p(temp)

        print(na)

        #retPH = pipePacketAtt.from_address(self.lib.pyRunWrapper2(len(temp),ctypes.c_char_p(temp), self.data.ctypes.data_as(ctypes.POINTER(ctypes.c_double))))
        retPH = pipePacketAtt.from_address(self.lib.pyRunWrapper2(argc, na))

        print("Total Boundaries", retPH.size)
        # print(retPH.ident)

        # Reconstruct the boundary table array from the address?

        bettiBoundaryTableEntries = type("array", (ctypes.Structure, ), {
            # data members
            "_fields_": [("arr", bettiBoundaryTableEntry * retPH.size)]
        })

        retBounds = bettiBoundaryTableEntries.from_address(retPH.bettiTable)

        # for i in range(retPH.size):
            # print(retBounds.arr[i].dim,retBounds.arr[i].birth,retBounds.arr[i].death)

        ###############################################################################
        # Reconstruct the inputData array from the address
        inputDataEntries = type("array", (ctypes.Structure, ), {
            # data members
            "_fields_": [("arr", ctypes.c_double * (retPH.LHF_size * retPH.LHF_dim))]
        })

        retinputData = inputDataEntries.from_address(retPH.inputData)

        # for i in range(retPH.LHF_size * retPH.LHF_dim):
        #     print(i, ": ", retinputData.arr[i])

        inputData = [[0]*retPH.LHF_dim]*retPH.LHF_size
        sizof = 0
        for i in range(retPH.LHF_size):
            for j in range(retPH.LHF_dim):
                # print(retinputData.arr[sizof])
                inputData[i][j] = retinputData.arr[sizof]
                sizof = sizof + 1

        # print(inputData)
        ###############################################################################

        distMatrixEntries = type("array", (ctypes.Structure, ), {
            # data members
            "_fields_": [("arr", ctypes.c_double * (retPH.LHF_size * retPH.LHF_size))]
        })

        retdistMatrix = distMatrixEntries.from_address(retPH.distMatrix)

        distMatrix = [[0]*retPH.LHF_size]*retPH.LHF_size
        sizof = 0
        for i in range(retPH.LHF_size):
            for j in range(retPH.LHF_size):
                # print(retdistMatrix.arr[sizof])
                distMatrix[i][j] = retdistMatrix.arr[sizof]
                sizof = sizof + 1

        # print(distMatrix)
        ##################################################################################

        centroidLabelsEntries = type("array", (ctypes.Structure, ), {
            # data members
            "_fields_": [("arr", ctypes.c_double * (retPH.LHF_size)*1)]
        })

        retcentroidLabels = centroidLabelsEntries.from_address(retPH.centroidLabels)

        # print(type(retcentroidLabels.arr[0]))
        if isinstance(retcentroidLabels.arr[0], int):
            # for i in range(retPH.LHF_size * 1):
            #     print(i, ": ", retcentroidLabels.arr[i])

            centroidLabels = [[0]*1]*retPH.LHF_size
            sizof = 0
            for i in range(retPH.LHF_size):
                for j in range(1):
                    # print(retdistMatrix.arr[sizof])
                    centroidLabels[i][j] = retcentroidLabels.arr[sizof]
                    sizof = sizof + 1

            # print(centroidLabels)
        else:
            centroidLabels = ""
        # for i in range(retPH.LHF_size * 1):
        #     print(i, ": ", retcentroidLabels.arr[i])

        ####################################################################################

        workDataEntries = type("array", (ctypes.Structure, ), {
            # data members
            "_fields_": [("arr", ctypes.c_double * (retPH.workData_size * retPH.LHF_dim))]
        })

        retworkData = workDataEntries.from_address(retPH.workData)

        workData = [[0]*retPH.LHF_dim]*retPH.LHF_size
        sizof = 0
        for i in range(retPH.LHF_size):
            for j in range(retPH.LHF_dim):
                # print(retdistMatrix.arr[sizof])
                workData[i][j] = retworkData.arr[sizof]
                sizof = sizof + 1

        # print(type(retworkData.arr[0]))
        # print(retworkData.arr[0])
        # for i in range(retPH.workData_size * retPH.LHF_dim):
        #     print(i, ": ", retworkData.arr[i])

        #####################################################################################

        pystats = str(retPH.stats).replace(
            '\\n', '\n').replace('b\'', '').replace('\'', '')
        # pystats = retPH.stats.decode('UTF-8')
        # print(pystats)

        pyrunLog = str(retPH.runLog).replace(
            '\\n', '\n').replace('b\'', '').replace('\'', '')
        # pyrunLog = retPH.runLog.decode('UTF-8')
        # print(pyrunLog)

        pyident = str(retPH.ident).replace(
            '\\n', '\n').replace('b\'', '').replace('\'', '')
        # pyident = retPH.ident.decode('UTF-8') ##does not work when c++ returns None
        # print(pyident)

        ########################################################################################
        ### Reconstruct the return
        dim = []
        birth = []
        death = []
        for i in range(retPH.size):
            # print(retBounds.arr[i].dim,retBounds.arr[i].birth,retBounds.arr[i].death)
            dim.append(retBounds.arr[i].dim)
            birth.append(retBounds.arr[i].birth)
            death.append(retBounds.arr[i].death)


        
        result = BettiTable(np.asarray(dim),np.asarray(birth),np.asarray(death))
        result2 = pipePackett(result,np.asarray(inputData),np.asarray(distMatrix),np.asarray(workData),np.asarray(centroidLabels),pystats,pyrunLog,pyident)


        return result2
