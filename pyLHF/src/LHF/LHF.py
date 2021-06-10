import ctypes
import numpy as np
import site 
import os.path


class pipePackett:
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
    def __init__(self,arr1,arr2,arr3):
        self.dim = arr1
        self.birth = arr2
        self.death = arr3

#("dim", ctypes.POINTER(ctypes.c_int)),
class bettiBoundaryTableEntry(ctypes.Structure):
    _fields_ = [  # ("dim", ctypes.POINTER(ctypes.c_int)),
        ("dim", ctypes.c_int),
        ("birth", ctypes.c_double),
        ("death", ctypes.c_double)]


class bettiBoundaryTable(ctypes.Structure):
    _fields_ = [("size", ctypes.c_int),
                ("bettis", ctypes.c_void_p)]


class pipePacketAtt(ctypes.Structure):
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
    
    # Get the actual path of shared library
    
    dllpath = site.getusersitepackages() + os.path.sep + "LHF" + os.path.sep + "libLHFlib.so"
    
    
    # Use RTLD_LAZY mode due to undefined symbols
       
    lib = ctypes.CDLL(dllpath, mode=1)
    args = {"reductionPercentage": "10", "maxSize": "2000", "threads": "30", "threshold": "250", "scalar": "2.0", "mpi": "0", "mode": "standard", "dimensions": "1", "iterations": "250", "pipeline": "", "inputFile": "None",
            "outputFile": "output", "epsilon": "5", "lambda": ".25", "debug": "0", "complexType": "simplexArrayList", "clusters": "20", "preprocessor": "", "upscale": "false", "seed": "-1", "twist": "false", "collapse": "false"}
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

    def __init__(self, data):
        self.args["datadim"] = len(data[0])
        self.args["datasize"] = len(data)

        self.data = data

        self.lib.testFunc.argtypes = [ctypes.c_int, ctypes.c_char_p]
        self.lib.testFunc.restype = None

        self.lib.pyRunWrapper.argtypes = [
            ctypes.c_int, ctypes.c_char_p, ctypes.POINTER(ctypes.c_double)]
        self.lib.pyRunWrapper.restype = None

        # self.lib.pyRunWrapper2.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.POINTER(ctypes.c_double)]
        # self.lib.pyRunWrapper2.restype = ctypes.c_void_p

        LP_c_char = ctypes.POINTER(ctypes.c_char)
        LP_LP_c_char = ctypes.POINTER(LP_c_char)

        self.lib.pyRunWrapper2.argtypes = [
            ctypes.c_int, LP_LP_c_char, ctypes.POINTER(ctypes.c_double)]
        self.lib.pyRunWrapper2.restype = ctypes.c_void_p

    def allocation(size):
        class pybettiBoundaryTableEntry(ctypes.Structure):
            _fields_ = [("dim", ctypes.c_int),
                        #("Bettidim", ctypes.POINTER(ctypes.c_uint)),
                        ("Bettidim", ctypes.c_double * size),
                        ("birth", ctypes.c_double * size),
                        ("death", types.c_double * size)]

    def args2string(self, inList):
        ret = ""

        for a in inList:
            # print(a)
            ret += a+" "
            ret += str(inList[a])+" "
            # print(ret)

        return ret.encode('utf-8')

    def runPH(self):
        # Create char* for passing to C++
        temp = self.args2string(self.args)

        return self.lib.pyRunWrapper(len(temp), ctypes.c_char_p(temp), self.data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))

    def testFunc(self, num, st):
        return self.lib.testFunc(num, ctypes.c_char_p(st.encode('utf-8')))

    def runPH3(self, cmd_input):
        # Create char* for passing to C++
        # att = {"inputFile" : str(cmd_input[1])}
        # print(self.args)
        temp = self.args2string(self.args)

        argv = cmd_input
        argc = len(argv)

        p = ((ctypes.POINTER(ctypes.c_char))*len(argv))()
        for i, arg in enumerate(argv):  # not sys.argv, but argv!!!
            # print(i)

            enc_arg = arg.encode('utf-8')
            p[i] = ctypes.create_string_buffer(enc_arg)

        na = ctypes.cast(p, ctypes.POINTER(ctypes.POINTER(ctypes.c_char)))

        #retPH = pipePacketAtt.from_address(self.lib.pyRunWrapper2(len(temp),ctypes.c_char_p(temp), self.data.ctypes.data_as(ctypes.POINTER(ctypes.c_double))))
        retPH = pipePacketAtt.from_address(self.lib.pyRunWrapper2(argc, na, self.data.ctypes.data_as(ctypes.POINTER(ctypes.c_double))))

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
