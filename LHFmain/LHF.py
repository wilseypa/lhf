import ctypes

class bettiBoundaryTableEntry(ctypes.Structure):
	_fields_ = [("Bettidim", ctypes.c_uint),
			("birth", ctypes.c_double),
			("death", ctypes.c_double),
			("boundaryPoints", ctypes.POINTER(ctypes.c_uint))] ##dynamic

class LHF:
    lib = ctypes.cdll.LoadLibrary("./libLHFlib.so")
    args = {"reductionPercentage":"10","maxSize":"2000","threads":"30","threshold":"250","scalar":"2.0","mpi": "0","mode": "standard","dimensions":"1","iterations":"250","pipeline":"","inputFile":"None","outputFile":"output","epsilon":"5","lambda":".25","debug":"0","complexType":"simplexArrayList","clusters":"20","preprocessor":"","upscale":"false","seed":"-1","twist":"false","collapse":"false"}    
    data = []

    ## Some notes here:
    ##
    ##  These are returned by C++ in a C function; need to be wrapped into a C structure in LHF 
    ##      This indicates we can return only necessary (C-format) entries - TBD
    ##
    ##  Some entries have known sizes (based on LHF class):
    ##          bettiBoundaryTable = (? x 4) - betti boundary struct
    ##          workData = ((? <= LHF.size) x LHF.dim) - centroids
    ##          centroidLabels = (LHF.size x 1) - centroid labels
    ##          inputData = (LHF.size x LHF.dim) - original data
    ##          distMatrix = (LHF.size x LHF.size) - distance matrix
    ##          
    ##
    class pipePacket(ctypes.Structure):
        _fields_ = [("bettiBoundaryTableEntry", ctypes.POINTER(bettiBoundaryTableEntry)), ##dynamic
                    ("ident", ctypes.c_char_p),
                    ("stats", ctypes.c_char_p),
                    ("runLog", ctypes.c_char_p),
                    ("workData", ctypes.POINTER(ctypes.c_double)), ###dynamic
                    ("centroidLabels", ctypes.POINTER(ctypes.c_double)), ##dynamic
                    ("inputData", ctypes.POINTER(ctypes.c_uint)), ##dynamic
                    ("distMatrix", ctypes.POINTER(ctypes.c_double)), ##dynamic
                    ("boundaries", ctypes.POINTER(ctypes.c_uint)), ##dynamic
                    ("weights", ctypes.POINTER(ctypes.c_double)), ##dynamic
                    ("bettiOutput", ctypes.c_char_p)] 
    
    
    def __init__(self, data):  
        self.args["datadim"] = len(data[0])
        self.args["datasize"] = len(data)
        
        self.data = data;
        
        self.lib.testFunc.argtypes = [ctypes.c_int, ctypes.c_char_p]
        self.lib.testFunc.restype = None

        self.lib.pyRunWrapper.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.POINTER(ctypes.c_double)]
        self.lib.pyRunWrapper.restype = None
   
    def args2string(self, inList):
        ret = ""
        
        for a in inList:
            ret += a+" "
            ret += str(inList[a])+" "
            
        return ret.encode('utf-8')
        
    def runPH(self):
        #Create char* for passing to C++
        temp = self.args2string(self.args)
        
        return self.lib.pyRunWrapper(len(temp),ctypes.c_char_p(temp), self.data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))


    def testFunc(self, num, st):
        return self.lib.testFunc(num, ctypes.c_char_p(st.encode('utf-8')))


