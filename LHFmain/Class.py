import ctypes

class point2d(ctypes.Structure):
    _fields_ = ("x0", ctypes.c_double), ("x1", ctypes.c_double)

class setCppU(ctypes.Structure):
    _fields_ = [("setU", ctypes.c_uint)]

class setCppD(ctypes.Structure):
    _fields_ = [("setD", ctypes.c_double)]

class vectordouble(ctypes.Structure):
    _fields_ = [("vectorD", ctypes.c_double)]

class vectorUnsigned(ctypes.Structure):
    _fields_ = [("vectorU", ctypes.c_uint)]

class weights(ctypes.Structure):
    _fields_ = [("setD", ctypes.c_double), ("greater", ctypes.c_double)]

class array2d(ctypes.Structure):
    _fields_ = [("a", ctypes.c_int),
                ("b", ctypes.c_float),
                ("arr", point2d * 10)] ##dynamic



########################################################

# class simplexBase(ctypes.Structure): ###needed? 
#     _fields_ = 

class bettiBoundaryTableEntry(ctypes.Structure):
    _fields_ = [("Bettidim", ctypes.c_uint),
                ("birth", ctypes.c_double),
                ("death", ctypes.c_double),
                ("boundaryPoints", setCppU * 10)] ##dynamic


class pipePacket(ctypes.Structure):
    _fields_ = [("bettiBoundaryTableEntry", bettiBoundaryTableEntry * 10), ##dynamic
                ("ident", ctypes.c_char_p),
                ("stats", ctypes.c_char_p),
                ("runLog", ctypes.c_char_p),
                ("workData", vectordouble * 10), ###dynamic
                ("centroidLabels", vectorUnsigned * 10), ##dynamic
                ("inputData", vectordouble * 10), ##dynamic
                ("distMatrix", vectordouble * 10), ##dynamic
                #("simplexBase", simplexBase * 10),
                ("boundaries", vectorUnsigned * 10), ##dynamic
                ("weights", weights * 10), ##dynamic
                ("bettiOutput", ctypes.c_char_p)] 
                #getSize?
                #getStats?
class map(ctypes.Structure):
    _fields_ = [("string 1", ctypes.c_char_p),
                ("string2", ctypes.c_char_p)]

class LHF(object):
    lib = ctypes.cdll.LoadLibrary("./libLHFlib.so")
    
    def ___init__(self):    
        self.lib.testFunc.argtypes = [ctypes.c_int]
        self.lib.testFunc.restype = None

        self.lib.pyRunWrapper.argtypes = [map, array2d]
        self.lib.pyRunWrapper.restype = None

    def pyRunWrapper(self, map, array2d):
        return self.lib.pyRunWrapper(map, array2d)

    def testFunc(self, int):
        return self.lib.testFunc(int)

a = LHF()
a.testFunc(1)
a.pyRunWrapper(map(), array2d())

