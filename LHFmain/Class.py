import ctypes

dimension = 3

class point2d(ctypes.Structure):
    _fields_ = ("x0", ctypes.c_double), ("x1", ctypes.c_double)
    

class array2d(ctypes.Structure):
    _fields_ = [("a", ctypes.c_int),
                ("b", ctypes.c_float),
                ("arr", point2d * 10)]


class map(ctypes.Structure):
    _fields_ = [("string 1", ctypes.c_char_p),
                ("string2", ctypes.c_char_p)]

class LHF:
    lib = ctypes.cdll.LoadLibrary("./libLHFlib.so")
    dim = 0
    size = 0
    
    def __init__(self, dim, size):   
        self.dim = dim 
        self.size = size
        
        class point(ctypes.Structure):
            _files = [("x", ctypes.c_double) for x in range(dimension)];
        
        class array2d(ctypes.Structure):
            _fields_ = [("a", ctypes.c_int),
                        ("b", ctypes.c_float),
                        ("arr", point * size)]
        
        self.lib.testFunc.argtypes = [ctypes.c_int]
        self.lib.testFunc.restype = None

        self.lib.pyRunWrapper.argtypes = [map, ctypes.Structure]
        self.lib.pyRunWrapper.restype = None
        

    def pyRunWrapper(self, args, data):
        return self.lib.pyRunWrapper(args, data)

    def testFunc(self, num):
        return self.lib.testFunc(num)

a = LHF(3,2)
a.testFunc(1)
a.pyRunWrapper(map(), array2d())

