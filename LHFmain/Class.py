import ctypes

class point2d(ctypes.Structure):
    _fields_ = ("x0", ctypes.c_double), ("x1", ctypes.c_double)

class array2d(ctypes.Structure):
    _fields_ = [("a", ctypes.c_int),
                ("b", ctypes.c_float),
                ("arr", point2d * 10)]


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
        return self.lib.testFunc(map, array2d)

    def testFunc(self, int):
        return self.lib.pyRunWrapper(int)

a = LHF()
a.testFunc(1)
a.pyRunWrapper(map(), array2d())

