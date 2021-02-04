import ctypes

lib = ctypes.cdll.LoadLibrary("../lhf/LHFmain/libLHFlib.so")

class pipepacket(Structure):
    ##
    ##

class map(Structure):
    _fields_ = [("string 1", c_string),
                ("string2", c_string)]

class LHF(object):

    lib.outputBettis.argtypes = [map pipepacket] ###pipepacket??
    lib.outputBettis.restype = None

    lib.runPipeLine.argtypes = [map, pipepacket] ###pipepacket??
    lib.runPipeLine.restype = None

    lib.processDataWrapper.argtypes = [map, pipepacket] ###pipepacket???
    lib.processDataWrapper.restype = None

    lib.processIterupScale.argtypes = [map, pipepacket] ###pipepacket???
    lib.processIterupScale.restype = None ##pointer

    lib.processUpscaleWrapper.argtypes = [map, pipepacket] ###pipepacket
    lib.processUpscaleWrapper.restype = None ##pointer

    lib.new_vector.restype = c_void_p
    lib.new_vector.argtypes = []

    def _init_(self):
        self.vector = LHF.lib.new_vector() ##pointer to the BettiTable vector

    def outputBettis(map, pipepacket):
        LHF.lib.outputBettis(map, pipepacket)
    
    def runPipeLine(map, pipepacket):
        LHF.lib.runPipeLine(map, pipepacket)

    def processDataWrapper(map, pipepacket):
        LHF.lib.processDataWrapper(map, pipepacket)

    def processIterupScale(map, pipepacket):
        LHF.lib.processIterupScale(self.vector, map, pipepacket)

    def processUpscaleWrapper(map, pipepacket):
        LHF.lib.processUpscaleWrapper(self.vector, map, pipepacket)