AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = $(ACLOCAL_FLAGS)
AM_CXXFLAGS = -fpic -Wall -Wextra -std=c++17

bin_PROGRAMS = TDA_Cplusplus
TDA_Cplusplus_SOURCES = main.cpp
TDA_Cplusplus_LDADD = libPipePacket.la libReadInput.la libWriteOutput.la libArgParser.la libBasePipe.la libDistMatrixPipe.la libNeighGraphPipe.la

lib_LTLIBRARIES = libPipePacket.la libReadInput.la libWriteOutput.la libArgParser.la libBasePipe.la libDistMatrixPipe.la libNeighGraphPipe.la
libPipePacket_la_LDFLAGS = -version-info 0:0:0
libPipePacket_la_SOURCES = pipePacket.cpp pipePacket.hpp

libReadInput_la_LDFLAGS = -version-info 0:0:0
libReadInput_la_SOURCES = readInput.cpp readInput.hpp

libWriteOutput_la_LDFLAGS = -version-info 0:0:0
libWriteOutput_la_SOURCES = writeOutput.cpp writeOutput.hpp

libArgParser_la_LDFLAGS = -version-info 0:0:0
libArgParser_la_SOURCES = argParser.cpp argParser.hpp

libBasePipe_la_LDFLAGS = -version-info 0:0:0
libBasePipe_la_SOURCES = basePipe.cpp basePipe.hpp

libDistMatrixPipe_la_LDFLAGS = -verion-info 0:0:0
libDistMatrixPipe_la_SOURCES = distMatrixPipe.cpp distMatrixPipe.hpp

libNeighGraphPipe_la_LDFLAGS = -version-info 0:0:0
libNeighGraphPipe_la_SOURCES = neighGraphPipe.cpp neighGraphPipe.hpp
