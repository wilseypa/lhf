AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = $(ACLOCAL_FLAGS)

bin_PROGRAMS = TDA_Cplusplus
TDA_Cplusplus_SOURCES = main.cpp
TDA_Cplusplus_LDADD = libReadInput.la

lib_LTLIBRARIES = libReadInput.la
libReadInput_la_LDFLAGS = -version-info 0:0:0
libReadInput_la_SOURCES = readInput.cpp readInput.hpp
