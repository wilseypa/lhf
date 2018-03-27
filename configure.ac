## --AutoConf--
# Process with autoconf to produce configuration script

AC_PREREQ([2.69])
AC_INIT([LHF], [0.0.1], [wilseypa@mail.uc.edu])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([m4])
AC_LANG([C++])

# Check for programs
AC_PROG_CXX
AC_PROG_CC

# Check for libraries
LT_INIT
AC_ENABLE_SHARED
AC_DISABLE_STATIC
AC_PROG_LIBTOOL(libtool)

# Check for typedefs, structs, and compiler characteristics
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Output Makefile
AC_CONFIG_FILES([Makefile test/Makefile src/Makefile])
AC_OUTPUT
