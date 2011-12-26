############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   October 2, 1995
#
#  Modified:	    February 4, 1997  Version 1.0
#		    November 15, 1997 Version 1.1
#		    September 1, 1999 Version 2.0
#                   March    10, 2007 (Only this file) 
#
############################################################################
#
#
#  The name of the libraries to be created/linked to
#
TMGLIB       = libtmglib.a
SUPERLULIB   = libsuperlu.a
#SUPERLULIB   = libsuperlu_3.0.a
BLASLIB      = -lblas

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = /usr/bin/ar
ARCHFLAGS    = cr
RANLIB       = ranlib

CC           = gcc
CFLAGS       = -g -O3
FORTRAN	     = g77
FFLAGS       = -O2
F77_LIBS     =  -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/../../../../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/../../.. -L/lib/../lib64 -L/usr/lib/../lib64 -lfrtbegin -lg2c -lm -lgcc_s
LOADER       = gcc
LOADOPTS     = -xO3
#
#  C preprocessor defs for compilation (-DNoChange, -DAdd_, or -DUpCase)
#
CDEFS        = -DAdd_
#
# The directory in which Matlab is installed
#
MATLAB	     = /usr/sww/pkg/matlab
