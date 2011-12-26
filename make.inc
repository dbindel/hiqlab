INCLUDE = -I$(TOOLS)/UFconfig/ \
	-I$(TOOLS)/AMD/Include -I$(TOOLS)/UMFPACK/Include \
	-I$(TOOLS)/lua/include -I$(TOOLS)/tolua++/include 

CC	 = gcc
CXX	 = g++
CFLAGS	 = -g -O2   $(INCLUDE)
CXXFLAGS = -g -O2 $(INCLUDE)
CYGPATH  = echo
HIQINT32 = int
F77      = gfortran-mp-4.4
FFLAGS   = -g -O2

AR	 = /usr/bin/ar
RANLIB	 = ranlib
AWK	 = gawk

MEX	 = /Applications/MATLAB_R2011a.app/bin/mex -largeArrayDims
MEXEXT	 = mexmac
MEXOCT   = mkoctfile --mex
#MEXEXT   = mex

SRCDIR   = /Users/dbindel/work/hiqlab-new/src/
BINDIR   = /Users/dbindel/work/hiqlab-new/bin/
LIBDIR   = /Users/dbindel/work/hiqlab-new/lib/
TOOLS    = /Users/dbindel/work/hiqlab-new/tools/
MWRAP	 = $(TOOLS)/mwrap/mwrap
MATEXPR	 = $(TOOLS)/matexpr/matexpr
TOLUA	 = $(TOOLS)/tolua++/bin/tolua
DSBWEB   = $(TOOLS)/dsbweb/dsbweb

LIBQLAB    = $(LIBDIR)/libqlab.a
HIQ        = $(BINDIR)/hiqlab
HIQ_HYBRID = $(SRCDIR)/lua_hybrid/hiqlab

FLTK_CONFIG   = 
FLTK_CXXFLAGS = 
FLTK_LIBS     = 
#MESHVFLTK_BUILD = meshvfltk.build

GD_CFLAGS = 
GD_LIBS   =  

READLINE_LIBS = -lreadline
LUALIBS	   = -L$(TOOLS)/tolua++/lib -ltolua \
		-L$(TOOLS)/lua/lib -llua -llualib

UMF_LIBS = $(TOOLS)/UMFPACK/Lib/libumfpack.a \
	   $(TOOLS)/AMD/Lib/libamd.a

LAPACK_LIBS = 
BLAS_LIBS   = -Wl,-framework -Wl,vecLib
ARLIBS      = $(TOOLS)/ARPACK/libarpack.a

F77_LIBS    =  -L/opt/local/lib/gcc44/gcc/x86_64-apple-darwin11/4.4.6 -L/opt/local/lib/gcc44/gcc/x86_64-apple-darwin11/4.4.6/../../.. -lgfortranbegin -lgfortran
MF77_LIBS   =  -L/opt/local/lib/gcc44/gcc/x86_64-apple-darwin11/4.4.6 -L/opt/local/lib/gcc44/gcc/x86_64-apple-darwin11/4.4.6/../../.. -lgfortranbegin -lgfortran
OTHER_LIBS  = -lm 

LALIBS      = $(ARLIBS) $(UMF_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS)
MLALIBS     = $(ARLIBS) $(UMF_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS)

LIBS       = $(LUALIBS) $(LALIBS)  $(F77_LIBS)  $(OTHER_LIBS)
MEXLIBS    = $(LUALIBS) $(MLALIBS) $(MF77_LIBS) $(OTHER_LIBS)

VERSION	   = 0.4


# -- Specifically for compiling MEX files with the mingw32 toolchain

WMEXINC    = $(TOOLS)/mex-mingw32/mexinc
WMEXLIBDIR = $(TOOLS)/mex-mingw32/mexlib

WMEXLIBS   = $(WMEXLIBDIR)/fixup.o \
	$(WMEXLIBDIR)/mexlib1.lib \
	$(WMEXLIBDIR)/mexlib2.lib \
	$(WMEXLIBDIR)/mexlib3.lib 

WMEXFLAGS=\
        -DMATLAB_MEX_FILE -I $(WMEXINC) -O3 \
        -malign-double -fno-exceptions -mcpu=pentium 
