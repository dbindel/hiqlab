include ../../make.inc

TRI_DIR=/usr/local/trilinos/7.0.5
TRI_VER=linux_mpi1_gnu
EPETRA_DIR=${TRI_DIR}/${TRI_VER}/packages/epetra
EPETRA_CON=${EPETRA_DIR}/src
EPETRA_INC=${TRI_DIR}/packages/epetra/src
EPETRA_LIB=${EPETRA_DIR}/src

DIDASKO_DIR=${TRI_DIR}/${TRI_VER}/packages/didasko
DIDASKO_CON=${DIDASKO_DIR}/src
DIDASKO_INC=${TRI_DIR}/packages/didasko/src
DIDASKO_LIB=${DIDASKO_DIR}/src

AMESOS_DIR=${TRI_DIR}/${TRI_VER}/packages/amesos
AMESOS_CON=${AMESOS_DIR}/src
AMESOS_INC=${TRI_DIR}/packages/amesos/src
AMESOS_LIB=${AMESOS_DIR}/src

TEUCHOS_DIR=${TRI_DIR}/${TRI_VER}/packages/teuchos
TEUCHOS_CON=${TEUCHOS_DIR}/src
TEUCHOS_INC=${TRI_DIR}/packages/teuchos/src
TEUCHOS_LIB=${TEUCHOS_DIR}/src

TRIUTILS_DIR=${TRI_DIR}/${TRI_VER}/packages/triutils
TRIUTILS_CON=${TRIUTILS_DIR}/src
TRIUTILS_INC=${TRI_DIR}/packages/triutils/src
TRIUTILS_LIB=${TRIUTILS_DIR}/src

EPETRAEXT_DIR=${TRI_DIR}/${TRI_VER}/packages/epetraext
EPETRAEXT_CON=${EPETRAEXT_DIR}/src
EPETRAEXT_INC=${TRI_DIR}/packages/epetraext/src
EPETRAEXT_LIB=${EPETRAEXT_DIR}/src

MUMPS_INC = /usr/local/mumps/4.6.4/include/

LIB_NP    = /usr/local/mumps/4.6.4/lib/libdmumps.a \
            /usr/local/mumps/4.6.4/lib/libzmumps.a \
            /usr/local/mumps/4.6.4/lib/libpord.a \
            /usr/local/parmetis/3.1/libparmetis.a \
            /usr/local/parmetis/3.1/libmetis.a \
            /usr/local/scalapack/1.7.5/libscalapack.a \
            /usr/local/mpiblacs/BLACS/LIB/libblacs.a \
            /usr/local/mpiblacs/BLACS/LIB/libblacsF77init.a \
            /usr/local/mpiblacs/BLACS/LIB/libblacs.a

LIB_GFORTRAN= /usr/lib/gcc/x86_64-redhat-linux/4.1.0/libgfortranbegin.a \
              /usr/lib/gcc/x86_64-redhat-linux/4.1.0/libgfortran.a


TRI_LIBS= -L${DIDASKO_LIB} \
          ${AMESOS_LIB}/libamesos.a \
          ${TRIUTILS_LIB}/libtriutils.a \
          ${TEUCHOS_LIB}/libteuchos.a  \
          ${EPETRAEXT_LIB}/libepetraext.a \
          ${EPETRA_LIB}/libepetra.a 

LIBQLAB=../corelib/libqlab.a
LIBQLAB_TL = ../corelib_trilinos/libqlab_trilinos.a

OBJS= 

all: ex9.exe ex12.exe ex1.exe

.cc.o: 
	$(CXX) -DHAVE_CONFIG_H $(CXXFLAGS) -I../corelib \
               -I${EPETRA_CON} -I${EPETRA_INC}  \
               -I${DIDASKO_CON} -I${DIDASKO_INC}  \
               -I${AMESOS_CON} -I${AMESOS_INC}  \
               -I${TEUCHOS_CON} -I${TEUCHOS_INC}  \
               -I${TRIUTILS_CON} -I${TRIUTILS_INC}  \
               -I${EPETRAEXT_CON} -I${EPETRAEXT_INC}  \
               -I${MUMPS_INC}  \
               -I. -I.. -I../corelib_trilinos -c $<

ex12.o: ex12.cc
ex9.o: ex9.cc
trilinos_super_matrix.o: trilinos_super_matrix.cc trilinos_super_matrix.h
trilinos_super_vector.o: trilinos_super_vector.cc trilinos_super_vector.h

ex1.exe: ex1.o 
	$(CXX) -DHAVE_CONFIG_H -o ex1.exe ex1.o \
                $(OBJS) $(LIBQLAB_TL) $(LIBQLAB)\
                $(TRI_LIBS) $(LIB_NP) $(LIB_GFORTRAN) \
                $(LIBQLAB_TL) $(LIBQLAB_TL) \
		$(READLINE_LIBS) $(FLTK_LIBS) $(LIBS)

ex9.exe: ex9.o 
	$(CXX) -DHAVE_CONFIG_H -o ex9.exe ex9.o \
                $(OBJS) $(LIBQLAB_TL) $(LIBQLAB)\
                $(TRI_LIBS) \
		$(READLINE_LIBS) $(FLTK_LIBS) $(LIBS)

ex12.exe: ex12.o trilinos_super_matrix.o trilinos_super_vector.o
	$(CXX) -DHAVE_CONFIG_H -o ex12.exe ex12.o trilinos_super_matrix.o \
                trilinos_super_vector.o \
                $(OBJS) $(LIBQLAB_TL) $(LIBQLAB)\
                $(TRI_LIBS) \
		$(READLINE_LIBS) $(FLTK_LIBS) $(LIBS)

clean:
	rm -f *lua.* *.o *~ ex12.exe ex9.exe
