include ../../make.inc

PETSC_DIR=/usr/local/petsc/2.3.2-p6
PETSC_ARCH=linux_mpi1_cxx_hypre
LIBSUPERLU=/usr/local/superlu/serial_3.0/libsuperlu.a

CFLAGS	         := $(CFLAGS) -I../corelib -I../lua 
CXXFLAGS	 := $(CXXFLAGS) -I../corelib -I../lua 
FFLAGS	         = 
CPPFLAGS         =
FPPFLAGS         =
NP               = 1
DISPLAY          = :0.0
PETSC_LIBS       = -L/usr/local/petsc/2.3.2-p6/lib/linux_mpi_cxx_prom/ -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc 


include ${PETSC_DIR}/bmake/common/base

.cc.o :
	$(CXX) -c $(CXXFLAGS) ${CCPPFLAGS} -o $@ $< 

hiqlab.cc: qpetsc.o qpetsclua.o qpassemblylua.o

hiqlab: hiqlab.o qpetsclua.o qpetsc.o qpassemblylua.o qpassembly.o $(LIBQLAB)
	$(CXX) -o hiqlab hiqlab.o \
		qpetsclua.o qpetsc.o qpassemblylua.o qpassembly.o $(LIBQLAB) \
		${PETSC_LIBS}  \
		${PETSC_KSP_LIB} \
		${LIBSUPERLU} \
                $(READLINE_LIBS) $(FLTK_LIBS) $(LIBS)

qpetsclua.o: qpetsclua.cc
qpetsc.o: qpetsc.cc
qpassemblylua.o: qpassemblylua.cc
qpassembly.o: qpassembly.cc
hiqlab.o: hiqlab.cc

qpetsclua.cc: qpetsc.pkg
	$(TOLUA) -o qpetsclua.cc -H qpetsclua.h qpetsc.pkg

qpassemblylua.cc: qpassembly.pkg
	$(TOLUA) -o qpassemblylua.cc -H qpassemblylua.h qpassembly.pkg

clean:
	${RM} hiqlab qpetsclua* qpetsclua* qpassemblylua* *.o
