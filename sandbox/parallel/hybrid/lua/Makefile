include ../../../make.inc
include $(PETSC_DIR)/bmake/$(PETSC_ARCH)/petscconf
include ../../../make.inc
include ../../petsc/make.inc
include ../../trilinos/corelib/make_trilinos.inc

LIBQLAB       = $(LIBQLAB)
LIBQLAB_TL    = $(LIBDIR)/libqlab_trilinos.a
LIBQLAB_PETSC = $(LIBDIR)/libqlab_petsc.a
LIBQLAB_HB    = $(LIBDIR)/libqlab_hybrid.a

OBJS= qpetsc_trilinoslua.o	

QINCL= -I. -I../../corelib -I../../lua \
	-I../../trilinos/corelib \
	-I../../petsc/corelib \
	-I../../hybrid/corelib \
	-I../../trilinos/lua \
	-I../../petsc/lua 

all: lib exe
lib: $(LIBQLAB_HB)
exe: hiqlab

$(LIBQLAB_HB): $(OBJS)
	$(AR) rv $(LIBQLAB_HB) $?
	$(RANLIB) $(LIBQLAB_HB)

hiqlab: hiqlab.o $(LIBQLAB_HB) $(LIBQLAB_TL) $(LIBQLAB_PETSC) $(LIBQLAB)
	$(CXX) -o hiqlab hiqlab.o \
		$(LIBQLAB_HB) $(LIBQLAB_TL) $(LIBQLAB_PETSC) \
                $(LIBQLAB) $(LIB_TL) \
		${SLEPC_LIBS} ${PETSC_LIBS} ${PETSC_KSP_LIB} \
		${HYPRE_LIBS} \
		${PROMETHEUS_LIBS} ${METIS_LIBS} \
		${SUPERLU_LIBS} \
		${X_LIBS} \
		$(LIBS) $(F90LIBS) $(READLINE_LIBS) $(FLTK_LIBS)

hiqlab.o: hiqlab.cc $(OBJS)
	$(CXX) $(CXXFLAGS) $(FLTK_CXXFLAGS) $(QINCL) \
		 $(SLEPC_INCLUDE) $(PETSC_INCLUDE) $(INCLUDE_TL) -c hiqlab.cc

.cc.o: 
	$(CXX) -DHAVE_CONFIG_H $(CXXFLAGS) $(QINCL) \
		$(INCLUDE_TL) $(PETSC_INCLUDE) -c $<

qpetsc_trilinoslua.cc: qpetsc_trilinos.pkg
	$(TOLUA) -o qpetsc_trilinoslua.cc \
		-H qpetsc_trilinoslua.h qpetsc_trilinos.pkg

depend: ${OBJS:.o=.cc}
	makedepend -DHAVE_CONFIG_H $(QINCL) -Y ${OBJS:.o=.cc}

clean:
	rm -f *lua.* *.o hiqlab *~ Makefile.bak

# DO NOT DELETE

qpetsc_trilinoslua.o: ../../hybrid/corelib/PetscOO.h ../../corelib/mesh.h
qpetsc_trilinoslua.o: ../../corelib/qcomplex.h ../../corelib/element.h
qpetsc_trilinoslua.o: ../../corelib/qassembly.h ../../corelib/qmatrix.h
qpetsc_trilinoslua.o: ../../corelib/coordmatrix.h ../../corelib/cscmatrix.h
qpetsc_trilinoslua.o: ../../corelib/coordmatrix_g.h
qpetsc_trilinoslua.o: ../../petsc/corelib/qpetsc_pc.h
qpetsc_trilinoslua.o: ../../petsc/corelib/hiqlab_petsc_config.h
qpetsc_trilinoslua.o: ../../hybrid/corelib/Petsc_Operator.h
qpetsc_trilinoslua.o: ../../trilinos/corelib/trilinos_epetra_vector.h
qpetsc_trilinoslua.o: ../../trilinos/corelib/trilinos_epetra_matrix.h
qpetsc_trilinoslua.o: ../../trilinos/corelib/trilinos_epetra_operator.h
qpetsc_trilinoslua.o: ../../trilinos/corelib/trilinos_epetra_linearproblem.h
qpetsc_trilinoslua.o: ../../trilinos/corelib/trilinos_super_matrix.h
qpetsc_trilinoslua.o: ../../trilinos/corelib/trilinos_indexmap.h
qpetsc_trilinoslua.o: ../../trilinos/corelib/trilinos_super_vector.h
