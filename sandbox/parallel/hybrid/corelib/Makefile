include ../../../make.inc
include $(PETSC_DIR)/bmake/$(PETSC_ARCH)/petscconf
include ../../../make.inc
include ../../petsc/make.inc
include ../../trilinos/corelib/make_trilinos.inc

LIBQLAB_HB= $(LIBDIR)/libqlab_hybrid.a

OBJS=	Petsc_Operator.o \
	PetscOO.o	

QINCL= -I../../corelib -I../../trilinos/corelib -I../../petsc/corelib -I. 

all: lib
lib: $(LIBQLAB_HB)

$(LIBQLAB_HB): $(OBJS) 
	$(AR) rv $(LIBQLAB_HB) $?
	$(RANLIB) $(LIBQLAB_HB)

.cc.o: 
	$(CXX) -DHAVE_CONFIG_H $(CXXFLAGS) $(QINCL) \
		$(INCLUDE_TL) $(INCLUDE_NP) $(PETSC_INCLUDE) -c $<

depend:
	makedepend -DHAVE_CONFIG_H $(QINCL) ${OBJS:.o=.cc}

clean:
	rm -f *lua.* *.o *~ libqlab_hybrid.a Makefile.bak

# DO NOT DELETE

Petsc_Operator.o: Petsc_Operator.h
Petsc_Operator.o: ../../trilinos/corelib/trilinos_epetra_vector.h
Petsc_Operator.o: ../../corelib/qcomplex.h
Petsc_Operator.o: ../../trilinos/corelib/trilinos_epetra_matrix.h
Petsc_Operator.o: ../../trilinos/corelib/trilinos_epetra_operator.h
Petsc_Operator.o: ../../trilinos/corelib/trilinos_epetra_linearproblem.h
Petsc_Operator.o: ../../trilinos/corelib/trilinos_super_matrix.h
Petsc_Operator.o: ../../trilinos/corelib/trilinos_indexmap.h
Petsc_Operator.o: ../../trilinos/corelib/trilinos_super_vector.h
Petsc_Operator.o: ../../corelib/mesh.h ../../corelib/element.h
Petsc_Operator.o: ../../corelib/qassembly.h ../../corelib/qmatrix.h
Petsc_Operator.o: ../../corelib/coordmatrix.h ../../corelib/cscmatrix.h
Petsc_Operator.o: ../../corelib/coordmatrix_g.h
Petsc_Operator.o: ../../petsc/corelib/qpetsc_pc.h
Petsc_Operator.o: ../../petsc/corelib/hiqlab_petsc_config.h
PetscOO.o: PetscOO.h ../../corelib/mesh.h ../../corelib/qcomplex.h
PetscOO.o: ../../corelib/element.h ../../corelib/qassembly.h
PetscOO.o: ../../corelib/qmatrix.h ../../corelib/coordmatrix.h
PetscOO.o: ../../corelib/cscmatrix.h ../../corelib/coordmatrix_g.h
PetscOO.o: ../../petsc/corelib/qpetsc_pc.h
PetscOO.o: ../../petsc/corelib/hiqlab_petsc_config.h
