include ../../../make.inc
include make_trilinos.inc
include make.inc

LIBQLAB_TL= libqlab_trilinos.a

OBJS=	qpassembly_trilinos.o	\
	pareigs.o \
	trilinos_arpack.o \
	trilinos_anasazi.o \
	trilinos_epetra_linearproblem.o \
	trilinos_epetra_matrix.o \
	trilinos_epetra_vector.o \
	Amesos_Operator.o \
	Amesos_Operator_Complex.o \
        AnasaziEpetraMultiVectorComplexAdapter.o \
        AnasaziEpetraOperatorComplexAdapter.o \
        trilinos_amesos_base.o  \
	trilinos_super_matrix.o \
	trilinos_super_vector.o \
	trilinos_komplex.o \
        AztecOO_Operator_Komplex.o \
	trilinos_mesh.o \
	$(USE_MUMPS_OBJS) \
	$(USE_SUPERLU_OBJS)

MUMPS_OBJS=Amesos_Mumps_Complex.o
SUPERLU_OBJS=Amesos_Superludist_Complex.o

QINCL=-I../../corelib -I../../numeric -I.

all: lib
lib: $(LIBQLAB_TL)

$(LIBQLAB_TL): $(OBJS) 
	$(AR) rv $(LIBQLAB_TL) $?
	$(RANLIB) $(LIBQLAB_TL)

.cc.o: 
	$(CXX) -DHAVE_CONFIG_H $(CXXFLAGS) \
		$(QINCL) $(INCLUDE_TL) $(INCLUDE_NP) -c $<

depend:
	makedepend -DHAVE_CONFIG_H $(QINCL) -Y ${OBJS:.o=.cc} \
		Amesos_Mumps_Complex.cc Amesos_Superludist_Complex.cc

clean:
	rm -f *lua.* *.o *~ libqlab_trilinos.a Makefile.bak

# DO NOT DELETE

qpassembly_trilinos.o: ../../corelib/qcomplex.h ../../corelib/qassembly.h
qpassembly_trilinos.o: ../../corelib/qmatrix.h qpassembly_trilinos.h
qpassembly_trilinos.o: trilinos_epetra_vector.h trilinos_epetra_matrix.h
qpassembly_trilinos.o: trilinos_epetra_operator.h ../../corelib/coordmatrix.h
qpassembly_trilinos.o: ../../corelib/cscmatrix.h
qpassembly_trilinos.o: ../../corelib/coordmatrix_g.h ../../corelib/mesh.h
qpassembly_trilinos.o: ../../corelib/element.h
pareigs.o: ../../corelib/qcomplex.h pareigs.h ../../numeric/areigs.h
pareigs.o: ../../numeric/arinfo.cc
trilinos_arpack.o: ../../corelib/qcomplex.h ../../corelib/cscmatrix.h
trilinos_arpack.o: ../../corelib/mesh.h ../../corelib/element.h
trilinos_arpack.o: ../../corelib/qassembly.h ../../corelib/qmatrix.h
trilinos_arpack.o: ../../corelib/coordmatrix.h ../../corelib/coordmatrix_g.h
trilinos_arpack.o: ../../numeric/pmlmode.h pareigs.h ../../numeric/areigs.h
trilinos_arpack.o: qpassembly_trilinos.h trilinos_epetra_vector.h
trilinos_arpack.o: trilinos_epetra_matrix.h trilinos_epetra_operator.h
trilinos_arpack.o: trilinos_arpack.h Amesos_Operator.h
trilinos_arpack.o: Amesos_Operator_Complex.h trilinos_epetra_linearproblem.h
trilinos_arpack.o: trilinos_amesos_base.h
trilinos_anasazi.o: ../../corelib/mesh.h ../../corelib/qcomplex.h
trilinos_anasazi.o: ../../corelib/element.h ../../corelib/qassembly.h
trilinos_anasazi.o: ../../corelib/qmatrix.h ../../corelib/coordmatrix.h
trilinos_anasazi.o: ../../corelib/cscmatrix.h ../../corelib/coordmatrix_g.h
trilinos_anasazi.o: qpassembly_trilinos.h trilinos_epetra_vector.h
trilinos_anasazi.o: trilinos_epetra_matrix.h trilinos_epetra_operator.h
trilinos_anasazi.o: trilinos_anasazi.h Amesos_Operator.h
trilinos_anasazi.o: AnasaziEpetraMultiVectorComplexAdapter.hpp
trilinos_anasazi.o: AnasaziEpetraOperatorComplexAdapter.hpp
trilinos_anasazi.o: Amesos_Operator_Complex.h trilinos_epetra_linearproblem.h
trilinos_anasazi.o: trilinos_amesos_base.h
trilinos_epetra_linearproblem.o: trilinos_epetra_linearproblem.h
trilinos_epetra_linearproblem.o: trilinos_epetra_matrix.h
trilinos_epetra_linearproblem.o: trilinos_epetra_vector.h
trilinos_epetra_linearproblem.o: ../../corelib/qcomplex.h
trilinos_epetra_linearproblem.o: trilinos_epetra_operator.h
trilinos_epetra_matrix.o: trilinos_epetra_matrix.h trilinos_epetra_vector.h
trilinos_epetra_matrix.o: ../../corelib/qcomplex.h trilinos_epetra_operator.h
trilinos_epetra_vector.o: trilinos_epetra_vector.h ../../corelib/qcomplex.h
Amesos_Operator.o: Amesos_Operator.h
Amesos_Operator_Complex.o: Amesos_Operator_Complex.h trilinos_epetra_vector.h
Amesos_Operator_Complex.o: ../../corelib/qcomplex.h trilinos_epetra_matrix.h
Amesos_Operator_Complex.o: trilinos_epetra_operator.h
Amesos_Operator_Complex.o: trilinos_epetra_linearproblem.h
Amesos_Operator_Complex.o: trilinos_amesos_base.h
AnasaziEpetraMultiVectorComplexAdapter.o: AnasaziEpetraMultiVectorComplexAdapter.hpp
AnasaziEpetraMultiVectorComplexAdapter.o: trilinos_epetra_vector.h
AnasaziEpetraMultiVectorComplexAdapter.o: ../../corelib/qcomplex.h
AnasaziEpetraOperatorComplexAdapter.o: AnasaziEpetraOperatorComplexAdapter.hpp
AnasaziEpetraOperatorComplexAdapter.o: trilinos_epetra_operator.h
AnasaziEpetraOperatorComplexAdapter.o: ../../corelib/qcomplex.h
AnasaziEpetraOperatorComplexAdapter.o: AnasaziEpetraMultiVectorComplexAdapter.hpp
AnasaziEpetraOperatorComplexAdapter.o: trilinos_epetra_vector.h
trilinos_amesos_base.o: hiqlab_trilinos_config.h trilinos_amesos_base.h
trilinos_amesos_base.o: trilinos_epetra_linearproblem.h
trilinos_amesos_base.o: trilinos_epetra_matrix.h trilinos_epetra_vector.h
trilinos_amesos_base.o: ../../corelib/qcomplex.h trilinos_epetra_operator.h
trilinos_super_matrix.o: trilinos_super_matrix.h trilinos_epetra_vector.h
trilinos_super_matrix.o: ../../corelib/qcomplex.h trilinos_epetra_operator.h
trilinos_super_matrix.o: trilinos_indexmap.h
trilinos_super_vector.o: trilinos_super_vector.h trilinos_epetra_vector.h
trilinos_super_vector.o: ../../corelib/qcomplex.h trilinos_epetra_operator.h
trilinos_super_vector.o: trilinos_indexmap.h
trilinos_komplex.o: trilinos_komplex.h trilinos_epetra_matrix.h
trilinos_komplex.o: trilinos_epetra_vector.h ../../corelib/qcomplex.h
trilinos_komplex.o: trilinos_epetra_operator.h trilinos_super_matrix.h
trilinos_komplex.o: trilinos_indexmap.h trilinos_super_vector.h
AztecOO_Operator_Komplex.o: AztecOO_Operator_Komplex.h
AztecOO_Operator_Komplex.o: trilinos_epetra_vector.h ../../corelib/qcomplex.h
AztecOO_Operator_Komplex.o: trilinos_epetra_matrix.h
AztecOO_Operator_Komplex.o: trilinos_epetra_operator.h
AztecOO_Operator_Komplex.o: trilinos_epetra_linearproblem.h
AztecOO_Operator_Komplex.o: trilinos_super_matrix.h trilinos_indexmap.h
AztecOO_Operator_Komplex.o: trilinos_super_vector.h trilinos_komplex.h
trilinos_mesh.o: trilinos_mesh.h ../../corelib/mesh.h
trilinos_mesh.o: ../../corelib/qcomplex.h ../../corelib/element.h
trilinos_mesh.o: ../../corelib/qassembly.h ../../corelib/qmatrix.h
trilinos_mesh.o: ../../corelib/coordmatrix.h ../../corelib/cscmatrix.h
trilinos_mesh.o: ../../corelib/coordmatrix_g.h
Amesos_Mumps_Complex.o: Amesos_Mumps_Complex.h ../../corelib/qcomplex.h
Amesos_Mumps_Complex.o: trilinos_epetra_vector.h trilinos_epetra_matrix.h
Amesos_Mumps_Complex.o: trilinos_epetra_operator.h trilinos_amesos_base.h
Amesos_Mumps_Complex.o: trilinos_epetra_linearproblem.h
Amesos_Superludist_Complex.o: Amesos_Superludist_Complex.h
Amesos_Superludist_Complex.o: trilinos_amesos_base.h
Amesos_Superludist_Complex.o: trilinos_epetra_linearproblem.h
Amesos_Superludist_Complex.o: trilinos_epetra_matrix.h
Amesos_Superludist_Complex.o: trilinos_epetra_vector.h
Amesos_Superludist_Complex.o: ../../corelib/qcomplex.h
Amesos_Superludist_Complex.o: trilinos_epetra_operator.h
