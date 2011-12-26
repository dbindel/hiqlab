# HiQLab
# Copyright (c): Regents of the University of California

include ../../make.inc

OBJS=	material_model.o \
	couple_em2d.matexp.o	\
	elastic2d.o     \
	pmlfunc.o	\
	pmlelement.o	\
	pmlscalar1d.o	\
	pmlscalar2d.o	\
	pmlscalar3d.o	\
	pmlscalarax.o	\
	pmlelastic2d.matexp.o	\
	pmlelastic3d.matexp.o	\
	pmlelasticax.matexp.o	\
	pmlelastictax.o	\
        pmlelastic2d_te.matexp.o \
        pmlelasticax_te.matexp.o \
        pmlelastic3d_te.matexp.o \
        pmlelastic3d_pz.matexp.o \
        pmlelastic2d_pz.matexp.o \
        pmlelastic2hd_pz.matexp.o \
        tie_field.o \
        tie_field2.o \
        model_circuit.o

QINCL=-I../corelib

all: lib
lib: $(LIBQLAB)

$(LIBQLAB): $(OBJS)
	$(AR) rv $(LIBQLAB) $?
	$(RANLIB) $(LIBQLAB)

.cc.o: 
	$(CXX) $(CXXFLAGS) $(QINCL) -c $<

%.matexp.cc: %.cc
	$(MATEXPR) -line $< > $@

clean:
	rm -f *.o *~ Makefile.bak *.matexp.cc

depend:
	makedepend $(QINCL) -Y ${OBJS:.o=.cc}

# DO NOT DELETE

material_model.o: material_model.h
couple_em2d.o: ../corelib/qmatrix.h ../corelib/gaussquad.h
couple_em2d.o: ../corelib/qcomplex.h ../corelib/shapes.h couple_em2d.h
couple_em2d.o: pmlelement.h pmlfunc.h ../corelib/element.h
couple_em2d.o: ../corelib/qassembly.h ../corelib/fieldeval.h
couple_em2d.o: ../corelib/mesh.h
elastic2d.o: ../corelib/qmatrix.h ../corelib/gaussquad.h
elastic2d.o: ../corelib/qcomplex.h ../corelib/shapes.h elastic2d.h
elastic2d.o: ../corelib/element.h ../corelib/qassembly.h
elastic2d.o: ../corelib/fieldeval.h ../corelib/mesh.h material_model.h
pmlfunc.o: ../corelib/luasupport.h pmlfunc.h
pmlelement.o: ../corelib/luasupport.h ../corelib/qmatrix.h pmlelement.h
pmlelement.o: pmlfunc.h ../corelib/element.h ../corelib/qassembly.h
pmlelement.o: ../corelib/qcomplex.h ../corelib/fieldeval.h ../corelib/mesh.h
pmlscalar1d.o: ../corelib/qmatrix.h ../corelib/gaussquad.h
pmlscalar1d.o: ../corelib/qcomplex.h ../corelib/shapes.h pmlscalar1d.h
pmlscalar1d.o: pmlelement.h pmlfunc.h ../corelib/element.h
pmlscalar1d.o: ../corelib/qassembly.h ../corelib/fieldeval.h
pmlscalar1d.o: ../corelib/mesh.h
pmlscalar2d.o: ../corelib/qmatrix.h ../corelib/gaussquad.h
pmlscalar2d.o: ../corelib/qcomplex.h ../corelib/shapes.h pmlscalar2d.h
pmlscalar2d.o: pmlelement.h pmlfunc.h ../corelib/element.h
pmlscalar2d.o: ../corelib/qassembly.h ../corelib/fieldeval.h
pmlscalar2d.o: ../corelib/mesh.h
pmlscalar3d.o: ../corelib/qmatrix.h ../corelib/gaussquad.h
pmlscalar3d.o: ../corelib/qcomplex.h ../corelib/shapes.h pmlscalar3d.h
pmlscalar3d.o: pmlelement.h pmlfunc.h ../corelib/element.h
pmlscalar3d.o: ../corelib/qassembly.h ../corelib/fieldeval.h
pmlscalar3d.o: ../corelib/mesh.h
pmlscalarax.o: ../corelib/qmatrix.h ../corelib/gaussquad.h
pmlscalarax.o: ../corelib/qcomplex.h ../corelib/shapes.h pmlscalarax.h
pmlscalarax.o: pmlelement.h pmlfunc.h ../corelib/element.h
pmlscalarax.o: ../corelib/qassembly.h ../corelib/fieldeval.h
pmlscalarax.o: ../corelib/mesh.h
pmlelastic3d.o: ../corelib/qmatrix.h ../corelib/gaussquad.h
pmlelastic3d.o: ../corelib/qcomplex.h ../corelib/shapes.h pmlelastic3d.h
pmlelastic3d.o: pmlelement.h pmlfunc.h ../corelib/element.h
pmlelastic3d.o: ../corelib/qassembly.h ../corelib/fieldeval.h
pmlelastic3d.o: ../corelib/mesh.h material_model.h
pmlelastictax.o: ../corelib/qmatrix.h ../corelib/gaussquad.h
pmlelastictax.o: ../corelib/qcomplex.h ../corelib/shapes.h pmlelastictax.h
pmlelastictax.o: pmlelement.h pmlfunc.h ../corelib/element.h
pmlelastictax.o: ../corelib/qassembly.h ../corelib/fieldeval.h
pmlelastictax.o: ../corelib/mesh.h
pmlelastic2d_te.o: ../corelib/qmatrix.h ../corelib/gaussquad.h
pmlelastic2d_te.o: ../corelib/qcomplex.h ../corelib/shapes.h
pmlelastic2d_te.o: pmlelastic2d_te.h pmlelement.h pmlfunc.h
pmlelastic2d_te.o: ../corelib/element.h ../corelib/qassembly.h
pmlelastic2d_te.o: ../corelib/fieldeval.h material_model.h ../corelib/mesh.h
pmlelasticax_te.o: ../corelib/qmatrix.h ../corelib/gaussquad.h
pmlelasticax_te.o: ../corelib/qcomplex.h ../corelib/shapes.h
pmlelasticax_te.o: pmlelasticax_te.h pmlelement.h pmlfunc.h
pmlelasticax_te.o: ../corelib/element.h ../corelib/qassembly.h
pmlelasticax_te.o: ../corelib/fieldeval.h ../corelib/mesh.h material_model.h
pmlelastic3d_te.o: ../corelib/qmatrix.h ../corelib/gaussquad.h
pmlelastic3d_te.o: ../corelib/qcomplex.h ../corelib/shapes.h
pmlelastic3d_te.o: pmlelastic3d_te.h pmlelement.h pmlfunc.h
pmlelastic3d_te.o: ../corelib/element.h ../corelib/qassembly.h
pmlelastic3d_te.o: ../corelib/fieldeval.h ../corelib/mesh.h material_model.h
pmlelastic3d_pz.o: ../corelib/qmatrix.h ../corelib/gaussquad.h
pmlelastic3d_pz.o: ../corelib/qcomplex.h ../corelib/shapes.h
pmlelastic3d_pz.o: pmlelastic3d_pz.h pmlelement.h pmlfunc.h
pmlelastic3d_pz.o: ../corelib/element.h ../corelib/qassembly.h
pmlelastic3d_pz.o: ../corelib/fieldeval.h ../corelib/mesh.h material_model.h
pmlelastic2d_pz.o: ../corelib/qmatrix.h ../corelib/gaussquad.h
pmlelastic2d_pz.o: ../corelib/qcomplex.h ../corelib/shapes.h
pmlelastic2d_pz.o: pmlelastic2d_pz.h pmlelement.h pmlfunc.h
pmlelastic2d_pz.o: ../corelib/element.h ../corelib/qassembly.h
pmlelastic2d_pz.o: ../corelib/fieldeval.h material_model.h ../corelib/mesh.h
pmlelastic2hd_pz.o: ../corelib/qmatrix.h ../corelib/gaussquad.h
pmlelastic2hd_pz.o: ../corelib/qcomplex.h ../corelib/shapes.h
pmlelastic2hd_pz.o: pmlelastic2hd_pz.h pmlelement.h pmlfunc.h
pmlelastic2hd_pz.o: ../corelib/element.h ../corelib/qassembly.h
pmlelastic2hd_pz.o: ../corelib/fieldeval.h material_model.h ../corelib/mesh.h
tie_field.o: ../corelib/luasupport.h ../corelib/qmatrix.h tie_field.h
tie_field.o: ../corelib/element.h ../corelib/qassembly.h
tie_field.o: ../corelib/qcomplex.h ../corelib/fieldeval.h ../corelib/mesh.h
tie_field2.o: ../corelib/luasupport.h ../corelib/qmatrix.h tie_field2.h
tie_field2.o: ../corelib/element.h ../corelib/qassembly.h
tie_field2.o: ../corelib/qcomplex.h ../corelib/fieldeval.h ../corelib/mesh.h
model_circuit.o: ../corelib/qmatrix.h ../corelib/qcomplex.h model_circuit.h
model_circuit.o: ../corelib/element.h ../corelib/qassembly.h
model_circuit.o: ../corelib/fieldeval.h ../corelib/mesh.h
