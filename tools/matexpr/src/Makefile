include ../make.inc

# === Primary targets ===

OBJS= \
	matexpr.o lex.yy.o matexpr-gen.o \
	matexpr-ir.o \
	matexpr-typecheck.o \
	matexpr-ast2const.o \
	matexpr-ast2ir.o \
	matexpr-ir2c.o

.cc.o:
	$(CXX) -c $*.cc

../matexpr: $(OBJS)
	$(CXX) -o ../matexpr $(OBJS)

matexpr.cc: matexpr.y
	$(BISON) -d -v matexpr.y -o matexpr.cc

lex.yy.o: lex.yy.c matexpr.cc
	$(CC) -c lex.yy.c

lex.yy.c: matexpr.l
	$(FLEX) matexpr.l

matexpr.o: matexpr.cc matexpr-ast.h
matexpr-ir.o: matexpr-ir.cc matexpr-ir.h
matexpr-typecheck.o: matexpr-typecheck.cc matexpr-typecheck.h \
  matexpr-ast.h
matexpr-ast2const: matexpr-ast.h
matexpr-ast2ir.o: matexpr-ast2ir.cc matexpr-ast2ir.h matexpr-ast.h \
  matexpr-ir.h
matexpr-ir2c.o: matexpr-ir2c.cc matexpr-ir2c.h matexpr-ir.h
matexpr-gen.o: matexpr-gen.cc matexpr-ast.h matexpr-typecheck.h \
  matexpr-ast2ir.h matexpr-ir.h matexpr-ir2c.h


# === Test targets ===

testx: matexpr test.in
	./matexpr test.in > test.cc
	$(CXX) -o testx test.cc

# === Clean-up targets ===

clean:
	rm -f matexpr.output
	rm -f matexpr *.o *~
	rm -f test.cc testx

realclean: clean
	rm -f lex.yy.c matexpr.cc matexpr.hh 

