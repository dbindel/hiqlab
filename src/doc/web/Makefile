include ../../../make.inc

notes.pdf: notes.tex mesh.tex element.tex fieldeval.tex
	pdflatex notes.tex
	pdflatex notes.tex

force:
	pdflatex notes.tex

mesh.tex: ../../corelib/mesh.h ../../corelib/mesh.cc
	$(DSBWEB) ../../corelib/mesh.h ../../corelib/mesh.cc -o mesh.tex

element.tex: ../../corelib/element.h ../../corelib/element.cc
	$(DSBWEB) ../../corelib/element.h ../../corelib/element.cc \
		-o element.tex

fieldeval.tex: ../../corelib/fieldeval.h ../../corelib/fieldeval.cc
	$(DSBWEB) ../../corelib/fieldeval.h ../../corelib/fieldeval.cc \
		-o fieldeval.tex

clean:
	rm -f notes.aux notes.log notes.toc
	rm -f mesh.tex element.tex fieldeval.tex

realclean: clean
	rm -f notes.pdf
