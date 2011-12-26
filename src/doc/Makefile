all: manual.force pmanual.force
	(cd umanual; make)
	(cd emanual; make)

manual.pdf: manual.tex
	pdflatex manual.tex

pmanual.pdf: pmanual.tex
	pdflatex pmanual.tex

pmanual.force:
	pdflatex pmanual.tex
	#bibtex pmanual
	pdflatex pmanual.tex
	pdflatex pmanual.tex

manual.force:
	pdflatex manual.tex
	#bibtex manual
	pdflatex manual.tex
	pdflatex manual.tex

clean:
	rm -f *.log *.aux *.bbl *.blg *.toc 
	(cd umanual; make clean)
	(cd emanual; make clean)

cleanpdf:
	rm -f *.pdf *.ps
	(cd umanual; make remove_pdf)
	(cd emanual; make remove_pdf)
