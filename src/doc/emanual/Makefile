all:
	pdflatex emanual.tex
	pdflatex emanual.tex
	pdflatex emanual.tex

emanual.pdf: emanual.tex
	pdflatex emanual.tex

clean:
	rm -rf *.log *.dvi *.aux *.toc

remove_pdf:
	rm -rf *.pdf
