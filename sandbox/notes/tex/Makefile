notes.pdf: notes.tex
	pdflatex notes.tex

all:
	pdflatex notes.tex
	bibtex notes
	pdflatex notes.tex
	pdflatex notes.tex

clean:
	rm -f *~ *.aux *.log *.blg *.bbl notes.pdf

