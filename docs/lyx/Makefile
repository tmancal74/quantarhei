all:

html:
	latex quantarhei_theory.tex
	latex quantarhei_theory.tex
	latex2html -local_icons quantarhei_theory.tex

pdf: 
	latex quantarhei_theory.tex
	latex quantarhei_theory.tex
	dvips -o quantarhei_theory.ps quantarhei_theory.dvi
	ps2pdf quantarhei_theory.ps

clean:
	rm -r -f *.ps *.toc *.dvi *.aux *.log *.idx

delete: clean
	rm -r -f quantarhei_theory/
	rm -r -f quantarhei_theory.pdf

