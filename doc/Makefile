# Makefile
# GNU makefile for MMSP documentation
# Questions/comments to gruberja@gmail.com (Jason Gruber)

chapters = MMSP.chapter1.tex \
           MMSP.chapter2.tex \
           MMSP.chapter3.tex \
           MMSP.chapter4.tex \
           MMSP.chapter4.tex

manual: MMSP.manual.tex $(chapters)
	pdflatex MMSP.manual
	pdflatex MMSP.manual

clean:
	rm -f MMSP.pdf *.aux *.log *.toc
