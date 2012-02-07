PDFFILE = manual.pdf crispr_spec.pdf
if HAVE_PDFLATEX
manual: pdf

doc: pdf

pdf: $(PDFFILE)

%.pdf: %.tex
		${PDFLATEX} $?; ${PDFLATEX} $?;
endif
