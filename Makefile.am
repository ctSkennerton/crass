
AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src man 
dist_doc_DATA = resources/crass.dtd man/crass.1 doc/manual.tex doc/crispr_spec.tex


if HAVE_PDFLATEX
manual: pdf

doc: pdf

#pdf: %.pdf
#
#%.pdf: %.tex
#		cd doc; ${PDFLATEX} $<; ${PDFLATEX} $<; mv $@ ${top_srcdir};
pdf: manual.pdf crispr_spec.pdf

manual.pdf: 
		cd doc; ${PDFLATEX} manual.tex; ${PDFLATEX} manual.tex; mv $@ ${top_srcdir};
crispr_spec.pdf: 
		cd doc; ${PDFLATEX} crispr_spec.tex; ${PDFLATEX} crispr_spec.tex; mv $@ ${top_srcdir};
endif
