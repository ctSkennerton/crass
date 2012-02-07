
AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src man 
dist_doc_DATA = resources/crass.dtd man/crass.1 doc/manual.tex doc/crispr_spec.tex

if HAVE_PDFLATEX
manual: pdf

doc: pdf

pdf:
		cd doc; make manual;
endif
