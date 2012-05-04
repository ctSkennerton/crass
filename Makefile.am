AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src man doc
dist_doc_DATA =  man/crass.1 man/crisprtools.1
EXTRA_DIST = doc/manual.tex doc/crispr_spec.tex resources/crass.dtd autogen.sh
if HAVE_PDFLATEX
manual: pdf

doc: pdf

pdf:
		cd doc; make manual;
endif
