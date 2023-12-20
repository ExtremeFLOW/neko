# See INPUT and FILE_PATTERNS in Doxyfile.in
INPUT := $(shell find @top_srcdir@/doc/pages @top_srcdir@/src \
	@top_srcdir@/contrib/rea2nbin \
	! -name '.*' -a '(' -iname '*.f90' -o -name '*.c' -o -name '*.h' -o \
	-name '*.hip' -o -name '*.cu' -o -name '*.md' ')' )

doxygen.stamp: Doxyfile $(INPUT)
	$(DOXYGEN) $^
	date > $@

CLEANFILES = doxygen.stamp

clean-local:
	rm -fr html 

html-local: doxygen.stamp

EXTRA_DIST = pages/code-style.md\
	     pages/publication-list.md\
	     pages/accelerators.md\
	     pages/case-file.md\
	     pages/statistics.md\
             pages/simcomps.md\
             pages/point-zones.md\
             pages/testing.md\
             pages/appendix.md\
             pages/installation.md\
             pages/dev_patterns.md\
	     pages/index.md
