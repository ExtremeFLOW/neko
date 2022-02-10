# See INPUT and FILE_PATTERNS in Doxyfile.in
INPUT := $(shell find @top_srcdir@/doc/pages @top_srcdir@/src \
	@top_srcdir@/contrib/rea2nbin \
	-iname '*.f90' -or -name '*.c' -or -name '*.h' -or -name '*.hip' \
	-or -name '*.cu' -or -name '*.md')

doxygen.stamp: Doxyfile $(INPUT)
	$(DOXYGEN) $^
	date > $@

CLEANFILES = doxygen.stamp

clean-local:
	rm -fr html 

html-local: doxygen.stamp
