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

EXTRA_DIST = pages/appendices/governing-equations.md \
			pages/appendices/publications.md \
			pages/developer-guide \
			pages/developer-guide/accelerators.md \
			pages/developer-guide/run_time_selectable_types.md \
			pages/developer-guide/important_types.md \
			pages/developer-guide/code-style.md \
			pages/developer-guide/contributing.md \
			pages/developer-guide/dev_patterns.md \
			pages/developer-guide/testing.md \
			pages/user-guide \
			pages/user-guide/case-file.md \
			pages/user-guide/installation.md \
			pages/user-guide/io.md \
			pages/user-guide/point-zones.md \
			pages/user-guide/simcomps.md \
			pages/user-guide/statistics-guide.md \
			pages/user-guide/user-file.md \
			pages/appendices.md \
			pages/developer-guide.md \
			pages/index.md \
			pages/user-guide.md
