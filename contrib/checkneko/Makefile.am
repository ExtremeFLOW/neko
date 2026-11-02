# checkneko has its own configure script, which Neko's configure runs. It is
# built with the plain Fortran compiler FC and needs only json-fortran. It
# reads the JSON schemas of the case file at run time, from
# share/neko/schemas next to its bin directory, where they are installed
# with it, see README.md for compiling the tool on its own.
ACLOCAL_AMFLAGS = -I ../../m4

if BUILD_CHECKNEKO
bin_PROGRAMS = checkneko
checkneko_SOURCES = checkneko.F90
checkneko_LDADD = $(JSON_Fortran_LIBS)
AM_FCFLAGS = $(JSON_Fortran_CFLAGS)

man_MANS = checkneko.1

# The schemas of doc/schemas, with their subdirectories
schemadir = $(datadir)/neko/schemas
schema_src = $(top_srcdir)/../../doc/schemas

install-data-local:
	@list=`cd $(schema_src) && find . -name '*.schema.json' | sed 's|^\./||'`; \
	for f in $$list; do \
	  d=`dirname "$$f"`; \
	  $(MKDIR_P) "$(DESTDIR)$(schemadir)/$$d" || exit 1; \
	  echo " $(INSTALL_DATA) $(schema_src)/$$f '$(DESTDIR)$(schemadir)/$$d'"; \
	  $(INSTALL_DATA) "$(schema_src)/$$f" "$(DESTDIR)$(schemadir)/$$d" \
	    || exit 1; \
	done

uninstall-local:
	@list=`cd $(schema_src) && find . -name '*.schema.json' | sed 's|^\./||'`; \
	for f in $$list; do \
	  echo " rm -f '$(DESTDIR)$(schemadir)/$$f'"; \
	  rm -f "$(DESTDIR)$(schemadir)/$$f"; \
	done
endif

EXTRA_DIST = README.md

CLEANFILES = checkneko_util.mod checkneko_regex.mod checkneko_schema.mod \
	     checkneko_checks.mod
