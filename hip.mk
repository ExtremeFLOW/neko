.hip.o:
	$(HIPCC) -I@top_builddir@/src -I@top_srcdir@/src -O3 -o $@ -c $<
