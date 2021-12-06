.hip.o:
	$(HIPCC) -O3 -I@top_builddir@/src -o $@ -c $<
