.cu.o:
	$(NVCC) -I@top_builddir@/src -I@top_srcdir@/src -O3 -o $@ -c $<
