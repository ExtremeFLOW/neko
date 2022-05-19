.cu.o:
	$(NVCC) -arch sm_60 -I@top_builddir@/src -I@top_srcdir@/src -O3 -o $@ -c $<
