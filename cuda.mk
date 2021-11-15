.cu.o:
	$(NVCC) -I@top_builddir@/src -O3 -o $@ -c $<
