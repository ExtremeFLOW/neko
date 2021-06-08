.cu.o:
	$(NVCC) -O3 -o $@ -c $<
