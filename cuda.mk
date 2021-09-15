.cu.o:
	$(NVCC) -O3 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_80,code=sm_80 -o $@ -c $<
