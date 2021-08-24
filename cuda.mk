.cu.o:
	$(NVCC) -O3 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_61,code=sm_61 -o $@ -c $<
