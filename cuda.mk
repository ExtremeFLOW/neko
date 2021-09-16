.cu.o:
	$(NVCC) -I@top_builddir@/src -O3 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_75,code=sm_75  -gencode arch=compute_80,code=sm_80 -o $@ -c $<
