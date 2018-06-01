rm *.so

nvcc --ptxas-options=-v --compiler-options '-fPIC' -o ac3d_cuda.so --shared ac3d_cuda.cu
