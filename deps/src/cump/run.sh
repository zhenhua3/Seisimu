rm *.so

nvcc --ptxas-options=-v --compiler-options '-fPIC' -o el3d_cuda.so --shared el3d_cuda.cu

gcc -fPIC -shared -fopenmp -o3 -o el3d_openmp.so el3d_openmp.c
