rm *.so

gcc -fPIC -shared -fopenmp -o3 -o el3d_openmp.so el3d_openmp.c
