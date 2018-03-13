rm *.so

gcc -fPIC -shared -fopenmp -o3 -o ac2d_openmp.so ac2d_openmp.c
gcc -fPIC -shared -fopenmp -o3 -o el2d_openmp.so el2d_openmp.c
gcc -fPIC -shared -fopenmp -o3 -o ac3d_openmp.so ac3d_openmp.c
gcc -fPIC -shared -fopenmp -o3 -o el3d_openmp.so el3d_openmp.c
