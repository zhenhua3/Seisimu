rm /home/lzh/Dropbox/Zhenhua/Ongoing/Seisimu/deps/builds/*openmp.so

gcc -fPIC -shared -fopenmp -o3 -o /home/lzh/Dropbox/Zhenhua/Ongoing/Seisimu/deps/builds/ac2d_openmp.so ac2d_openmp.c
gcc -fPIC -shared -fopenmp -o3 -o /home/lzh/Dropbox/Zhenhua/Ongoing/Seisimu/deps/builds/el2d_openmp.so el2d_openmp.c
gcc -fPIC -shared -fopenmp -o3 -o /home/lzh/Dropbox/Zhenhua/Ongoing/Seisimu/deps/builds/ac3d_openmp.so ac3d_openmp.c
gcc -fPIC -shared -fopenmp -o3 -o /home/lzh/Dropbox/Zhenhua/Ongoing/Seisimu/deps/builds/el3d_openmp.so el3d_openmp.c
