extern "C"{
  #include<stdio.h>
  #include<unistd.h>
  #include<sys/stat.h>
  #include<sys/mman.h>
  #include<fcntl.h>
  #include<string.h>
  #include<stdlib.h>
  #include"mmap_snapshot.c"
  #include"HostRegister.cu"
  #include"GetDevicePointer.cu"
  #include"ackernel.cu"

  //double *fdc/
  //double *intvl: dz, dx, dt
  //double *modelsize : BDnDZ, BDnHX, nDZ, nHX
  void ac3d_cuda(double *vx, int BD_nx_vx, int BD_ny_vx, int BD_nz_vx,
  double *pvxbtpp,
  double *vy, int BD_nx_vy, int BD_ny_vy, int BD_nz_vy,
  double *pvybtpp,
  double *vz, int BD_nx_vz, int BD_ny_vz, int BD_nz_vz,
  double *pvzbtpp,
  double *tpp, double *ptppbvx, double *ptppbvy, double *ptppbvz,
  int BD_nx_tpp, int BD_ny_tpp, int BD_nz_tpp,
  double *rho, double *lambda, double *fdc,
  int nT, double dt, double dx, double dy, double dz, int ext,
  double *bhalf, double *ahalf, double *bfull, double *afull,
  char *snp_path, long long dim, int intvl,
  long long *threadim, long long *blockdim)
  {

    int t,nsnp=0;

    double *snp_ptr = mmap_snapshot(snp_path, dim);

    // enable zero copy access at device end
    cudaError_t err = cudaSetDeviceFlags(cudaDeviceMapHost);
    if(err!= cudaSuccess){printf("%s\n",cudaGetErrorString(err));}

    // copy finite difference coefficient to constant memory
    err = cudaMemcpyToSymbol(dev_fdc,fdc,4*sizeof(double));
    if(err!= cudaSuccess){printf("%s\n",cudaGetErrorString(err));}

    // number of threads and blocks used
    dim3 threads(*threadim,*(threadim+1),*(threadim+2));
    dim3 blocks(*blockdim,*(blockdim+1),*(blockdim+2));

    // pinned host memory for faster transfer between host and device mem
    acHostRegister(vx, BD_nx_vx, BD_ny_vx, BD_nz_vx,
    vy, BD_nx_vy, BD_ny_vy, BD_nz_vy,
    vz, BD_nx_vz, BD_ny_vz, BD_nz_vz,
    tpp, BD_nx_tpp, BD_ny_tpp, BD_nz_tpp,
    rho, lambda);

    // device
    double *dev_vx, *dev_vy, *dev_vz;
    double *dev_tpp;
    double *dev_rho, *dev_lambda;

    //get device pointer
    // acGetDevicePointer(
    //   vx, vy, vz, tpp, rho, lambda,
    //   dev_vx, dev_vy, dev_vz, dev_tpp, dev_rho, dev_lambda);
      err = cudaHostGetDevicePointer((void **)&dev_vx, (void *)vx,0);
      if(err!= cudaSuccess){printf("%s\n",cudaGetErrorString(err));}
      err = cudaHostGetDevicePointer((void **)&dev_vy, (void *)vy,0);
      if(err!= cudaSuccess){printf("%s\n",cudaGetErrorString(err));}
      err = cudaHostGetDevicePointer((void **)&dev_vz, (void *)vz,0);
      if(err!= cudaSuccess){printf("%s\n",cudaGetErrorString(err));}
      err = cudaHostGetDevicePointer((void **)&dev_tpp, (void *)tpp,0);
      if(err!= cudaSuccess){printf("%s\n",cudaGetErrorString(err));}
      err = cudaHostGetDevicePointer((void **)&dev_rho, (void *)rho,0);
      if(err!= cudaSuccess){printf("%s\n",cudaGetErrorString(err));}
      err = cudaHostGetDevicePointer((void **)&dev_lambda, (void *)lambda,0);
      if(err!= cudaSuccess){printf("%s\n",cudaGetErrorString(err));}
    //************ time iteration ************//
    for(t=0;t<nT;t++)
    {

      // particle velocities
      ackernel_v<<<blocks,threads>>>
      (dev_vx, BD_nx_vx, BD_ny_vx, BD_nz_vx,
        dev_vy, BD_nx_vy, BD_ny_vy, BD_nz_vy,
        dev_vz, BD_nx_vz, BD_ny_vz, BD_nz_vz,
        dev_tpp, BD_nx_tpp, BD_ny_tpp, BD_nz_tpp,
        dev_rho, dx, dy, dz, dt);

      cudaDeviceSynchronize();

      // stress
      ackernel_tau<<<blocks,threads>>>
      (dev_vx, BD_nx_vx, BD_ny_vx, BD_nz_vx,
        dev_vy, BD_nx_vy, BD_ny_vy, BD_nz_vy,
        dev_vz, BD_nx_vz, BD_ny_vz, BD_nz_vz,
        dev_tpp, BD_nx_tpp, BD_ny_tpp, BD_nz_tpp,
        dev_lambda, dx, dy, dz, dt);

      cudaDeviceSynchronize();

      if(t%intvl==0){
        memcpy(snp_ptr+nsnp*BD_nx_vz*BD_ny_vz*BD_nz_vz,vz,BD_nx_vz*BD_ny_vz*BD_nz_vz*sizeof(double));
        nsnp++;
      }
    }

    acHostUnRegister(vx, vy, vz, tpp, rho, lambda);
    err = cudaDeviceReset();
    if(err!= cudaSuccess){printf("%s\n",cudaGetErrorString(err));}
    munmap_snapshot(snp_ptr, snp_path);
  }
}
