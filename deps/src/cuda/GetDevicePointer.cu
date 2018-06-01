void acGetDevicePointer(
  double *Host_vx, double *Host_vy, double *Host_vz,
  double *Host_tpp, double *Host_rho, double *Host_lambda,
  double *Dev_vx, double *Dev_vy, double *Dev_vz,
  double *Dev_tpp, double *Dev_rho, double *Dev_lambda)
  {
    int err = cudaHostGetDevicePointer((void **)&Dev_vx, (void *)Host_vx,0);
    if(err!= cudaSuccess){printf("%s\n","error 9");}
    err = cudaHostGetDevicePointer((void **)&Dev_vy, (void *)Host_vy,0);
    if(err!= cudaSuccess){printf("%s\n","error 10");}
    err = cudaHostGetDevicePointer((void **)&Dev_vz, (void *)Host_vz,0);
    if(err!= cudaSuccess){printf("%s\n","error 11");}
    err = cudaHostGetDevicePointer((void **)&Dev_tpp, (void *)Host_tpp,0);
    if(err!= cudaSuccess){printf("%s\n","error 12");}
    err = cudaHostGetDevicePointer((void **)&Dev_rho, (void *)Host_rho,0);
    if(err!= cudaSuccess){printf("%s\n","error 13");}
    err = cudaHostGetDevicePointer((void **)&Dev_lambda, (void *)Host_lambda,0);
    if(err!= cudaSuccess){printf("%s\n","error 14");}
  }
