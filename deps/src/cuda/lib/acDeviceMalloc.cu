__host__ void acDeviceMalloc(double *dev_vx[], int BD_nx_vx, int BD_nz_vx,
  double *dev_vy[], int BD_nx_vy, int BD_nz_vy,
  double *dev_vz[], int BD_nx_vz, int BD_nz_vz,
  double *dev_tpp[], int BD_nx_tpp, int BD_nz_tpp,
  double *dev_rho[], double *dev_lambda[],
  int AssignedStreamNum, int RegStreamDim)
  {
    int nstream, i;
    cudaError_t err[6];
    for(nstream=0;nstream<AssignedStreamNum;nstream++){
      err[0] = cudaMalloc((void**)&dev_vx[nstream], (RegStreamDim+3)*BD_nx_vx*BD_nz_vx*sizeof(double));
      err[1] = cudaMalloc((void**)&dev_vy[nstream], (RegStreamDim+3)*BD_nx_vy*BD_nz_vy*sizeof(double));
      err[2] = cudaMalloc((void**)&dev_vz[nstream], (RegStreamDim+3)*BD_nx_vz*BD_nz_vz*sizeof(double));
      err[3] = cudaMalloc((void**)&dev_tpp[nstream], (RegStreamDim+3)*BD_nx_tpp*BD_nz_tpp*sizeof(double));
      err[4] = cudaMalloc((void**)&dev_lambda[nstream], (RegStreamDim+3)*BD_nx_tpp*BD_nz_tpp*sizeof(double));
      err[5] = cudaMalloc((void**)&dev_rho[nstream], (RegStreamDim+3)*BD_nx_tpp*BD_nz_tpp*sizeof(double));
      for(i=0;i<6;i++)
      {
        if(err[i]!= cudaSuccess)
        {printf("Device Memory Allocation Error No. %i: %s\n",i, cudaGetErrorString(err[i]));}
      }
    }
  }
