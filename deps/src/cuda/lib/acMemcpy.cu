__host__ void acMemcpyHToDforParticleVel(
  double *dev_vx[], double *vx, int BD_nx_vx, int BD_nz_vx,
  double *dev_vy[], double *vy, int BD_nx_vy, int BD_nz_vy,
  double *dev_vz[], double *vz, int BD_nx_vz, int BD_nz_vz,
  double *dev_tpp[], double *tpp, int BD_nx_tpp, int BD_nz_tpp,
  double *dev_rho[], double *rho,
  int TotalStreamNum, int AssignedStreamNum, cudaStream_t stream[],
  int *vx_PV_start, int *vy_PV_start, int *vz_PV_start,
  int *tpp_PV_start, int *rho_PV_start,
  int *vx_PV_offset, int *vy_PV_offset, int *vz_PV_offset,
  int *tpp_PV_offset, int *rho_PV_offset)
{
  int nstream;
  cudaError_t err;
  for(nstream=0 ; nstream < TotalStreamNum ; nstream++){

    err = cudaMemcpyAsync(dev_vx[nstream%AssignedStreamNum],
      vx+*(vx_PV_start+nstream)*BD_nz_vx*BD_nx_vx,
      *(vx_PV_offset+nstream)*BD_nz_vx*BD_nx_vx*sizeof(double),
      cudaMemcpyHostToDevice, stream[nstream%AssignedStreamNum]);
      if(err!= cudaSuccess){printf("1 %i %s\n",nstream,cudaGetErrorString(err));}

    err = cudaMemcpyAsync(dev_vy[nstream%AssignedStreamNum],
      vy+*(vy_PV_start+nstream)*BD_nz_vy*BD_nx_vy,
      *(vy_PV_offset+nstream)*BD_nz_vy*BD_nx_vy*sizeof(double),
      cudaMemcpyHostToDevice, stream[nstream%AssignedStreamNum]);
      if(err!= cudaSuccess){printf("2 %i %s\n",nstream,cudaGetErrorString(err));}

    err = cudaMemcpyAsync(dev_vz[nstream%AssignedStreamNum],
      vz+*(vz_PV_start+nstream)*BD_nz_vz*BD_nx_vz,
      *(vz_PV_offset+nstream)*BD_nz_vz*BD_nx_vz*sizeof(double),
      cudaMemcpyHostToDevice, stream[nstream%AssignedStreamNum]);
      if(err!= cudaSuccess){printf("3 %i %s\n",nstream,cudaGetErrorString(err));}

    err = cudaMemcpyAsync(dev_rho[nstream%AssignedStreamNum],
      rho+*(rho_PV_start+nstream)*BD_nz_tpp*BD_nx_tpp,
      *(rho_PV_offset+nstream)*BD_nz_tpp*BD_nx_tpp*sizeof(double),
      cudaMemcpyHostToDevice, stream[nstream%AssignedStreamNum]);
      if(err!= cudaSuccess){printf("4 %i %s\n",nstream,cudaGetErrorString(err));}

    err = cudaMemcpyAsync(dev_tpp[nstream%AssignedStreamNum],
      tpp+*(tpp_PV_start+nstream)*BD_nz_tpp*BD_nx_tpp,
      *(tpp_PV_offset+nstream)*BD_nz_tpp*BD_nx_tpp*sizeof(double),
      cudaMemcpyHostToDevice, stream[nstream%AssignedStreamNum]);
      if(err!= cudaSuccess){printf("5 %i %s\n",nstream,cudaGetErrorString(err));}


  }
}

__host__ void acMemcpyDToHforParticleVel(
  double *dev_vx[], double *vx, int BD_nx_vx, int BD_nz_vx,
  double *dev_vy[], double *vy, int BD_nx_vy, int BD_nz_vy,
  double *dev_vz[], double *vz, int BD_nx_vz, int BD_nz_vz,
  int TotalStreamNum, int AssignedStreamNum, cudaStream_t stream[],
  int *vx_PV_start, int *vy_PV_start, int *vz_PV_start,
  int *vx_PV_offset, int *vy_PV_offset, int *vz_PV_offset)
{
  int nstream;
  cudaError_t err;
  for(nstream=0 ; nstream < TotalStreamNum ; nstream++){

    err = cudaMemcpyAsync(vx+*(vx_PV_start+nstream)*BD_nz_vx*BD_nx_vx,
      dev_vx[nstream%AssignedStreamNum],
      *(vx_PV_offset+nstream)*BD_nz_vx*BD_nx_vx*sizeof(double),
      cudaMemcpyDeviceToHost, stream[nstream%AssignedStreamNum]);
      if(err!= cudaSuccess){printf("6 %i %s\n",nstream,cudaGetErrorString(err));}

    err = cudaMemcpyAsync(vy+*(vy_PV_start+nstream)*BD_nz_vy*BD_nx_vy,
      dev_vy[nstream%AssignedStreamNum],
      *(vy_PV_offset+nstream)*BD_nz_vy*BD_nx_vy*sizeof(double),
      cudaMemcpyDeviceToHost, stream[nstream%AssignedStreamNum]);
      if(err!= cudaSuccess){printf("7 %i %s\n",nstream,cudaGetErrorString(err));}


    err = cudaMemcpyAsync(vz+*(vz_PV_start+nstream)*BD_nz_vz*BD_nx_vz,
      dev_vz[nstream%AssignedStreamNum],
      *(vz_PV_offset+nstream)*BD_nz_vz*BD_nx_vz*sizeof(double),
      cudaMemcpyDeviceToHost, stream[nstream%AssignedStreamNum]);
      if(err!= cudaSuccess){printf("8 %i %s\n",nstream,cudaGetErrorString(err));}

  }
}


__host__ void acMemcpyHToDforStress(
  double *dev_vx[], double *vx, int BD_nx_vx, int BD_nz_vx,
  double *dev_vy[], double *vy, int BD_nx_vy, int BD_nz_vy,
  double *dev_vz[], double *vz, int BD_nx_vz, int BD_nz_vz,
  double *dev_tpp[], double *tpp, int BD_nx_tpp, int BD_nz_tpp,
  double *dev_lambda[], double *lambda,
  int TotalStreamNum, int AssignedStreamNum, cudaStream_t stream[],
  int *vx_SS_start, int *vy_SS_start, int *vz_SS_start,
  int *tpp_SS_start, int *lambda_SS_start,
  int *vx_SS_offset, int *vy_SS_offset, int *vz_SS_offset,
  int *tpp_SS_offset, int *lambda_SS_offset)
  {
    int nstream;
    cudaError_t err;
    for(nstream=0 ; nstream < TotalStreamNum ; nstream++){

    err = cudaMemcpyAsync(dev_vx[nstream%AssignedStreamNum],
        vx+*(vx_SS_start+nstream)*BD_nz_vx*BD_nx_vx,
        *(vx_SS_offset+nstream)*BD_nz_vx*BD_nx_vx*sizeof(double),
        cudaMemcpyHostToDevice, stream[nstream%AssignedStreamNum]);
        if(err!= cudaSuccess){printf("9 %s\n",cudaGetErrorString(err));}

    err = cudaMemcpyAsync(dev_vy[nstream%AssignedStreamNum],
        vy+*(vy_SS_start+nstream)*BD_nz_vy*BD_nx_vy,
        *(vy_SS_offset+nstream)*BD_nz_vy*BD_nx_vy*sizeof(double),
        cudaMemcpyHostToDevice, stream[nstream%AssignedStreamNum]);
        if(err!= cudaSuccess){printf("10 %s\n",cudaGetErrorString(err));}

    err = cudaMemcpyAsync(dev_vz[nstream%AssignedStreamNum],
        vz+*(vz_SS_start+nstream)*BD_nz_vz*BD_nx_vz,
        *(vz_SS_offset+nstream)*BD_nz_vz*BD_nx_vz*sizeof(double),
        cudaMemcpyHostToDevice, stream[nstream%AssignedStreamNum]);
        if(err!= cudaSuccess){printf("11 %s\n",cudaGetErrorString(err));}

    err = cudaMemcpyAsync(dev_tpp[nstream%AssignedStreamNum],
        tpp+*(tpp_SS_start+nstream)*BD_nz_tpp*BD_nx_tpp,
        *(tpp_SS_offset+nstream)*BD_nz_tpp*BD_nx_tpp*sizeof(double),
        cudaMemcpyHostToDevice, stream[nstream%AssignedStreamNum]);
        if(err!= cudaSuccess){printf("12 %s\n",cudaGetErrorString(err));}

    err = cudaMemcpyAsync(dev_lambda[nstream%AssignedStreamNum],
        lambda+*(tpp_SS_start+nstream)*BD_nz_tpp*BD_nx_tpp,
        *(lambda_SS_offset+nstream)*BD_nz_tpp*BD_nx_tpp*sizeof(double),
        cudaMemcpyHostToDevice, stream[nstream%AssignedStreamNum]);
        if(err!= cudaSuccess){printf("13 %s\n",cudaGetErrorString(err));}
      }
    }



__host__ void acMemcpyDToHforStress(
  double *dev_tpp[], double *tpp, int BD_nx_tpp, int BD_nz_tpp,
  int TotalStreamNum, int AssignedStreamNum, cudaStream_t stream[],
  int *tpp_SS_start, int *tpp_SS_offset)
  {
    int nstream;
    cudaError_t err;
    for(nstream=0 ; nstream < TotalStreamNum ; nstream++){

    err = cudaMemcpyAsync(tpp+*(tpp_SS_start+nstream)*BD_nz_tpp*BD_nx_tpp,
        dev_tpp[nstream%AssignedStreamNum],
        *(tpp_SS_offset+nstream)*BD_nz_tpp*BD_nx_tpp*sizeof(double),
        cudaMemcpyDeviceToHost, stream[nstream%AssignedStreamNum]);
        if(err!= cudaSuccess){printf("14 %s\n",cudaGetErrorString(err));}
      }
    }
