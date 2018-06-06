__host__ void acDeviceFree(
  double *dev_vx[], double *dev_vy[],
  double *dev_vz[], double *dev_tpp[],
  double *dev_lambda[], double *dev_rho[],
  int AssignedStreamNum, cudaStream_t stream[])
  {
    int nstream;
    for(nstream=0 ; nstream<AssignedStreamNum ; nstream++)
    {
      cudaFree(dev_vx[nstream]);
      cudaFree(dev_vy[nstream]);
      cudaFree(dev_vz[nstream]);
      cudaFree(dev_tpp[nstream]);
      cudaFree(dev_rho[nstream]);
      cudaFree(dev_lambda[nstream]);
      cudaStreamDestroy(stream[nstream]);
    }
  }
