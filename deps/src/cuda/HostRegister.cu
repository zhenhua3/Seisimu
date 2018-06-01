void acHostRegister(double *vx, int BD_nx_vx, int BD_ny_vx, int BD_nz_vx,
double *vy, int BD_nx_vy, int BD_ny_vy, int BD_nz_vy,
double *vz, int BD_nx_vz, int BD_ny_vz, int BD_nz_vz,
double *tpp, int BD_nx_tpp, int BD_ny_tpp, int BD_nz_tpp,
double *rho, double *lambda)
{
  int err = cudaHostRegister(vx, BD_nx_vx*BD_ny_vx*BD_nz_vx*sizeof(double),cudaHostRegisterMapped);
  if(err!= cudaSuccess){printf("%s\n","error 3");}
  err = cudaHostRegister(vy, BD_nx_vy*BD_ny_vy*BD_nz_vy*sizeof(double),cudaHostRegisterMapped);
  if(err!= cudaSuccess){printf("%s\n","error 4");}
  err = cudaHostRegister(vz, BD_nx_vz*BD_ny_vz*BD_nz_vz*sizeof(double),cudaHostRegisterMapped);
  if(err!= cudaSuccess){printf("%s\n","error 5");}
  err = cudaHostRegister(tpp, BD_nx_tpp*BD_ny_tpp*BD_nz_tpp*sizeof(double),cudaHostRegisterMapped);
  if(err!= cudaSuccess){printf("%s\n","error 6");}
  err = cudaHostRegister(rho, BD_nx_tpp*BD_ny_tpp*BD_nz_tpp*sizeof(double),cudaHostRegisterMapped);
  if(err!= cudaSuccess){printf("%s\n","error 7");}
  err = cudaHostRegister(lambda, BD_nx_tpp*BD_ny_tpp*BD_nz_tpp*sizeof(double),cudaHostRegisterMapped);
  if(err!= cudaSuccess){printf("%s\n","error 8");}
}

void acHostUnRegister(double *vx, double *vy, double *vz, double *tpp, double *rho, double *lambda)
{
  cudaHostUnregister(vx);
  cudaHostUnregister(vy);
  cudaHostUnregister(vz);
  cudaHostUnregister(tpp);
  cudaHostUnregister(rho);
  cudaHostUnregister(lambda);
}

// void elHostRegister(double *vx, int BD_nx_vx, int BD_ny_vx, int BD_nz_vx,
// double *vy, int BD_nx_vy, int BD_ny_vy, int BD_nz_vy,
// double *vz, int BD_nx_vz, int BD_ny_vz, int BD_nz_vz,
// double *tpp, int BD_nx_tpp, int BD_ny_tpp, int BD_nz_tpp,
// double *rho, double *lambda)
// {
//   // cudaHostRegister(vx, BD_nx_vx*BD_ny_vx*BD_nz_vx*sizeof(double),cudaHostRegisterMapped);
//   // cudaHostRegister(vy, BD_nx_vy*BD_ny_vy*BD_nz_vy*sizeof(double),cudaHostRegisterMapped);
//   // cudaHostRegister(vz, BD_nx_vz*BD_ny_vz*BD_nz_vz*sizeof(double),cudaHostRegisterMapped);
//   // cudaHostRegister(tpp, BD_nx_tpp*BD_ny_tpp*BD_nz_tpp*sizeof(double),cudaHostRegisterMapped);
//   // cudaHostRegister(rho, BD_nx_tpp*BD_ny_tpp*BD_nz_tpp*sizeof(double),cudaHostRegisterMapped);
//   // cudaHostRegister(lambda, BD_nx_tpp*BD_ny_tpp*BD_nz_tpp*sizeof(double),cudaHostRegisterMapped);
// }
