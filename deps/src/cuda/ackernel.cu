__constant__ double dev_fdc[4];

////////////// acoustic 3D kernel //////////////////
__global__ void ackernel_v(double *dev_vx, int BD_nx_vx, int BD_ny_vx, int BD_nz_vx,
  double *dev_vy, int BD_nx_vy, int BD_ny_vy, int BD_nz_vy,
  double *dev_vz, int BD_nx_vz, int BD_ny_vz, int BD_nz_vz,
  double *dev_tpp, int BD_nx_tpp, int BD_ny_tpp, int BD_nz_tpp,
  double *dev_rho, double dx, double dy, double dz, double dt)
{
  int tX = threadIdx.x + blockIdx.x*blockDim.x;
  int tY = threadIdx.y + blockIdx.y*blockDim.y;
  int tZ = threadIdx.z + blockIdx.z*blockDim.z;

  int tid_vx, tid_vy, tid_vz;
  int tid_tpp_0, tid_tpp_1, tid_tpp_2, tid_tpp_3;
  int tid_rho_0, tid_rho_1;
  double tmp_rho;
  //***************** vx ********************//
  tid_vx = tZ + tX*BD_nz_vx + tY*BD_nx_vx*BD_nz_vx;
  tid_rho_0 = tZ + (tX+0)*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
  tid_rho_1 = tZ + (tX+1)*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;

  //***************** vx by tpp ********************//
  if(tX < BD_nx_vx-1 && tX > 0 && tY < BD_ny_vx && tZ < BD_nz_vx)
  {
    tmp_rho = (0.5*(*(dev_rho + tid_rho_0) + *(dev_rho + tid_rho_1)))/dt;

    tid_tpp_0 = tZ + (tX-1)*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
    tid_tpp_1 = tZ + (tX+0)*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
    tid_tpp_2 = tZ + (tX+1)*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
    tid_tpp_3 = tZ + (tX+2)*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
    //
    *(dev_vx+tid_vx) = *(dev_vx+tid_vx)
    + ((*(dev_tpp+tid_tpp_0)* *(dev_fdc)
    + *(dev_tpp+tid_tpp_1)* *(dev_fdc+1)
    + *(dev_tpp+tid_tpp_2)* *(dev_fdc+2)
    + *(dev_tpp+tid_tpp_3)* *(dev_fdc+3)) / dx)/ tmp_rho;
  }

  //***************** vy ********************//
  tid_vy = tZ + tX*BD_nz_vy + tY*BD_nx_vy*BD_nz_vy;
  tid_rho_0 = tZ + tX*BD_nz_tpp + (tY+0)*BD_nx_tpp*BD_nz_tpp;
  tid_rho_1 = tZ + tX*BD_nz_tpp + (tY+1)*BD_nx_tpp*BD_nz_tpp;

  //***************** vy by tpp********************//
  if(tX < BD_nx_vy && tY > 0 && tY < BD_ny_vy-1 && tZ < BD_nz_vy)
  {
    tmp_rho = (0.5*(*(dev_rho + tid_rho_0) + *(dev_rho + tid_rho_1)))/dt;

    tid_tpp_0 = tZ + tX*BD_nz_tpp + (tY-1)*BD_nx_tpp*BD_nz_tpp;
    tid_tpp_1 = tZ + tX*BD_nz_tpp + (tY+0)*BD_nx_tpp*BD_nz_tpp;
    tid_tpp_2 = tZ + tX*BD_nz_tpp + (tY+1)*BD_nx_tpp*BD_nz_tpp;
    tid_tpp_3 = tZ + tX*BD_nz_tpp + (tY+2)*BD_nx_tpp*BD_nz_tpp;

    *(dev_vy+tid_vy) = *(dev_vy+tid_vy)
    + ((*(dev_tpp+tid_tpp_0)* *(dev_fdc)
    + *(dev_tpp+tid_tpp_1)* *(dev_fdc+1)
    + *(dev_tpp+tid_tpp_2)* *(dev_fdc+2)
    + *(dev_tpp+tid_tpp_3)* *(dev_fdc+3)) / dy) / tmp_rho;
  }

  //***************** vz ********************//
  tid_vz = tZ + tX*BD_nz_vz + tY*BD_nx_vz*BD_nz_vz;
  tid_rho_0 = (tZ+0) + tX*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
  tid_rho_1 = (tZ+1) + tX*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;

  //***************** vz by tpp ********************//
  if(tX < BD_nx_vz && tY < BD_ny_vz && tZ < BD_nz_vz-1 && tZ > 0)
  {
    tmp_rho = (0.5*(*(dev_rho + tid_rho_0) + *(dev_rho + tid_rho_1)))/dt;

    tid_tpp_0 = (tZ-1) + tX*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
    tid_tpp_1 = (tZ+0) + tX*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
    tid_tpp_2 = (tZ+1) + tX*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
    tid_tpp_3 = (tZ+2) + tX*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;

    *(dev_vz+tid_vz) = *(dev_vz+tid_vz) + ((*(dev_tpp+tid_tpp_0)* *(dev_fdc)
    + *(dev_tpp+tid_tpp_1)* *(dev_fdc+1)
    + *(dev_tpp+tid_tpp_2)* *(dev_fdc+2)
    + *(dev_tpp+tid_tpp_3)* *(dev_fdc+3)) / dz) / tmp_rho;

  }
}

__global__ void ackernel_tau(double *dev_vx, int BD_nx_vx, int BD_ny_vx, int BD_nz_vx,
  double *dev_vy, int BD_nx_vy, int BD_ny_vy, int BD_nz_vy,
  double *dev_vz, int BD_nx_vz, int BD_ny_vz, int BD_nz_vz,
  double *dev_tpp, int BD_nx_tpp, int BD_ny_tpp, int BD_nz_tpp,
  double *dev_lambda,
  double dx, double dy, double dz, double dt)
  {
    int tX = threadIdx.x + blockIdx.x*blockDim.x;
    int tY = threadIdx.y + blockIdx.y*blockDim.y;
    int tZ = threadIdx.z + blockIdx.z*blockDim.z;

    int tid_tpp, tid_lambda_tpp;
    int tid_vx_0, tid_vx_1, tid_vx_2, tid_vx_3;
    int tid_vy_0, tid_vy_1, tid_vy_2, tid_vy_3;
    int tid_vz_0, tid_vz_1, tid_vz_2, tid_vz_3;
    double tmp_vx, tmp_vy, tmp_vz;

    //**************** tau_pp ********************//
    tid_tpp = tZ + tX*BD_nz_tpp + tY*BD_nz_tpp*BD_nx_tpp;
    tid_lambda_tpp = tZ + tX*BD_nz_tpp + tY*BD_nz_tpp*BD_nx_tpp;

    //**************** tau_pp by vx ********************//
    if(tX > 1 && tX < BD_nx_tpp-2 && tY < BD_ny_tpp && tZ < BD_nz_tpp)
    {
      tid_vx_0 = tZ + (tX-2)*BD_nz_vx + tY*BD_nz_vx*BD_nx_vx;
      tid_vx_1 = tZ + (tX-1)*BD_nz_vx + tY*BD_nz_vx*BD_nx_vx;
      tid_vx_2 = tZ + (tX+0)*BD_nz_vx + tY*BD_nz_vx*BD_nx_vx;
      tid_vx_3 = tZ + (tX+1)*BD_nz_vx + tY*BD_nz_vx*BD_nx_vx;

      tmp_vx = ((*(dev_vx+tid_vx_0)* *(dev_fdc)
      + *(dev_vx+tid_vx_1)* *(dev_fdc+1)
      + *(dev_vx+tid_vx_2)* *(dev_fdc+2)
      + *(dev_vx+tid_vx_3)* *(dev_fdc+3)) / dx)*dt;

      *(dev_tpp+tid_tpp) = *(dev_tpp+tid_tpp)
      + *(dev_lambda+tid_lambda_tpp) * tmp_vx;
    }
    //**************** tau_pp by vy ********************//
    if(tX < BD_nx_tpp && tY > 1 && tY < BD_ny_tpp-2  && tZ < BD_nz_tpp)
    {
      tid_vy_0 = tZ + tX*BD_nz_vy + (tY-2)*BD_nz_vy*BD_nx_vy;
      tid_vy_1 = tZ + tX*BD_nz_vy + (tY-1)*BD_nz_vy*BD_nx_vy;
      tid_vy_2 = tZ + tX*BD_nz_vy + (tY+0)*BD_nz_vy*BD_nx_vy;
      tid_vy_3 = tZ + tX*BD_nz_vy + (tY+1)*BD_nz_vy*BD_nx_vy;

      tmp_vy = ((*(dev_vy+tid_vy_0)* *(dev_fdc)
      + *(dev_vy+tid_vy_1)* *(dev_fdc+1)
      + *(dev_vy+tid_vy_2)* *(dev_fdc+2)
      + *(dev_vy+tid_vy_3)* *(dev_fdc+3)) / dy)*dt;

      *(dev_tpp+tid_tpp) = *(dev_tpp+tid_tpp)
      + *(dev_lambda+tid_lambda_tpp) * tmp_vy;
    }
    //**************** tau_pp by vz ********************//
    if(tX < BD_nx_tpp && tY < BD_ny_tpp && tZ > 1 && tZ < BD_nz_tpp-2)
    {
      tid_vz_0 = (tZ-2) + tX*BD_nz_vz + tY*BD_nz_vz*BD_nx_vz;
      tid_vz_1 = (tZ-1) + tX*BD_nz_vz + tY*BD_nz_vz*BD_nx_vz;
      tid_vz_2 = (tZ+0) + tX*BD_nz_vz + tY*BD_nz_vz*BD_nx_vz;
      tid_vz_3 = (tZ+1) + tX*BD_nz_vz + tY*BD_nz_vz*BD_nx_vz;

      tmp_vz = ((*(dev_vz+tid_vz_0)* *(dev_fdc)
      + *(dev_vz+tid_vz_1)* *(dev_fdc+1)
      + *(dev_vz+tid_vz_2)* *(dev_fdc+2)
      + *(dev_vz+tid_vz_3)* *(dev_fdc+3)) / dz)*dt;

      *(dev_tpp+tid_tpp) = *(dev_tpp+tid_tpp)
      + *(dev_lambda+tid_lambda_tpp) * tmp_vz;
    }
  }
