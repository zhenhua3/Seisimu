//////////////////// elastic 3D kernel //////////////////
  __global__ void kernel_v(double *dev_vx, int BD_nx_vx, int BD_nz_vx,
    double *dev_vy, int BD_nx_vy, int BD_nz_vy,
    double *dev_vz, int BD_nx_vz, int BD_nz_vz,
    double *dev_txx, double *dev_tyy, double *dev_tzz,
    int BD_nx_tpp, int BD_nz_tpp,
    double *dev_txy, int BD_nx_txy, int BD_nz_txy,
    double *dev_tyz, int BD_nx_tyz, int BD_nz_tyz,
    double *dev_txz, int BD_nx_txz, int BD_nz_txz,
    double *dev_rho, double *dev_fdc,
    double dx, double dy, double dz, double dt,
    int chunk_full, int chunk_half)
  {
    int tX = threadIdx.x + blockIdx.x*blockDim.x;
    int tY = threadIdx.y + blockIdx.y*blockDim.y;
    int tZ = threadIdx.z + blockIdx.z*blockDim.z;

    int tid_vx, tid_vy, tid_vz;
    int tid_tpp_0, tid_tpp_1, tid_tpp_2, tid_tpp_3;
    int tid_txy_0, tid_txy_1, tid_txy_2, tid_txy_3;
    int tid_tyz_0, tid_tyz_1, tid_tyz_2, tid_tyz_3;
    int tid_txz_0, tid_txz_1, tid_txz_2, tid_txz_3;
    int tid_rho_0, tid_rho_1;
    double tmp, tmp_rho;
    //***************** vx ********************//
    tid_vx = tZ + tX*BD_nz_vx + tY*BD_nx_vx*BD_nz_vx;
    tid_rho_0 = tZ + (tX+0)*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
    tid_rho_1 = tZ + (tX+1)*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;

    //***************** vx by txx ********************//
    if(tX < BD_nx_vx-1 && tX > 0 && tY < chunk_full && tZ < BD_nz_vx)
    {
      tmp_rho = (0.5*(*(dev_rho + tid_rho_0) + *(dev_rho + tid_rho_1)))/dt;

      tid_tpp_0 = tZ + (tX-1)*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
      tid_tpp_1 = tZ + (tX+0)*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
      tid_tpp_2 = tZ + (tX+1)*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
      tid_tpp_3 = tZ + (tX+2)*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;

      *(dev_vx+tid_vx) = *(dev_vx+tid_vx)
      + ((*(dev_txx+tid_tpp_0)* *(dev_fdc)
      + *(dev_txx+tid_tpp_1)* *(dev_fdc+1)
      + *(dev_txx+tid_tpp_2)* *(dev_fdc+2)
      + *(dev_txx+tid_tpp_3)* *(dev_fdc+3)) / dx)/ tmp_rho;
    }

    //***************** vx by txy ********************//
    if(tX < BD_nx_vx && tX > 0 && tY < chunk_full && tZ < BD_nz_vx)
    {
      tmp_rho = (0.5*(*(dev_rho + tid_rho_0) + *(dev_rho + tid_rho_1)))/dt;

      tid_txy_0 = tZ + tX*BD_nz_txy + (tY+0)*BD_nx_txy*BD_nz_txy;
      tid_txy_1 = tZ + tX*BD_nz_txy + (tY+1)*BD_nx_txy*BD_nz_txy;
      tid_txy_2 = tZ + tX*BD_nz_txy + (tY+2)*BD_nx_txy*BD_nz_txy;
      tid_txy_3 = tZ + tX*BD_nz_txy + (tY+3)*BD_nx_txy*BD_nz_txy;

      *(dev_vx+tid_vx) = *(dev_vx+tid_vx)
      + ((*(dev_txy+tid_txy_0)* *(dev_fdc)
      + *(dev_txy+tid_txy_1)* *(dev_fdc+1)
      + *(dev_txy+tid_txy_2)* *(dev_fdc+2)
      + *(dev_txy+tid_txy_3)* *(dev_fdc+3)) / dy) / tmp_rho;
    }

    // ***************** vx by txz ********************//
    if(tX < BD_nx_vx && tY < chunk_full && tZ < BD_nz_vx-2 && tZ > 1)
    {
      tmp_rho = (0.5*(*(dev_rho + tid_rho_0) + *(dev_rho + tid_rho_1)))/dt;

      tid_txz_0 = (tZ-2) + tX*BD_nz_txz + tY*BD_nx_txz*BD_nz_txz;
      tid_txz_1 = (tZ-1) + tX*BD_nz_txz + tY*BD_nx_txz*BD_nz_txz;
      tid_txz_2 = (tZ+0) + tX*BD_nz_txz + tY*BD_nx_txz*BD_nz_txz;
      tid_txz_3 = (tZ+1) + tX*BD_nz_txz + tY*BD_nx_txz*BD_nz_txz;

      *(dev_vx+tid_vx) = *(dev_vx+tid_vx)
      + ((*(dev_txz+tid_txz_0)* *(dev_fdc)
      + *(dev_txz+tid_txz_1)* *(dev_fdc+1)
      + *(dev_txz+tid_txz_2)* *(dev_fdc+2)
      + *(dev_txz+tid_txz_3)* *(dev_fdc+3)) / dz) / tmp_rho;
    }

    //***************** vy ********************//
    tid_vy = tZ + tX*BD_nz_vy + tY*BD_nx_vy*BD_nz_vy;
    tid_rho_0 = tZ + tX*BD_nz_tpp + (tY+0)*BD_nx_tpp*BD_nz_tpp;
    tid_rho_1 = tZ + tX*BD_nz_tpp + (tY+1)*BD_nx_tpp*BD_nz_tpp;

    //***************** vy by txy********************//
    if(tX < BD_nx_vy-2 && tX > 1 && tY < chunk_half && tZ < BD_nz_vy)
    {
      tmp_rho = (0.5*(*(dev_rho + tid_rho_0) + *(dev_rho + tid_rho_1)))/dt;

      tid_txy_0 = tZ + (tX-2)*BD_nz_txy + (tY+1)*BD_nx_txy*BD_nz_txy;
      tid_txy_1 = tZ + (tX-1)*BD_nz_txy + (tY+1)*BD_nx_txy*BD_nz_txy;
      tid_txy_2 = tZ + (tX+0)*BD_nz_txy + (tY+1)*BD_nx_txy*BD_nz_txy;
      tid_txy_3 = tZ + (tX+1)*BD_nz_txy + (tY+1)*BD_nx_txy*BD_nz_txy;

      *(dev_vy+tid_vy) = *(dev_vy+tid_vy)
      + ((*(dev_txy+tid_txy_0)* *(dev_fdc)
      + *(dev_txy+tid_txy_1)* *(dev_fdc+1)
      + *(dev_txy+tid_txy_2)* *(dev_fdc+2)
      + *(dev_txy+tid_txy_3)* *(dev_fdc+3)) / dx) / tmp_rho;
    }

    //***************** vy by tyy********************//
    if(tX < BD_nx_vy && tY < chunk_half && tZ < BD_nz_vy)
    {
      tmp_rho = (0.5*(*(dev_rho + tid_rho_0) + *(dev_rho + tid_rho_1)))/dt;

      tid_tpp_0 = tZ + tX*BD_nz_tpp + (tY+0)*BD_nx_tpp*BD_nz_tpp;
      tid_tpp_1 = tZ + tX*BD_nz_tpp + (tY+1)*BD_nx_tpp*BD_nz_tpp;
      tid_tpp_2 = tZ + tX*BD_nz_tpp + (tY+2)*BD_nx_tpp*BD_nz_tpp;
      tid_tpp_3 = tZ + tX*BD_nz_tpp + (tY+3)*BD_nx_tpp*BD_nz_tpp;

      *(dev_vy+tid_vy) = *(dev_vy+tid_vy)
      + ((*(dev_tyy+tid_tpp_0)* *(dev_fdc)
      + *(dev_tyy+tid_tpp_1)* *(dev_fdc+1)
      + *(dev_tyy+tid_tpp_2)* *(dev_fdc+2)
      + *(dev_tyy+tid_tpp_3)* *(dev_fdc+3)) / dy) / tmp_rho;
    }

    //***************** vy by tyz********************//
    if(tX < BD_nx_vy && tY < chunk_half && tZ < BD_nz_vy - 2 && tZ > 1)
    {
      tmp_rho = (0.5*(*(dev_rho + tid_rho_0) + *(dev_rho + tid_rho_1)))/dt;

      tid_tyz_0 = (tZ-2) + tX*BD_nz_tyz + (tY+1)*BD_nx_tyz*BD_nz_tyz;
      tid_tyz_1 = (tZ-1) + tX*BD_nz_tyz + (tY+1)*BD_nx_tyz*BD_nz_tyz;
      tid_tyz_2 = (tZ+0) + tX*BD_nz_tyz + (tY+1)*BD_nx_tyz*BD_nz_tyz;
      tid_tyz_3 = (tZ+1) + tX*BD_nz_tyz + (tY+1)*BD_nx_tyz*BD_nz_tyz;

      tmp = ((*(dev_tyz+tid_tyz_0)* *(dev_fdc)
      + *(dev_tyz+tid_tyz_1)* *(dev_fdc+1)
      + *(dev_tyz+tid_tyz_2)* *(dev_fdc+2)
      + *(dev_tyz+tid_tyz_3)* *(dev_fdc+3)) / dz) / tmp_rho;
      *(dev_vy+tid_vy) = *(dev_vy+tid_vy) + tmp;
    }

    //***************** vz ********************//
    tid_vz = tZ + tX*BD_nz_vz + tY*BD_nx_vz*BD_nz_vz;
    tid_rho_0 = (tZ+0) + tX*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
    tid_rho_1 = (tZ+1) + tX*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;


    //***************** vz by txz ********************//
    if(tX < BD_nx_vz-2 && tX > 1 && tY < chunk_full && tZ < BD_nz_vz)
    {
      tmp_rho = (0.5*(*(dev_rho + tid_rho_0) + *(dev_rho + tid_rho_1)))/dt;

      tid_txz_0 = tZ + (tX-2)*BD_nz_txz + tY*BD_nx_txz*BD_nz_txz;
      tid_txz_1 = tZ + (tX-1)*BD_nz_txz + tY*BD_nx_txz*BD_nz_txz;
      tid_txz_2 = tZ + (tX+0)*BD_nz_txz + tY*BD_nx_txz*BD_nz_txz;
      tid_txz_3 = tZ + (tX+1)*BD_nz_txz + tY*BD_nx_txz*BD_nz_txz;

      tmp = ((*(dev_txz+tid_txz_0)* *(dev_fdc)
      + *(dev_txz+tid_txz_1)* *(dev_fdc+1)
      + *(dev_txz+tid_txz_2)* *(dev_fdc+2)
      + *(dev_txz+tid_txz_3)* *(dev_fdc+3)) / dx) / tmp_rho;
      *(dev_vz+tid_vz) = *(dev_vz+tid_vz) + tmp;
    }

    //***************** vz by tyz ********************//
    if(tX < BD_nx_vz && tY < chunk_full && tZ < BD_nz_vz)
    {
      tmp_rho = (0.5*(*(dev_rho + tid_rho_0) + *(dev_rho + tid_rho_1)))/dt;

      tid_tyz_0 = tZ + tX*BD_nz_tyz + (tY+0)*BD_nx_tyz*BD_nz_tyz;
      tid_tyz_1 = tZ + tX*BD_nz_tyz + (tY+1)*BD_nx_tyz*BD_nz_tyz;
      tid_tyz_2 = tZ + tX*BD_nz_tyz + (tY+2)*BD_nx_tyz*BD_nz_tyz;
      tid_tyz_3 = tZ + tX*BD_nz_tyz + (tY+3)*BD_nx_tyz*BD_nz_tyz;

      tmp = ((*(dev_tyz+tid_tyz_0)* *(dev_fdc)
      + *(dev_tyz+tid_tyz_1)* *(dev_fdc+1)
      + *(dev_tyz+tid_tyz_2)* *(dev_fdc+2)
      + *(dev_tyz+tid_tyz_3)* *(dev_fdc+3)) / dy) / tmp_rho;
      *(dev_vz+tid_vz) = *(dev_vz+tid_vz) + tmp;
    }

    //***************** vz by tzz ********************//
    if(tX < BD_nx_vz && tY < chunk_full && tZ < BD_nz_vz-1 && tZ > 0)
    {
      tmp_rho = (0.5*(*(dev_rho + tid_rho_0) + *(dev_rho + tid_rho_1)))/dt;

      tid_tpp_0 = (tZ-1) + tX*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
      tid_tpp_1 = (tZ+0) + tX*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
      tid_tpp_2 = (tZ+1) + tX*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
      tid_tpp_3 = (tZ+2) + tX*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;

      tmp = ((*(dev_tzz+tid_tpp_0)* *(dev_fdc)
      + *(dev_tzz+tid_tpp_1)* *(dev_fdc+1)
      + *(dev_tzz+tid_tpp_2)* *(dev_fdc+2)
      + *(dev_tzz+tid_tpp_3)* *(dev_fdc+3)) / dz) / tmp_rho;
      *(dev_vz+tid_vz) = *(dev_vz+tid_vz) +  tmp;
    }
  }
  __global__ void kernel_tau(double *dev_vx, int BD_nx_vx, int BD_nz_vx,
    double *dev_vy, int BD_nx_vy, int BD_nz_vy,
    double *dev_vz, int BD_nx_vz, int BD_nz_vz,
    double *dev_txx, double *dev_tyy, double *dev_tzz,
    int BD_nx_tpp, int BD_nz_tpp,
    double *dev_txy, int BD_nx_txy, int BD_nz_txy,
    double *dev_tyz, int BD_nx_tyz, int BD_nz_tyz,
    double *dev_txz, int BD_nx_txz, int BD_nz_txz,
    double *dev_lambda, double *dev_mu, double *dev_fdc,
    double dx, double dy, double dz, double dt,
    int chunk_full, int chunk_half)
    {
      int tX = threadIdx.x + blockIdx.x*blockDim.x;
      int tY = threadIdx.y + blockIdx.y*blockDim.y;
      int tZ = threadIdx.z + blockIdx.z*blockDim.z;

      int tid_tpp, tid_txy, tid_txz, tid_tyz, tid_lambda_tpp, tid_mu_tpp;
      int tid_vx_0, tid_vx_1, tid_vx_2, tid_vx_3;
      int tid_vy_0, tid_vy_1, tid_vy_2, tid_vy_3;
      int tid_vz_0, tid_vz_1, tid_vz_2, tid_vz_3;
      int tid_mu_0, tid_mu_1, tid_mu_2, tid_mu_3;
      double tmp_vx, tmp_vy, tmp_vz, tmp_mu;

      //**************** tau_pp ********************//
      tid_tpp = tZ + tX*BD_nz_tpp + tY*BD_nz_tpp*BD_nx_tpp;
      tid_lambda_tpp = tZ + tX*BD_nz_tpp + tY*BD_nz_tpp*BD_nx_tpp;
      tid_mu_tpp = tZ + tX*BD_nz_tpp + (tY+1)*BD_nz_tpp*BD_nx_tpp;
      //**************** tau_pp by vx ********************//
      if(tX > 1 && tX < BD_nx_tpp-2 && tY < chunk_full && tZ < BD_nz_tpp)
      {
        tid_vx_0 = tZ + (tX-2)*BD_nz_vx + (tY+2)*BD_nz_vx*BD_nx_vx;
        tid_vx_1 = tZ + (tX-1)*BD_nz_vx + (tY+2)*BD_nz_vx*BD_nx_vx;
        tid_vx_2 = tZ + (tX+0)*BD_nz_vx + (tY+2)*BD_nz_vx*BD_nx_vx;
        tid_vx_3 = tZ + (tX+1)*BD_nz_vx + (tY+2)*BD_nz_vx*BD_nx_vx;

        tmp_vx = ((*(dev_vx+tid_vx_0)* *(dev_fdc)
        + *(dev_vx+tid_vx_1)* *(dev_fdc+1)
        + *(dev_vx+tid_vx_2)* *(dev_fdc+2)
        + *(dev_vx+tid_vx_3)* *(dev_fdc+3)) / dx)*dt;

        *(dev_txx+tid_tpp) = *(dev_txx+tid_tpp)
        + (*(dev_lambda+tid_lambda_tpp) + 2* *(dev_mu+tid_mu_tpp)) * tmp_vx;

        *(dev_tyy+tid_tpp) = *(dev_tyy+tid_tpp)
        + *(dev_lambda+tid_lambda_tpp) * tmp_vx;

        *(dev_tzz+tid_tpp) = *(dev_tzz+tid_tpp)
        + *(dev_lambda+tid_lambda_tpp) * tmp_vx;
      }
      //**************** tau_pp by vy ********************//
      if(tX < BD_nx_tpp && tY < chunk_full && tZ < BD_nz_tpp)
      {
        tid_vy_0 = tZ + tX*BD_nz_vy + (tY+0)*BD_nz_vy*BD_nx_vy;
        tid_vy_1 = tZ + tX*BD_nz_vy + (tY+1)*BD_nz_vy*BD_nx_vy;
        tid_vy_2 = tZ + tX*BD_nz_vy + (tY+2)*BD_nz_vy*BD_nx_vy;
        tid_vy_3 = tZ + tX*BD_nz_vy + (tY+3)*BD_nz_vy*BD_nx_vy;

        tmp_vy = ((*(dev_vy+tid_vy_0)* *(dev_fdc)
        + *(dev_vy+tid_vy_1)* *(dev_fdc+1)
        + *(dev_vy+tid_vy_2)* *(dev_fdc+2)
        + *(dev_vy+tid_vy_3)* *(dev_fdc+3)) / dy)*dt;

        *(dev_txx+tid_tpp) = *(dev_txx+tid_tpp)
        + *(dev_lambda+tid_lambda_tpp) * tmp_vy;

        *(dev_tyy+tid_tpp) = *(dev_tyy+tid_tpp)
        + (*(dev_lambda+tid_lambda_tpp) + 2* *(dev_mu+tid_mu_tpp)) * tmp_vy;

        *(dev_tzz+tid_tpp) = *(dev_tzz+tid_tpp)
        + *(dev_lambda+tid_lambda_tpp) * tmp_vy;
      }
      //**************** tau_pp by vz ********************//
      if(tX < BD_nx_tpp && tY < chunk_full && tZ > 1 && tZ < BD_nz_tpp-2)
      {
        tid_vz_0 = (tZ-2) + tX*BD_nz_vz + (tY+2)*BD_nz_vz*BD_nx_vz;
        tid_vz_1 = (tZ-1) + tX*BD_nz_vz + (tY+2)*BD_nz_vz*BD_nx_vz;
        tid_vz_2 = (tZ+0) + tX*BD_nz_vz + (tY+2)*BD_nz_vz*BD_nx_vz;
        tid_vz_3 = (tZ+1) + tX*BD_nz_vz + (tY+2)*BD_nz_vz*BD_nx_vz;

        tmp_vz = ((*(dev_vz+tid_vz_0)* *(dev_fdc)
        + *(dev_vz+tid_vz_1)* *(dev_fdc+1)
        + *(dev_vz+tid_vz_2)* *(dev_fdc+2)
        + *(dev_vz+tid_vz_3)* *(dev_fdc+3)) / dz)*dt;

        *(dev_txx+tid_tpp) = *(dev_txx+tid_tpp)
        + *(dev_lambda+tid_lambda_tpp) * tmp_vz;

        *(dev_tyy+tid_tpp) = *(dev_tyy+tid_tpp)
        + *(dev_lambda+tid_lambda_tpp) * tmp_vz;

        *(dev_tzz+tid_tpp) = *(dev_tzz+tid_tpp)
        + (*(dev_lambda+tid_lambda_tpp) + 2* *(dev_mu+tid_mu_tpp)) * tmp_vz;
      }
      //**************** tau_xy ********************//
      tid_txy = tZ + tX*BD_nz_txy + tY*BD_nz_txy*BD_nx_txy;
      tid_mu_0 = tZ + tX*BD_nz_tpp + tY*BD_nz_tpp*BD_nx_tpp;
      tid_mu_1 = tZ + tX*BD_nz_tpp + (tY+1)*BD_nz_tpp*BD_nx_tpp;
      tid_mu_2 = tZ + (tX+1)*BD_nz_tpp + tY*BD_nz_tpp*BD_nx_tpp;
      tid_mu_3 = tZ + (tX+1)*BD_nz_tpp + (tY+1)*BD_nz_tpp*BD_nx_tpp;

      //**************** tau_xy by vy ********************//
      if(tX > 0 && tX < BD_nx_txy-1 && tY < chunk_half && tZ < BD_nz_txy)
      {
        tmp_mu = (*(dev_mu+tid_mu_0) + *(dev_mu+tid_mu_1)
        + *(dev_mu+tid_mu_2) + *(dev_mu+tid_mu_3))/4;

        tid_vy_0 = tZ + (tX-1)*BD_nz_vy + (tY+1)*BD_nz_vy*BD_nx_vy;
        tid_vy_1 = tZ + (tX+0)*BD_nz_vy + (tY+1)*BD_nz_vy*BD_nx_vy;
        tid_vy_2 = tZ + (tX+1)*BD_nz_vy + (tY+1)*BD_nz_vy*BD_nx_vy;
        tid_vy_3 = tZ + (tX+2)*BD_nz_vy + (tY+1)*BD_nz_vy*BD_nx_vy;

        tmp_vy = ((*(dev_vy+tid_vy_0)* *(dev_fdc)
        + *(dev_vy+tid_vy_1)* *(dev_fdc+1)
        + *(dev_vy+tid_vy_2)* *(dev_fdc+2)
        + *(dev_vy+tid_vy_3)* *(dev_fdc+3)) / dx)*dt;
        *(dev_txy+tid_txy) = *(dev_txy+tid_txy) + tmp_vy * tmp_mu;
      }

      //**************** tau_xy by vx ********************//
      if(tX < BD_nx_txy && tY < chunk_full && tZ < BD_nz_txy)
      {
        tmp_mu = (*(dev_mu+tid_mu_0) + *(dev_mu+tid_mu_1)
        + *(dev_mu+tid_mu_2) + *(dev_mu+tid_mu_3))/4;

        tid_vx_0 = tZ + tX*BD_nz_vx + (tY+0)*BD_nz_vx*BD_nx_vx;
        tid_vx_1 = tZ + tX*BD_nz_vx + (tY+1)*BD_nz_vx*BD_nx_vx;
        tid_vx_2 = tZ + tX*BD_nz_vx + (tY+2)*BD_nz_vx*BD_nx_vx;
        tid_vx_3 = tZ + tX*BD_nz_vx + (tY+3)*BD_nz_vx*BD_nx_vx;

        tmp_vx = ((*(dev_vx+tid_vx_0)* *(dev_fdc)
        + *(dev_vx+tid_vx_1)* *(dev_fdc+1)
        + *(dev_vx+tid_vx_2)* *(dev_fdc+2)
        + *(dev_vx+tid_vx_3)* *(dev_fdc+3)) / dy)*dt;

        *(dev_txy+tid_txy) = *(dev_txy+tid_txy) + tmp_mu * tmp_vx;
      }

      //**************** tau_xz ********************//
      tid_txz = tZ + tX*BD_nz_txz + tY*BD_nz_txz*BD_nx_txz;
      tid_mu_0 = tZ + tX*BD_nz_tpp + tY*BD_nz_tpp*BD_nx_tpp;
      tid_mu_1 = (tZ+1) + tX*BD_nz_tpp + tY*BD_nz_tpp*BD_nx_tpp;
      tid_mu_2 = tZ + (tX+1)*BD_nz_tpp + tY*BD_nz_tpp*BD_nx_tpp;
      tid_mu_3 = (tZ+1) + (tX+1)*BD_nz_tpp + tY*BD_nz_tpp*BD_nx_tpp;

      //**************** tau_xz by vx ********************//
      if(tX < BD_nx_txz && tY < chunk_half && tZ < BD_nz_txz-1 && tZ > 0)
      {
        tmp_mu = (*(dev_mu+tid_mu_0) + *(dev_mu+tid_mu_1)
        + *(dev_mu+tid_mu_2) + *(dev_mu+tid_mu_3))/4;

        tid_vx_0 = (tZ-1) + tX*BD_nz_vx + (tY+2)*BD_nz_vx*BD_nx_vx;
        tid_vx_1 = (tZ+0) + tX*BD_nz_vx + (tY+2)*BD_nz_vx*BD_nx_vx;
        tid_vx_2 = (tZ+1) + tX*BD_nz_vx + (tY+2)*BD_nz_vx*BD_nx_vx;
        tid_vx_3 = (tZ+2) + tX*BD_nz_vx + (tY+2)*BD_nz_vx*BD_nx_vx;

        tmp_vx = ((*(dev_vx+tid_vx_0)* *(dev_fdc)
        + *(dev_vx+tid_vx_1)* *(dev_fdc+1)
        + *(dev_vx+tid_vx_2)* *(dev_fdc+2)
        + *(dev_vx+tid_vx_3)* *(dev_fdc+3)) / dz)*dt;

        *(dev_txz+tid_txz) = *(dev_txz+tid_txz) + tmp_mu * tmp_vx;
      }

      //**************** tau_xz by vz ********************//
      if(tX > 0 && tX < BD_nx_txz-1 && tY < chunk_half && tZ < BD_nz_txz)
      {
        tmp_mu = (*(dev_mu+tid_mu_0) + *(dev_mu+tid_mu_1)
        + *(dev_mu+tid_mu_2) + *(dev_mu+tid_mu_3))/4;

        tid_vz_0 = tZ + (tX-1)*BD_nz_vz + (tY+2)*BD_nz_vz*BD_nx_vz;
        tid_vz_1 = tZ + (tX+0)*BD_nz_vz + (tY+2)*BD_nz_vz*BD_nx_vz;
        tid_vz_2 = tZ + (tX+1)*BD_nz_vz + (tY+2)*BD_nz_vz*BD_nx_vz;
        tid_vz_3 = tZ + (tX+2)*BD_nz_vz + (tY+2)*BD_nz_vz*BD_nx_vz;

        tmp_vz = ((*(dev_vz+tid_vz_0)* *(dev_fdc)
        + *(dev_vz+tid_vz_1)* *(dev_fdc+1)
        + *(dev_vz+tid_vz_2)* *(dev_fdc+2)
        + *(dev_vz+tid_vz_3)* *(dev_fdc+3)) / dx)*dt;

        *(dev_txz+tid_txz) = *(dev_txz+tid_txz) + tmp_mu * tmp_vz;
      }

      //**************** tau_yz ********************//
      tid_tyz = tZ + tX*BD_nz_tyz + tY*BD_nz_tyz*BD_nx_tyz;
      tid_mu_0 = tZ + tX*BD_nz_tpp + tY*BD_nz_tpp*BD_nx_tpp;
      tid_mu_1 = (tZ+1) + tX*BD_nz_tpp + tY*BD_nz_tpp*BD_nx_tpp;
      tid_mu_2 = tZ + tX*BD_nz_tpp + (tY+1)*BD_nz_tpp*BD_nx_tpp;
      tid_mu_3 = (tZ+1) + tX*BD_nz_tpp + (tY+1)*BD_nz_tpp*BD_nx_tpp;

      //**************** tau_yz by vy ********************//
      if(tX < BD_nx_tyz && tY < chunk_full && tZ < BD_nz_tyz-1 && tZ > 0)
      {
        tmp_mu = (*(dev_mu+tid_mu_0) + *(dev_mu+tid_mu_1)
        + *(dev_mu+tid_mu_2) + *(dev_mu+tid_mu_3))/4;

        tid_vy_0 = (tZ-1) + tX*BD_nz_vy + (tY+1)*BD_nz_vy*BD_nx_vy;
        tid_vy_1 = (tZ+0) + tX*BD_nz_vy + (tY+1)*BD_nz_vy*BD_nx_vy;
        tid_vy_2 = (tZ+1) + tX*BD_nz_vy + (tY+1)*BD_nz_vy*BD_nx_vy;
        tid_vy_3 = (tZ+2) + tX*BD_nz_vy + (tY+1)*BD_nz_vy*BD_nx_vy;

        tmp_vy = ((*(dev_vy+tid_vy_0)* *(dev_fdc)
        + *(dev_vy+tid_vy_1)* *(dev_fdc+1)
        + *(dev_vy+tid_vy_2)* *(dev_fdc+2)
        + *(dev_vy+tid_vy_3)* *(dev_fdc+3)) / dz)*dt;

        *(dev_tyz+tid_tyz) = *(dev_tyz+tid_tyz) + tmp_mu * tmp_vy;
      }

      //**************** tau_yz by vz ********************//
      if(tX < BD_nx_tyz && tY < chunk_full && tZ < BD_nz_txz)
      {
        tmp_mu = (*(dev_mu+tid_mu_0) + *(dev_mu+tid_mu_1)
        + *(dev_mu+tid_mu_2) + *(dev_mu+tid_mu_3))/4;

        tid_vz_0 = tZ + tX*BD_nz_vz + (tY+0)*BD_nz_vz*BD_nx_vz;
        tid_vz_1 = tZ + tX*BD_nz_vz + (tY+1)*BD_nz_vz*BD_nx_vz;
        tid_vz_2 = tZ + tX*BD_nz_vz + (tY+2)*BD_nz_vz*BD_nx_vz;
        tid_vz_3 = tZ + tX*BD_nz_vz + (tY+3)*BD_nz_vz*BD_nx_vz;

        tmp_vz = ((*(dev_vz+tid_vz_0)* *(dev_fdc)
        + *(dev_vz+tid_vz_1)* *(dev_fdc+1)
        + *(dev_vz+tid_vz_2)* *(dev_fdc+2)
        + *(dev_vz+tid_vz_3)* *(dev_fdc+3)) / dy)*dt;

        *(dev_tyz+tid_tyz) = *(dev_tyz+tid_tyz) + tmp_mu * tmp_vz;
      }
    }
