#include<stdio.h>
#include"3d_lib.c"

// vx
void ac3d_openmp(double *vx, int BD_nx_vx, int BD_ny_vx, int BD_nz_vx,
  double *pvxbtpp,
  double *vy, int BD_nx_vy, int BD_ny_vy, int BD_nz_vy,
  double *pvybtpp,
  double *vz, int BD_nx_vz, int BD_ny_vz, int BD_nz_vz,
  double *pvzbtpp,
  double *tpp, double *ptppbvx, double *ptppbvy, double *ptppbvz,
  int BD_nx_tpp, int BD_ny_tpp, int BD_nz_tpp,
  double *rho, double *lambda, double *fdc,
  double dt, double dx, double dy, double dz, int ext,
  double *bhalf, double *ahalf, double *bfull, double *afull)
{
  int i,j,k;
//********************* V_X *********************//
  #pragma omp parallel
  {
    // vxbtxx
    #pragma omp for collapse(2) private(j) nowait
    for(k=0; k<BD_ny_vx; k++)
    {
      for(i=0; i<BD_nz_vx; i++)
      {
        for(j=1; j<ext; j++)
        {
          bound_x_3d(vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
            tpp, BD_nz_tpp, BD_nx_tpp, 1,
          2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
          2.0/(*(rho+i+(BD_nx_vx-1-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(BD_nx_vx-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
          dx, dt, pvxbtpp, bhalf, ahalf, ext, fdc);
        }
        for(j=ext; j<BD_nx_vx-ext; j++)
        {
          body_x_3d(vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
            tpp, BD_nz_tpp, BD_nx_tpp, 1,
          2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
          dx, dt, fdc);
        }
      }
    }

    //********************* V_Y *********************//
    // vybtyy
    #pragma omp for collapse(2) private(k) nowait
    for(j=0; j<BD_nx_vy; j++)
    {
      for(i=0; i<BD_nz_vy; i++)
      {
        for(k=1; k<ext; k++)
        {
          bound_y_3d(vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
            tpp, BD_nz_tpp, BD_nx_tpp, 1,
            2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
            2.0/(*(rho+i+j*BD_nz_tpp+(BD_ny_vy-1-k)*BD_nz_tpp*BD_nx_tpp)+*(rho+i+j*BD_nz_tpp+(BD_ny_vy-k)*BD_nz_tpp*BD_nx_tpp)),
            dy,dt, pvybtpp, bhalf, ahalf, ext, fdc);
        }
        for(k=ext; k<BD_ny_vy-ext; k++)
        {
          body_y_3d(vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
            tpp, BD_nz_tpp, BD_nx_tpp, 1,
            2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
            dy,dt, fdc);
        }
      }
    }

    //********************* V_Z *********************//
    // vzbtzz
    #pragma omp for collapse(2) private(i) nowait
    for(j=0; j<BD_nx_vz; j++)
    {
      for(k=0; k<BD_ny_vz; k++)
      {
        for(i=1;i<ext;i++)
        {
          unlimited_bound_z_3d(vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
            tpp, BD_nz_tpp, BD_nx_tpp, 1,
            2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            2.0/(*(rho+(BD_nz_vz-1-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(BD_nz_vz-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            dz,dt, pvzbtpp, bhalf, ahalf, ext, fdc);
        }
        for(i=ext; i<BD_nz_vz-ext; i++)
        {
          body_z_3d(vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
            tpp, BD_nz_tpp, BD_nx_tpp, 1,
            2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            dz,dt, fdc);
        }
      }
    }

    #pragma omp barrier

    //************ T_X ************//
    #pragma omp for collapse(2) private(j)
    for(k=0; k<BD_ny_tpp; k++)
    {
      for(i=0; i<BD_nz_tpp; i++)
      {
        for(j=2; j<ext; j++)
        {
          ac_bound_tpp_x_3d(tpp, BD_nz_tpp, BD_nx_tpp, BD_ny_tpp,
          i, j, k, vx, BD_nz_vx, BD_nx_vx, 2, lambda, dx, dt,
          ptppbvx, bfull, afull, ext, fdc);
        }
        for(j=ext; j<BD_nx_tpp-ext; j++)
        {
          ac_body_tpp_x_3d(tpp, BD_nz_tpp, BD_nx_tpp, BD_ny_tpp,
          i, j, k, vx, BD_nz_vx, BD_nx_vx, 2, lambda, dx, dt, fdc);
        }
      }
    }

    //************ T_Y ************//
    #pragma omp for collapse(2) private(k)
    for(i=0; i<BD_nz_tpp; i++)
    {
      for(j=0; j<BD_nx_tpp; j++)
      {
        for(k=2; k<ext; k++)
        {
          ac_bound_tpp_y_3d(tpp, BD_nz_tpp, BD_nx_tpp, BD_ny_tpp,
          i, j, k, vy, BD_nz_vy, BD_nx_vy, 2, lambda, dy, dt,
          ptppbvy, bfull, afull, ext, fdc);
        }
        for(k=ext; k<BD_ny_tpp-ext; k++)
        {
          ac_body_tpp_y_3d(tpp, BD_nz_tpp, BD_nx_tpp, BD_ny_tpp,
          i, j, k, vy, BD_nz_vy, BD_nx_vy, 2, lambda, dy, dt, fdc);
        }
      }
    }

    //************ T_Z ************//
    #pragma omp for collapse(2) private(i)
    for(j=0; j<BD_nx_tpp; j++)
    {
      for(k=0; k<BD_ny_tpp; k++)
      {
        for(i=2; i<ext; i++)
        {
          ac_unlimited_bound_tpp_z_3d(tpp, BD_nz_tpp, BD_nx_tpp, BD_ny_tpp,
          i, j, k, vz, BD_nz_vz, BD_nx_vz, 2, lambda, dz, dt,
          ptppbvz, bfull, afull, ext, fdc);
        }
        for(i=ext; i<BD_nz_tpp-ext; i++)
        {
          ac_body_tpp_z_3d(tpp, BD_nz_tpp, BD_nx_tpp, BD_ny_tpp,
          i, j, k, vz, BD_nz_vz, BD_nx_vz, 2, lambda, dz, dt, fdc);
        }
      }
    }
  }
}
