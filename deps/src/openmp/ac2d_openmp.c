#include<stdio.h>
#include"2d_lib.c"

// vx
void ac2d_openmp(double *vx, int BD_nx_vx, int BD_nz_vx,
  double *pvxbtpp,
  double *vz, int BD_nx_vz, int BD_nz_vz,
  double *pvzbtpp,
  double *tpp, double *ptppbvx, double *ptppbvz,
  int BD_nx_tpp, int BD_nz_tpp,
  double *rho, double *lambda, double *fdc,
  double dt, double dx, double dz, int ext,
  double *bhalf, double *ahalf, double *bfull, double *afull)
{
  int i,j;

  //********************* V_X *********************//
  #pragma omp parallel
  {
    // vxbtxx
    #pragma omp for private(j)
      for(i=2; i<BD_nz_vx-2; i++)
      {
        for(j=1; j<ext; j++)
        {
          bound_x_2d(vx, BD_nz_vx, BD_nx_vx, i, j,
            tpp, BD_nz_tpp, BD_nx_tpp, 1,
          2.0/(*(rho+i+j*BD_nz_tpp)+*(rho+i+(j+1)*BD_nz_tpp)),
          2.0/(*(rho+i+(BD_nx_vx-1-j)*BD_nz_tpp)+*(rho+i+(BD_nx_vx-j)*BD_nz_tpp)),
          dx, dt, pvxbtpp, bhalf, ahalf, ext, fdc);
        }
        for(j=ext; j<BD_nx_vx-ext; j++)
        {
          body_x_2d(vx, BD_nz_vx, BD_nx_vx, i, j,
            tpp, BD_nz_tpp, BD_nx_tpp, 1,
          2.0/(*(rho+i+j*BD_nz_tpp)+*(rho+i+(j+1)*BD_nz_tpp)),
          dx, dt, fdc);
        }
      }

  //********************* V_Z *********************//
    // vzbtzz
    #pragma omp for private(i)
    for(j=0; j<BD_nx_vz; j++)
    {
      for(i=1;i<ext;i++)
      {
        unlimited_bound_z_2d(vz, BD_nz_vz, BD_nx_vz, i, j,
          tpp, BD_nz_tpp, BD_nx_tpp, 1,
          2.0/(*(rho+i+j*BD_nz_tpp)+*(rho+(i+1)+j*BD_nz_tpp)),
          2.0/(*(rho+(BD_nz_vz-1-i)+j*BD_nz_tpp)+*(rho+(BD_nz_vz-i)+j*BD_nz_tpp)),
          dz,dt, pvzbtpp, bhalf, ahalf, ext, fdc);
      }
      for(i=ext; i<BD_nz_vz-ext; i++)
      {
        body_z_2d(vz, BD_nz_vz, BD_nx_vz, i, j,
          tpp, BD_nz_tpp, BD_nx_tpp, 1,
          2.0/(*(rho+i+j*BD_nz_tpp)+*(rho+(i+1)+j*BD_nz_tpp)),
          dz,dt, fdc);
      }
    }

  //************ T_X ************//
    #pragma omp for private(j)
    for(i=2; i<BD_nz_tpp-2; i++)
    {
      for(j=2; j<ext; j++)
      {
        ac_bound_tpp_x_2d(tpp, BD_nz_tpp, BD_nx_tpp,
        i, j, vx, BD_nz_vx, BD_nx_vx, 2, lambda, dx, dt,
        ptppbvx, bfull, afull, ext, fdc);
      }
      for(j=ext; j<BD_nx_tpp-ext; j++)
      {
        ac_body_tpp_x_2d(tpp, BD_nz_tpp, BD_nx_tpp,
        i, j, vx, BD_nz_vx, BD_nx_vx, 2, lambda, dx, dt, fdc);
      }
    }

  //************ T_Z ************//
    #pragma omp for private(i)
    for(j=0; j<BD_nx_tpp; j++)
    {
      for(i=2; i<ext; i++)
      {
        ac_unlimited_bound_tpp_z_2d(tpp, BD_nz_tpp, BD_nx_tpp,
        i, j, vz, BD_nz_vz, BD_nx_vz, 2, lambda, dz, dt,
        ptppbvz, bfull, afull, ext, fdc);
      }
      for(i=ext; i<BD_nz_tpp-ext; i++)
      {
        ac_body_tpp_z_2d(tpp, BD_nz_tpp, BD_nx_tpp,
        i, j, vz, BD_nz_vz, BD_nx_vz, 2, lambda, dz, dt, fdc);
      }
    }
  }
}
