#include<stdio.h>
#include"2d_lib.c"

// vx
void el2d_openmp(double *vx, int BD_nx_vx, int BD_nz_vx,
  double *pvxbtxx, double *pvxbtxz,
  double *vz, int BD_nx_vz, int BD_nz_vz,
  double *pvzbtxz, double *pvzbtzz,
  double *txx, double *ptxxbvx, double *ptxxbvz,
  double *tzz, double *ptzzbvx, double *ptzzbvz,
  int BD_nx_tpp, int BD_nz_tpp,
  double *txz, int BD_nx_txz, int BD_nz_txz,
  double *ptxzbvx, double *ptxzbvz,
  double *rho, double *lambda, double *mu, double *fdc,
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
          txx, BD_nz_tpp, BD_nx_tpp, 1,
        2.0/(*(rho+i+j*BD_nz_tpp)+*(rho+i+(j+1)*BD_nz_tpp)),
        2.0/(*(rho+i+(BD_nx_vx-1-j)*BD_nz_tpp)+*(rho+i+(BD_nx_vx-j)*BD_nz_tpp)),
        dx, dt, pvxbtxx, bhalf, ahalf, ext, fdc);
      }
      for(j=ext; j<BD_nx_vx-ext; j++)
      {
        body_x_2d(vx, BD_nz_vx, BD_nx_vx, i, j,
          txx, BD_nz_tpp, BD_nx_tpp, 1,
        2.0/(*(rho+i+j*BD_nz_tpp)+*(rho+i+(j+1)*BD_nz_tpp)),
        dx, dt, fdc);
      }
    }

    // vzbtxz
    #pragma omp for private(j)
    for(i=0; i<BD_nz_vz; i++)
    {
      for(j=2; j<ext; j++)
      {
        bound_x_2d(vz, BD_nz_vz, BD_nx_vz, i, j,
          txz, BD_nz_txz, BD_nx_txz, 2,
        2.0/(*(rho+i+j*BD_nz_tpp)+*(rho+(i+1)+j*BD_nz_tpp)),
        2.0/(*(rho+i+(BD_nx_vz-1-j)*BD_nz_tpp)+*(rho+(i+1)+(BD_nx_vz-1-j)*BD_nz_tpp)),
        dx,dt, pvzbtxz, bfull, afull, ext, fdc);
      }
      for(j=ext; j<BD_nx_vz-ext; j++)
      {
        body_x_2d(vz, BD_nz_vz, BD_nx_vz, i, j,
          txz, BD_nz_txz, BD_nx_txz, 2,
        2.0/(*(rho+i+j*BD_nz_tpp)+*(rho+(i+1)+j*BD_nz_tpp)),
        dx, dt, fdc);
      }
    }

  //********************* V_Z *********************//
    // vxbtxz
    #pragma omp for private(i)
    for(j=1; j<BD_nx_vx-1; j++)
    {
      for(i=2; i<ext; i++)
      {
        unlimited_bound_z_2d(vx, BD_nz_vx, BD_nx_vx, i, j,
          txz, BD_nz_txz, BD_nx_txz, 2,
        2.0/(*(rho+i+j*BD_nz_tpp)+*(rho+i+(j+1)*BD_nz_tpp)),
        2.0/(*(rho+(BD_nz_vx-1-i)+j*BD_nz_tpp)+*(rho+(BD_nz_vx-1-i)+(j+1)*BD_nz_tpp)),
        dz,dt, pvxbtxz, bfull, afull, ext, fdc);
      }
      for(i=ext; i<BD_nz_vx-ext; i++)
      {
        body_z_2d(vx, BD_nz_vx, BD_nx_vx, i, j,
          txz, BD_nz_txz, BD_nx_txz, 2,
          2.0/(*(rho+i+j*BD_nz_tpp)+*(rho+i+(j+1)*BD_nz_tpp)),
          dz,dt, fdc);
      }
    }

    // vzbtzz
    #pragma omp for private(i)
    for(j=0; j<BD_nx_vz; j++)
    {
      for(i=1;i<ext;i++)
      {
        unlimited_bound_z_2d(vz, BD_nz_vz, BD_nx_vz, i, j,
          tzz, BD_nz_tpp, BD_nx_tpp, 1,
          2.0/(*(rho+i+j*BD_nz_tpp)+*(rho+(i+1)+j*BD_nz_tpp)),
          2.0/(*(rho+(BD_nz_vz-1-i)+j*BD_nz_tpp)+*(rho+(BD_nz_vz-i)+j*BD_nz_tpp)),
          dz,dt, pvzbtzz, bhalf, ahalf, ext, fdc);
      }
      for(i=ext; i<BD_nz_vz-ext; i++)
      {
        body_z_2d(vz, BD_nz_vz, BD_nx_vz, i, j,
          tzz, BD_nz_tpp, BD_nx_tpp, 1,
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
        el_bound_tpp_x_2d(txx, tzz, BD_nz_tpp, BD_nx_tpp,
        i, j, vx, BD_nz_vx, BD_nx_vx, 2, lambda, mu, dx, dt,
        ptxxbvx, ptzzbvx, bfull, afull, ext, fdc);
      }
      for(j=ext; j<BD_nx_tpp-ext; j++)
      {
        el_body_tpp_x_2d(txx, tzz, BD_nz_tpp, BD_nx_tpp,
        i, j, vx, BD_nz_vx, BD_nx_vx, 2, lambda, mu, dx, dt, fdc);
      }
    }

    #pragma omp for private(j)
    for(i=0; i<BD_nz_txz; i++)
    {
      for(j=1; j<ext; j++)
      {
        bound_x_2d(txz, BD_nz_txz, BD_nx_txz, i, j,
        vz, BD_nz_vz, BD_nx_vz, 1,
        (*(mu+i+j*BD_nz_tpp) + *(mu+i+(j+1)*BD_nz_tpp)
        + *(mu+(i+1)+(j+1)*BD_nz_tpp) + *(mu+(i+1)+j*BD_nz_tpp))/4,
        (*(mu+i+(BD_nx_txz-1-j)*BD_nz_tpp) + *(mu+i+(BD_nx_txz-j)*BD_nz_tpp)
        + *(mu+(i+1)+(BD_nx_txz-j)*BD_nz_tpp) + *(mu+(i+1)+(BD_nx_txz-1-j)*BD_nz_tpp))/4,
        dx, dt, ptxzbvz, bhalf, ahalf, ext, fdc);
      }
      for(j=ext; j<BD_nx_txz-ext; j++)
      {
        body_x_2d(txz, BD_nz_txz, BD_nx_txz, i, j,
        vz, BD_nz_vz, BD_nx_vz, 1,
        (*(mu+i+j*BD_nz_tpp) + *(mu+i+(j+1)*BD_nz_tpp)
        + *(mu+(i+1)+(j+1)*BD_nz_tpp) + *(mu+(i+1)+j*BD_nz_tpp))/4,
        dx, dt, fdc);
      }
    }

  //************ T_Z ************//
    #pragma omp for private(i)
    for(j=0; j<BD_nx_tpp; j++)
    {
      for(i=2; i<ext; i++)
      {
        el_unlimited_bound_tpp_z_2d(txx, tzz, BD_nz_tpp, BD_nx_tpp,
        i, j, vz, BD_nz_vz, BD_nx_vz, 2, lambda, mu, dz, dt,
        ptxxbvz, ptzzbvz, bfull, afull, ext, fdc);
      }
      for(i=ext; i<BD_nz_tpp-ext; i++)
      {
        el_body_tpp_z_2d(txx, tzz, BD_nz_tpp, BD_nx_tpp,
        i, j, vz, BD_nz_vz, BD_nx_vz, 2, lambda, mu, dz, dt, fdc);
      }
    }
    #pragma omp for private(i)
    for(j=1; j<BD_nx_txz-1; j++)
    {
      for(i=1; i<ext; i++)
      {
        unlimited_bound_z_2d(txz, BD_nz_txz, BD_nx_txz, i, j,
        vx, BD_nz_vx, BD_nx_vx, 1,
        (*(mu+i+j*BD_nz_tpp) + *(mu+i+(j+1)*BD_nz_tpp)
        + *(mu+(i+1)+(j+1)*BD_nz_tpp) + *(mu+(i+1)+j*BD_nz_tpp))/4,
        (*(mu+(BD_nz_txz-1-i)+j*BD_nz_tpp) + *(mu+(BD_nz_txz-1-i)+(j+1)*BD_nz_tpp)
        + *(mu+(BD_nz_txz-i)+(j+1)*BD_nz_tpp) + *(mu+(BD_nz_txz-i)+j*BD_nz_tpp))/4,
        dz, dt, ptxzbvx, bhalf, ahalf, ext, fdc);
      }
      for(i=ext; i<BD_nz_txz-ext; i++)
      {
        body_z_2d(txz, BD_nz_txz, BD_nx_txz, i, j,
        vx, BD_nz_vx, BD_nx_vx, 1,
        (*(mu+i+j*BD_nz_tpp) + *(mu+i+(j+1)*BD_nz_tpp)
        + *(mu+(i+1)+(j+1)*BD_nz_tpp) + *(mu+(i+1)+j*BD_nz_tpp))/4,
        dz, dt, fdc);
      }
    }
  }
}
