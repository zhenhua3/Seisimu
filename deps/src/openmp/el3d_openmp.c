#include<stdio.h>
#include"3d_lib.c"

// vx
void el3d_openmp(double *vx, int BD_nx_vx, int BD_ny_vx, int BD_nz_vx,
  double *pvxbtxx, double *pvxbtxy, double *pvxbtxz,
  double *vy, int BD_nx_vy, int BD_ny_vy, int BD_nz_vy,
  double *pvybtxy, double *pvybtyy, double *pvybtyz,
  double *vz, int BD_nx_vz, int BD_ny_vz, int BD_nz_vz,
  double *pvzbtxz, double *pvzbtyz, double *pvzbtzz,
  double *txx, double *ptxxbvx, double *ptxxbvy, double *ptxxbvz,
  double *tyy, double *ptyybvx, double *ptyybvy, double *ptyybvz,
  double *tzz, double *ptzzbvx, double *ptzzbvy, double *ptzzbvz,
  int BD_nx_tpp, int BD_ny_tpp, int BD_nz_tpp,
  double *txy, int BD_nx_txy, int BD_ny_txy, int BD_nz_txy,
  double *ptxybvx, double *ptxybvy,
  double *tyz, int BD_nx_tyz, int BD_ny_tyz, int BD_nz_tyz,
  double *ptyzbvy, double *ptyzbvz,
  double *txz, int BD_nx_txz, int BD_ny_txz, int BD_nz_txz,
  double *ptxzbvx, double *ptxzbvz,
  double *rho, double *lambda, double *mu, double *fdc,
  double dt, double dx, double dy, double dz, int ext,
  double *bhalf, double *ahalf, double *bfull, double *afull)
{
  int i,j,k,t;

  for(t = 0; t<1; t++){
    //********************* V_X *********************//
    #pragma omp parallel
    {
      // vxbtxx
      #pragma omp for collapse(2) private(j)
      for(k=2; k<BD_ny_vx-2; k++)
      {
        for(i=2; i<BD_nz_vx-2; i++)
        {
          for(j=1; j<ext; j++)
          {
            bound_x_3d(vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              txx, BD_nz_tpp, BD_nx_tpp, 1,
            2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            2.0/(*(rho+i+(BD_nx_vx-1-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(BD_nx_vx-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            dx, dt, pvxbtxx, bhalf, ahalf, ext, fdc);
          }
          for(j=ext; j<BD_nx_vx-ext; j++)
          {
            body_x_3d(vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              txx, BD_nz_tpp, BD_nx_tpp, 1,
            2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            dx, dt, fdc);
          }
        }
      }
      // vybtxy
      #pragma omp for collapse(2) private(j)
      for(k=1; k<BD_ny_vy-1; k++)
      {
        for(i=2; i<BD_nz_vy-2; i++)
        {
          for(j=2; j<ext; j++)
          {
            bound_x_3d(vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              txy, BD_nz_txy, BD_nx_txy, 2,
            2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
            2.0/(*(rho+i+(BD_nx_vy-1-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(BD_nx_vy-1-j)*BD_nz_tpp+(BD_ny_vy-k)*BD_nz_tpp*BD_nx_tpp)),
            dx, dt, pvybtxy, bfull, afull, ext, fdc);
          }
          for(j=ext; j<BD_nx_vy-ext; j++)
          {
            body_x_3d(vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              txy, BD_nz_txy, BD_nx_txy, 2,
            2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
            dx, dt, fdc);
          }
        }
      }
      // vzbtxz
      #pragma omp for collapse(2) private(j)
      for(k=0; k<BD_ny_vz; k++)
      {
        for(i=0; i<BD_nz_vz; i++)
        {
          for(j=2; j<ext; j++)
          {
            bound_x_3d(vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              txz, BD_nz_txz, BD_nx_txz, 2,
            2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            2.0/(*(rho+i+(BD_nx_vz-1-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(i+1)+(BD_nx_vz-1-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            dx,dt, pvzbtxz, bfull, afull, ext, fdc);
          }
          for(j=ext; j<BD_nx_vz-ext; j++)
          {
            body_x_3d(vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              txz, BD_nz_txz, BD_nx_txz, 2,
            2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            dx, dt, fdc);
          }
        }
      }
    //********************* V_Y *********************//
      // vxbtxy
      #pragma omp for collapse(2) private(k)
      for(j=0; j<BD_nx_vx; j++)
      {
        for(i=0; i<BD_nz_vx; i++)
        {
          for(k=2; k<ext; k++)
          {
            bound_y_3d(vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              txy, BD_nz_txy, BD_nx_txy, 2,
              2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              2.0/(*(rho+i+j*BD_nz_tpp+(BD_ny_vx-1-k)*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(j+1)*BD_nz_tpp+(BD_ny_vx-1-k)*BD_nz_tpp*BD_nx_tpp)),
              dy,dt, pvxbtxy, bfull, afull, ext, fdc);
          }
          for(k=ext; k<BD_ny_vx-ext; k++)
          {
            body_y_3d(vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              txy, BD_nz_txy, BD_nx_txy, 2,
              2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dy,dt, fdc);
          }
        }
      }
      // vybtyy
      #pragma omp for collapse(2) private(k)
      for(j=0; j<BD_nx_vy; j++)
      {
        for(i=0; i<BD_nz_vy; i++)
        {
          for(k=1; k<ext; k++)
          {
            bound_y_3d(vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              tyy, BD_nz_tpp, BD_nx_tpp, 1,
              2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              2.0/(*(rho+i+j*BD_nz_tpp+(BD_ny_vy-1-k)*BD_nz_tpp*BD_nx_tpp)+*(rho+i+j*BD_nz_tpp+(BD_ny_vy-k)*BD_nz_tpp*BD_nx_tpp)),
              dy,dt, pvybtyy, bhalf, ahalf, ext, fdc);
          }
          for(k=ext; k<BD_ny_vy-ext; k++)
          {
            body_y_3d(vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              tyy, BD_nz_tpp, BD_nx_tpp, 1,
              2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              dy,dt, fdc);
          }
        }
      }
      // vzbtyz
      #pragma omp for collapse(2) private(k)
      for(j=0; j<BD_nx_vz; j++)
      {
        for(i=0; i<BD_nz_vz; i++)
        {
          for(k=2; k<ext; k++)
          {
            bound_y_3d(vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              tyz, BD_nz_tyz, BD_nx_tyz, 2,
              2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              2.0/(*(rho+i+j*BD_nz_tpp+(BD_ny_vz-1-k)*BD_nz_tpp*BD_nx_tpp)+*(rho+(i+1)+j*BD_nz_tpp+(BD_ny_vz-1-k)*BD_nz_tpp*BD_nx_tpp)),
              dy,dt, pvzbtyz, bfull, afull, ext, fdc);
          }
          for(k=ext; k<BD_ny_vz-ext; k++)
          {
            body_y_3d(vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              tyz, BD_nz_tyz, BD_nx_tyz, 2,
            2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            dy,dt, fdc);
          }
        }
      }

    //********************* V_Z *********************//
      // vxbtxz
      #pragma omp for collapse(2) private(i)
      for(j=1; j<BD_nx_vx-1; j++)
      {
        for(k=2; k<BD_ny_vx-2; k++)
        {
          for(i=2; i<ext; i++)
          {
            unlimited_bound_z_3d(vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              txz, BD_nz_txz, BD_nx_txz, 2,
            2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            2.0/(*(rho+(BD_nz_vx-1-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(BD_nz_vx-1-i)+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            dz,dt, pvxbtxz, bfull, afull, ext, fdc);
          }
          for(i=ext; i<BD_nz_vx-ext; i++)
          {
            body_z_3d(vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              txz, BD_nz_txz, BD_nx_txz, 2,
              2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dz,dt, fdc);
          }
        }
      }
      // vybtyz
      #pragma omp for collapse(2) private(i)
      for(j=0; j<BD_nx_vy; j++)
      {
        for(k=0; k<BD_ny_vy; k++)
        {
          for(i=2; i<ext; i++)
          {
            unlimited_bound_z_3d(vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              tyz, BD_nz_tyz, BD_nx_tyz, 2,
              2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              2.0/(*(rho+(BD_nz_vy-1-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(BD_nz_vy-1-i)+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              dz,dt, pvybtyz, bfull, afull, ext, fdc);
          }
          for(i=ext; i<BD_nz_vy-ext; i++)
          {
            body_z_3d(vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              tyz, BD_nz_tyz, BD_nx_tyz, 2,
              2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              dz,dt, fdc);
          }
        }
      }
      // vzbtzz
      #pragma omp for collapse(2) private(i)
      for(j=0; j<BD_nx_vz; j++)
      {
        for(k=0; k<BD_ny_vz; k++)
        {
          for(i=1;i<ext;i++)
          {
            unlimited_bound_z_3d(vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              tzz, BD_nz_tpp, BD_nx_tpp, 1,
              2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              2.0/(*(rho+(BD_nz_vz-1-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(BD_nz_vz-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dz,dt, pvzbtzz, bhalf, ahalf, ext, fdc);
          }
          for(i=ext; i<BD_nz_vz-ext; i++)
          {
            body_z_3d(vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              tzz, BD_nz_tpp, BD_nx_tpp, 1,
              2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dz,dt, fdc);
          }
        }
      }

    //************ T_X ************//
      #pragma omp for collapse(2) private(j)
      for(k=2; k<BD_ny_tpp-2; k++)
      {
        for(i=2; i<BD_nz_tpp-2; i++)
        {
          for(j=2; j<ext; j++)
          {
            el_bound_tpp_x_3d(txx, tyy, tzz, BD_nz_tpp, BD_nx_tpp, BD_ny_tpp,
            i, j, k, vx, BD_nz_vx, BD_nx_vx, 2, lambda, mu, dx, dt,
            ptxxbvx, ptyybvx, ptzzbvx, bfull, afull, ext, fdc);
          }
          for(j=ext; j<BD_nx_tpp-ext; j++)
          {
            el_body_tpp_x_3d(txx, tyy, tzz, BD_nz_tpp, BD_nx_tpp, BD_ny_tpp,
            i, j, k, vx, BD_nz_vx, BD_nx_vx, 2, lambda, mu, dx, dt, fdc);
          }
        }
      }
      #pragma omp for collapse(2) private(j)
      for(k=1; k<BD_ny_txy-1; k++)
      {
        for(i=0; i<BD_nz_txy; i++)
        {
          for(j=1; j<ext; j++)
          {
            bound_x_3d(txy, BD_nz_txy, BD_nx_txy, BD_ny_txy, i, j, k,
            vy, BD_nz_vy, BD_nx_vy, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+i+(j+1)*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp) + *(mu+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp))/4,
            (*(mu+i+(BD_nx_txy-1-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(BD_nx_txy-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+i+(BD_nx_txy-j)*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(BD_nx_txy-1-j)*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp))/4,
            dx, dt, ptxybvy, bhalf, ahalf, ext, fdc);
          }
          for(j=ext; j<BD_nx_txy-ext; j++)
          {
            body_x_3d(txy, BD_nz_txy, BD_nx_txy, BD_ny_txy, i, j, k,
            vy, BD_nz_vy, BD_nx_vy, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+i+(j+1)*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp) + *(mu+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp))/4,
            dx, dt, fdc);
          }
        }
      }
      #pragma omp for collapse(2) private(j)
      for(k=0; k<BD_ny_txz; k++)
      {
        for(i=0; i<BD_nz_txz; i++)
        {
          for(j=1; j<ext; j++)
          {
            bound_x_3d(txz, BD_nz_txz, BD_nx_txz, BD_ny_txz, i, j, k,
            vz, BD_nz_vz, BD_nx_vz, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(i+1)+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            (*(mu+i+(BD_nx_txz-1-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(BD_nx_txz-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(i+1)+(BD_nx_txz-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+(i+1)+(BD_nx_txz-1-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            dx, dt, ptxzbvz, bhalf, ahalf, ext, fdc);
          }
          for(j=ext; j<BD_nx_txz-ext; j++)
          {
            body_x_3d(txz, BD_nz_txz, BD_nx_txz, BD_ny_txz, i, j, k,
            vz, BD_nz_vz, BD_nx_vz, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(i+1)+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            dx, dt, fdc);
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
            el_bound_tpp_y_3d(txx, tyy, tzz, BD_nz_tpp, BD_nx_tpp, BD_ny_tpp,
            i, j, k, vy, BD_nz_vy, BD_nx_vy, 2, lambda, mu, dy, dt,
            ptxxbvy, ptyybvy, ptzzbvy, bfull, afull, ext, fdc);
          }
          for(k=ext; k<BD_ny_tpp-ext; k++)
          {
            el_body_tpp_y_3d(txx, tyy, tzz, BD_nz_tpp, BD_nx_tpp, BD_ny_tpp,
            i, j, k, vy, BD_nz_vy, BD_nx_vy, 2, lambda, mu, dy, dt, fdc);
          }
        }
      }
      #pragma omp for collapse(2) private(k)
      for(i=0; i<BD_nz_txy; i++)
      {
        for(j=1; j<BD_nx_txy-1; j++)
        {
          for(k=1; k<ext; k++)
          {
            bound_y_3d(txy, BD_nz_txy, BD_nx_txy, BD_ny_txy, i, j, k,
            vx, BD_nz_vx, BD_nx_vx, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+i+(j+1)*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp) + *(mu+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp))/4,
            (*(mu+i+j*BD_nz_tpp+(BD_ny_txy-k-1)*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(j+1)*BD_nz_tpp+(BD_ny_txy-k-1)*BD_nz_tpp*BD_nx_tpp)
            + *(mu+i+(j+1)*BD_nz_tpp+(BD_ny_txy-k)*BD_nz_tpp*BD_nx_tpp) + *(mu+i+j*BD_nz_tpp+(BD_ny_txy-k)*BD_nz_tpp*BD_nx_tpp))/4,
            dy, dt, ptxybvx, bhalf, ahalf, ext, fdc);
          }
          for(k=ext; k<BD_ny_txy-ext; k++)
          {
            body_y_3d(txy, BD_nz_txy, BD_nx_txy, BD_ny_txy, i, j, k,
            vx, BD_nz_vx, BD_nx_vx, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+i+(j+1)*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp) + *(mu+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp))/4,
            dy, dt, fdc);
          }
        }
      }
      #pragma omp for collapse(2) private(k)
      for(i=0; i<BD_nz_tyz; i++)
      {
        for(j=0; j<BD_nx_tyz; j++)
        {
          for(k=1; k<ext; k++)
          {
            bound_y_3d(tyz, BD_nz_tyz, BD_nx_tyz, BD_ny_tyz, i, j, k,
            vz, BD_nz_vz, BD_nx_vz, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(i+1)+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp) + *(mu+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            (*(mu+i+j*BD_nz_tpp+(BD_ny_tyz-1-k)*BD_nz_tpp*BD_nx_tpp) + *(mu+i+j*BD_nz_tpp+(BD_ny_tyz-k)*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(i+1)+j*BD_nz_tpp+(BD_ny_tyz-k)*BD_nz_tpp*BD_nx_tpp) + *(mu+(i+1)+j*BD_nz_tpp+(BD_ny_tyz-1-k)*BD_nz_tpp*BD_nx_tpp))/4,
            dy, dt, ptyzbvz, bhalf, ahalf, ext, fdc);
          }
          for(k=ext; k<BD_ny_tyz-ext; k++)
          {
            body_y_3d(tyz, BD_nz_tyz, BD_nx_tyz, BD_ny_tyz, i, j, k,
            vz, BD_nz_vz, BD_nx_vz, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(i+1)+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp) + *(mu+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            dy, dt, fdc);
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
            el_unlimited_bound_tpp_z_3d(txx, tyy, tzz, BD_nz_tpp, BD_nx_tpp, BD_ny_tpp,
            i, j, k, vz, BD_nz_vz, BD_nx_vz, 2, lambda, mu, dz, dt,
            ptxxbvz, ptyybvz, ptzzbvz, bfull, afull, ext, fdc);
          }
          for(i=ext; i<BD_nz_tpp-ext; i++)
          {
            el_body_tpp_z_3d(txx, tyy, tzz, BD_nz_tpp, BD_nx_tpp, BD_ny_tpp,
            i, j, k, vz, BD_nz_vz, BD_nx_vz, 2, lambda, mu, dz, dt, fdc);
          }
        }
      }
      #pragma omp for collapse(2) private(i)
      for(j=1; j<BD_nx_txz-1; j++)
      {
        for(k=0; k<BD_ny_txz; k++)
        {
          for(i=1; i<ext; i++)
          {
            unlimited_bound_z_3d(txz, BD_nz_txz, BD_nx_txz, BD_ny_txz, i, j, k,
            vx, BD_nz_vx, BD_nx_vx, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(i+1)+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            (*(mu+(BD_nz_txz-1-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+(BD_nz_txz-1-i)+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(BD_nz_txz-i)+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+(BD_nz_txz-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            dz, dt, ptxzbvx, bhalf, ahalf, ext, fdc);
          }
          for(i=ext; i<BD_nz_txz-ext; i++)
          {
            body_z_3d(txz, BD_nz_txz, BD_nx_txz, BD_ny_txz, i, j, k,
            vx, BD_nz_vx, BD_nx_vx, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(i+1)+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            dz, dt, fdc);
          }
        }
      }
      #pragma omp for collapse(2) private(i)
      for(j=0; j<BD_nx_tyz; j++)
      {
        for(k=0; k<BD_ny_tyz; k++)
        {
          for(i=1; i<ext; i++)
          {
            unlimited_bound_z_3d(tyz, BD_nz_tyz, BD_nx_tyz, BD_ny_tyz, i, j, k,
            vy, BD_nz_vy, BD_nx_vy, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(i+1)+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp) + *(mu+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            (*(mu+(BD_nz_tyz-1-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+(BD_nz_tyz-1-i)+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(BD_nz_tyz-i)+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp) + *(mu+(BD_nz_tyz-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            dz, dt, ptyzbvy, bhalf, ahalf, ext, fdc);
          }
          for(i=ext; i<BD_nz_tyz-ext; i++)
          {
            body_z_3d(tyz, BD_nz_tyz, BD_nx_tyz, BD_ny_tyz, i, j, k,
            vy, BD_nz_vy, BD_nx_vy, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(i+1)+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp) + *(mu+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            dz, dt, fdc);
          }
        }
      }
    }
  }
}
