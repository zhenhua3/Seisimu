#include<stdio.h>

// vx
void elextrap2domp(double *vx, double *tmp_vx,  int nvx1, int nvx2, double *rho_vx, double *pvxbtxx, double *pvxbtxz, double *pvxbtxy,
  double *vz, double *tmp_vz,  int nvz1, int nvz2, double *rho_vz, double *pvzbtzz, double *pvzbtxz,
  double *txx, double *tzz, double *tmp_tpp, int ntpp1, int ntpp2, double *lamu, double *lambda, double *mu,
  double *ptxxbvx, double *ptxxbvz, double *ptzzbvx, double *ptzzbvz,
  double *txz, double *tmp_txz, int ntxz1, int ntxz2, double *mu_txz, double *ptxzbvx, double *ptxzbvz,
  double *fdc, double dt, double dx, double dz, int ext,
  double *bhalf, double *ahalf, double *bfull, double *afull)
{
  int i,j;
  // vxbtxx
  #pragma omp parallel for collapse(2)
    for(i=0; i<nvx1+2*ext; i++)
    {
      for(j=1; j<ext; j++)
      {
        *(tmp_vx+i+j*(nvx1+2*ext)) = (*(txx+i+(j-1)*(ntpp1+2*ext))* *(fdc) + *(txx+i+j*(ntpp1+2*ext))* *(fdc+1)
        + *(txx+i+(j+1)*(ntpp1+2*ext))* *(fdc+2) + *(txx+i+(j+2)*(ntpp1+2*ext))* *(fdc+3));
        *(tmp_vx+i+j*(nvx1+2*ext)) = *(rho_vx+i+j*(nvx1+2*ext))* *(tmp_vx+i+j*(nvx1+2*ext))/dx;
        *(pvxbtxx+i+j*(nvx1+2*ext)) = *(bhalf+j) * *(pvxbtxx+i+j*(nvx1+2*ext)) + *(ahalf+j)* *(tmp_vx+i+j*(nvx1+2*ext));
        *(vx+i+j*(nvx1+2*ext)) = *(vx+i+j*(nvx1+2*ext)) + *(tmp_vx+i+j*(nvx1+2*ext)) + *(pvxbtxx+i+j*(nvx1+2*ext));
      }
    }
    #pragma omp parallel for collapse(2)
    for(i=0; i<nvx1+2*ext; i++)
    {
      for(j=ext; j<ext+nvx2; j++)
      {
      *(tmp_vx+i+j*(nvx1+2*ext)) = (*(txx+i+(j-1)*(ntpp1+2*ext))* *(fdc) + *(txx+i+j*(ntpp1+2*ext))* *(fdc+1)
      + *(txx+i+(j+1)*(ntpp1+2*ext))* *(fdc+2) + *(txx+i+(j+2)*(ntpp1+2*ext))* *(fdc+3));
      *(vx+i+j*(nvx1+2*ext)) = *(vx+i+j*(nvx1+2*ext)) + *(rho_vx+i+j*(nvx1+2*ext))**(tmp_vx+i+j*(nvx1+2*ext))/dx;
      }
    }
    #pragma omp parallel for collapse(2)
    for(i=0; i<nvx1+2*ext; i++)
    {
      for(j=ext+nvx2; j<nvx2+2*ext-1; j++)
      {
      *(tmp_vx+i+j*(nvx1+2*ext)) = (*(txx+i+(j-1)*(ntpp1+2*ext))* *(fdc) + *(txx+i+j*(ntpp1+2*ext))* *(fdc+1)
      + *(txx+(j+1)*(ntpp1+2*ext)+i)* *(fdc+2) + *(txx+(j+2)*(ntpp1+2*ext)+i)* *(fdc+3));
      *(tmp_vx+i+j*(nvx1+2*ext)) = *(rho_vx+i+j*(nvx1+2*ext))**(tmp_vx+i+j*(nvx1+2*ext))/dx;
     *(pvxbtxx+i+(j-nvx2)*(nvx1+2*ext)) = (*(bhalf+nvx2+2*ext-j-1) * *(pvxbtxx+i+(j-nvx2)*(nvx1+2*ext))
     + *(ahalf+nvx2+2*ext-j-1)* *(tmp_vx+i+j*(nvx1+2*ext)));
      *(vx+i+j*(nvx1+2*ext)) = *(vx+i+j*(nvx1+2*ext)) + *(tmp_vx+i+j*(nvx1+2*ext)) + *(pvxbtxx+i+(j-nvx2)*(nvx1+2*ext));
      }
    }

  // vxbtxz

  #pragma omp parallel for collapse(2)
  for(j=0; j<nvx2+2*ext; j++)
  {
    for(i=2; i<ext; i++)
    {
      *(tmp_vx+i+j*(nvx1+2*ext)) = (*(txz+i-2+j*(2*ext+ntxz1))* *(fdc) + *(txz+i-1+j*(2*ext+ntxz1))* *(fdc+1)
      + *(txz+i+j*(2*ext+ntxz1))* *(fdc+2) + *(txz+i+1+j*(2*ext+ntxz1))* *(fdc+3));
      *(tmp_vx+i+j*(nvx1+2*ext)) = *(rho_vx+i+j*(nvx1+2*ext))* *(tmp_vx+i+j*(nvx1+2*ext))/dz;
      *(pvxbtxz+i+j*2*ext) = *(bfull+i) * *(pvxbtxz+i+j*2*ext) + *(afull+i)* *(tmp_vx+i+j*(nvx1+2*ext));
      *(vx+i+j*(nvx1+2*ext)) = *(vx+i+j*(nvx1+2*ext)) + *(tmp_vx+i+j*(nvx1+2*ext)) + *(pvxbtxz+i+j*2*ext);
    }
  }
  #pragma omp parallel for collapse(2)
  for(j=0; j<nvx2+2*ext; j++)
  {
    for(i=ext; i<ext+nvx1; i++)
    {
      *(tmp_vx+i+j*(nvx1+2*ext)) = (*(txz+i-2+j*(2*ext+ntxz1))* *(fdc) + *(txz+i-1+j*(2*ext+ntxz1))* *(fdc+1)
      + *(txz+i+j*(2*ext+ntxz1))* *(fdc+2) + *(txz+i+1+j*(2*ext+ntxz1))* *(fdc+3));
      *(vx+i+j*(nvx1+2*ext)) = *(vx+i+j*(nvx1+2*ext)) + *(rho_vx+i+j*(nvx1+2*ext))* *(tmp_vx+i+j*(nvx1+2*ext))/dz;
    }
  }
  #pragma omp parallel for collapse(2)
  for(j=0; j<nvx2+2*ext; j++)
  {
    for(i=ext+nvx1; i<nvx1+2*ext-2; i++)
    {
      *(tmp_vx+i+j*(nvx1+2*ext)) = (*(txz+i-2+j*(2*ext+ntxz1))* *(fdc) + *(txz+i-1+j*(2*ext+ntxz1))* *(fdc+1)
      + *(txz+i+j*(2*ext+ntxz1))* *(fdc+2) + *(txz+i+1+j*(2*ext+ntxz1))* *(fdc+3));
      *(tmp_vx+i+j*(nvx1+2*ext)) = *(rho_vx+i+j*(nvx1+2*ext))* *(tmp_vx+i+j*(nvx1+2*ext))/dz;
      *(pvxbtxz+i-nvx1+j*2*ext) = (*(bfull+nvx1+2*ext-i-1) * *(pvxbtxz+i-nvx1+j*2*ext)
      + *(afull+nvx1+2*ext-i-1)* *(tmp_vx+i+j*(nvx1+2*ext)));
      *(vx+i+j*(nvx1+2*ext)) = *(vx+i+j*(nvx1+2*ext)) + *(tmp_vx+i+j*(nvx1+2*ext)) + *(pvxbtxz+i-nvx1+j*2*ext) ;
    }
  }

  // vzbtzz
  #pragma omp parallel for collapse(2)
  for(j=0; j<nvz2+2*ext; j++)
  {
    for(i=1; i<ext; i++)
    {
      *(tmp_vz+i+j*(nvz1+2*ext)) = (*(tzz+i-1+j*(ntpp1+2*ext))* *(fdc) + *(tzz+i+j*(ntpp1+2*ext))* *(fdc+1)
      + *(tzz+i+1+j*(ntpp1+2*ext))* *(fdc+2) + *(tzz+i+2+j*(ntpp1+2*ext))* *(fdc+3));
      *(tmp_vz+i+j*(nvz1+2*ext)) = *(rho_vz+i+j*(nvz1+2*ext))* *(tmp_vz+i+j*(nvz1+2*ext))/dz;
      *(pvzbtzz+i+j*2*ext) = *(bhalf+i) * *(pvzbtzz+i+j*2*ext) + *(ahalf+i)* *(tmp_vz+i+j*(nvz1+2*ext));
      *(vz+i+j*(nvz1+2*ext)) = *(vz+i+j*(nvz1+2*ext)) + *(tmp_vz+i+j*(nvz1+2*ext)) + *(pvzbtzz+i+j*2*ext);
    }
  }
  #pragma omp parallel for collapse(2)
  for(j=0; j<nvz2+2*ext; j++)
  {
    for(i=ext; i<ext+nvz1; i++)
    {
      *(tmp_vz+i+j*(nvz1+2*ext)) = (*(tzz+i-1+j*(ntpp1+2*ext))* *(fdc) + *(tzz+i+j*(ntpp1+2*ext))* *(fdc+1)
      + *(tzz+i+1+j*(ntpp1+2*ext))* *(fdc+2) + *(tzz+i+2+j*(ntpp1+2*ext))* *(fdc+3));
      *(vz+i+j*(nvz1+2*ext)) = *(vz+i+j*(nvz1+2*ext)) + *(rho_vz+i+j*(nvz1+2*ext))* *(tmp_vz+i+j*(nvz1+2*ext))/dz;
    }
  }
  #pragma omp parallel for collapse(2)
  for(j=0; j<nvz2+2*ext; j++)
  {
    for(i=ext+nvz1; i<nvz1+2*ext-1; i++)
    {
      *(tmp_vz+i+j*(nvz1+2*ext)) = (*(tzz+i-1+j*(ntpp1+2*ext))* *(fdc) + *(tzz+i+j*(ntpp1+2*ext))* *(fdc+1)
      + *(tzz+i+1+j*(ntpp1+2*ext))* *(fdc+2) + *(tzz+i+2+j*(ntpp1+2*ext))* *(fdc+3));
      *(tmp_vz+i+j*(nvz1+2*ext)) = *(rho_vz+i+j*(nvz1+2*ext))* *(tmp_vz+i+j*(nvz1+2*ext))/dz;
      *(pvzbtzz+i-nvz1+j*2*ext) = (*(bhalf+nvz1+2*ext-i-1) * *(pvzbtzz+i-nvz1+j*2*ext)
      + *(ahalf+nvz1+2*ext-i-1)* *(tmp_vz+i+j*(nvz1+2*ext)));
      *(vz+i+j*(nvz1+2*ext)) = *(vz+i+j*(nvz1+2*ext)) + *(tmp_vz+i+j*(nvz1+2*ext)) + *(pvzbtzz+i-nvz1+j*2*ext);
    }
  }
  // vzbtxz
  #pragma omp parallel for collapse(2)
  for(i=0; i<nvz1+2*ext; i++)
  {
    for(j=2; j<ext; j++)
    {
      *(tmp_vz+i+j*(nvz1+2*ext)) = (*(txz+i+(j-2)*(ntxz1+2*ext))* *(fdc) + *(txz+i+(j-1)*(ntxz1+2*ext))* *(fdc+1)
      + *(txz+i+j*(ntxz1+2*ext))* *(fdc+2) + *(txz+i+(j+1)*(ntxz1+2*ext))* *(fdc+3));
      *(tmp_vz+i+j*(nvz1+2*ext)) = *(rho_vz+i+j*(nvz1+2*ext))* *(tmp_vz+i+j*(nvz1+2*ext))/dx;
      *(pvzbtxz+i+j*(nvz1+2*ext)) = *(bfull+j) * *(pvzbtxz+i+j*(nvz1+2*ext)) + *(afull+j)* *(tmp_vz+i+j*(nvz1+2*ext));
      *(vz+i+j*(nvz1+2*ext)) = *(vz+i+j*(nvz1+2*ext)) + *(pvzbtxz+i+j*(nvz1+2*ext)) + *(tmp_vz+i+j*(nvz1+2*ext));
    }
  }
  #pragma omp parallel for collapse(2)
  for(i=0; i<nvz1+2*ext; i++)
  {
    for(j=ext; j<ext+nvz2; j++)
    {
      *(tmp_vz+i+j*(nvz1+2*ext)) = (*(txz+i+(j-2)*(ntxz1+2*ext))* *(fdc) + *(txz+i+(j-1)*(ntxz1+2*ext))* *(fdc+1)
      + *(txz+i+j*(ntxz1+2*ext))* *(fdc+2) + *(txz+i+(j+1)*(ntxz1+2*ext))* *(fdc+3));
      *(vz+i+j*(ntxz1+2*ext)) = *(vz+i+j*(ntxz1+2*ext)) + *(rho_vz+i+j*(ntxz1+2*ext))* *(tmp_vz+i+j*(nvz1+2*ext))/dx;
    }
  }
  #pragma omp parallel for collapse(2)
  for(i=0; i<nvz1+2*ext; i++)
  {
    for(j=ext+nvz2; j<nvz2+2*ext-2; j++)
    {
      *(tmp_vz+i+j*(nvz1+2*ext)) = (*(txz+i+(j-2)*(ntxz1+2*ext))* *(fdc) + *(txz+i+(j-1)*(ntxz1+2*ext))* *(fdc+1)
      + *(txz+i+j*(ntxz1+2*ext))* *(fdc+2) + *(txz+i+(j+1)*(ntxz1+2*ext))* *(fdc+3));
      *(tmp_vz+i+j*(nvz1+2*ext)) = *(rho_vz+i+j*(nvz1+2*ext))* *(tmp_vz+i+j*(nvz1+2*ext))/dx;
      *(pvzbtxz+i+(j-nvz2)*(nvz1+2*ext)) = (*(bfull+nvz2+2*ext-j-1) * *(pvzbtxz+i+(j-nvz2)*(nvz1+2*ext))
      + *(afull+nvz2+2*ext-j-1)* *(tmp_vz+i+j*(nvz1+2*ext)));
      *(vz+i+j*(nvz1+2*ext)) = *(vz+i+j*(nvz1+2*ext)) + *(pvzbtxz+i+(j-nvz2)*(nvz1+2*ext)) + *(tmp_vz+i+j*(nvz1+2*ext));
    }
  }

  // tzzbvz
  #pragma omp parallel for collapse(2)
  for(j=0; j<ntpp2+2*ext; j++)
  {
    for(i=2; i<ext; i++)
    {
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = (*(vz+i-2+j*(nvz1+2*ext))* *(fdc) + *(vz+i-1+j*(nvz1+2*ext))* *(fdc+1)
       + *(vz+i+j*(nvz1+2*ext))* *(fdc+2) + *(vz+i+1+j*(nvz1+2*ext))* *(fdc+3));
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = *(lamu+i+j*(ntpp1+2*ext))* *(tmp_tpp+i+j*(ntpp1+2*ext))/dz;
      *(ptzzbvz+i+j*2*ext) = *(bfull+i) * *(ptzzbvz+i+j*2*ext) + *(afull+i)* *(tmp_tpp+i+j*(ntpp1+2*ext));
      *(tzz+i+j*(ntpp1+2*ext)) = *(tzz+i+j*(ntpp1+2*ext)) + *(tmp_tpp+i+j*(ntpp1+2*ext)) + *(ptzzbvz+i+j*2*ext);
    }
  }
  #pragma omp parallel for collapse(2)
  for(j=0; j<ntpp2+2*ext; j++)
  {
    for(i=ext; i<ext+ntpp1; i++)
    {
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = (*(vz+i-2+j*(nvz1+2*ext))* *(fdc) + *(vz+i-1+j*(nvz1+2*ext))* *(fdc+1)
      + *(vz+i+j*(nvz1+2*ext))* *(fdc+2) + *(vz+i+1+j*(nvz1+2*ext))* *(fdc+3));
      *(tzz+i+j*(ntpp1+2*ext)) = *(tzz+i+j*(ntpp1+2*ext)) + *(lamu+i+j*(ntpp1+2*ext))* *(tmp_tpp+i+j*(ntpp1+2*ext))/dz;
    }
  }
  #pragma omp parallel for collapse(2)
  for(j=0; j<ntpp2+2*ext; j++)
  {
    for(i=ext+ntpp1; i<ntpp1+2*ext-2; i++)
    {
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = (*(vz+i-2+j*(nvz1+2*ext))* *(fdc) + *(vz+i-1+j*(nvz1+2*ext))* *(fdc+1)
      + *(vz+i+j*(nvz1+2*ext))* *(fdc+2) + *(vz+i+1+j*(nvz1+2*ext))* *(fdc+3));
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = *(lamu+i+j*(ntpp1+2*ext))* *(tmp_tpp+i+j*(ntpp1+2*ext))/dz;
      *(ptzzbvz+i-ntpp1+j*2*ext) = (*(bfull+ntpp1+2*ext-i-1) * *(ptzzbvz+i-ntpp1+j*2*ext) + *(afull+ntpp1+2*ext-i-1)* *(tmp_tpp+i+j*(ntpp1+2*ext)));
      *(tzz+i+j*(ntpp1+2*ext)) = *(tzz+i+j*(ntpp1+2*ext)) + *(tmp_tpp+i+j*(ntpp1+2*ext)) + *(ptzzbvz+i-ntpp1+j*2*ext);
    }
  }

  // tzzbvx
  #pragma omp parallel for collapse(2)
  for(i=0; i<ntpp1+2*ext; i++)
  {
    for(j=2; j<ext; j++)
    {
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = (*(vx+i+(j-2)*(nvx1+2*ext))* *(fdc) + *(vx+i+(j-1)*(nvx1+2*ext))* *(fdc+1)
      + *(vx+i+j*(nvx1+2*ext))* *(fdc+2) + *(vx+i+(j+1)*(nvx1+2*ext))* *(fdc+3));
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = *(lambda+i+j*(ntpp1+2*ext))*dt* *(tmp_tpp+i+j*(ntpp1+2*ext))/dx;
      *(ptzzbvx+i+j*(ntpp1+2*ext)) = *(bfull+j) * *(ptzzbvx+i+j*(ntpp1+2*ext)) + *(afull+j)* *(tmp_tpp+i+j*(ntpp1+2*ext));
      *(tzz+i+j*(ntpp1+2*ext)) = *(tzz+i+j*(ntpp1+2*ext)) + *(ptzzbvx+i+j*(ntpp1+2*ext)) + *(tmp_tpp+i+j*(ntpp1+2*ext));
    }
  }
  #pragma omp parallel for collapse(2)
  for(i=0; i<ntpp1+2*ext; i++)
  {
    for(j=ext; j<ext+ntpp2; j++)
    {
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = (*(vx+i+(j-2)*(nvx1+2*ext))* *(fdc) + *(vx+i+(j-1)*(nvx1+2*ext))* *(fdc+1)
      + *(vx+i+j*(nvx1+2*ext))* *(fdc+2) + *(vx+i+(j+1)*(nvx1+2*ext))* *(fdc+3));
      *(tzz+i+j*(ntpp1+2*ext)) = *(tzz+i+j*(ntpp1+2*ext)) + *(lambda+i+j*(ntpp1+2*ext))*dt* *(tmp_tpp+i+j*(ntpp1+2*ext))/dx;
    }
  }
  #pragma omp parallel for collapse(2)
  for(i=0; i<ntpp1+2*ext; i++)
  {
    for(j=ext+ntpp2; j<ntpp2+2*ext-2; j++)
    {
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = (*(vx+i+(j-2)*(nvx1+2*ext))* *(fdc) + *(vx+i+(j-1)*(nvx1+2*ext))* *(fdc+1)
      + *(vx+i+j*(nvx1+2*ext))* *(fdc+2) + *(vx+i+(j+1)*(nvx1+2*ext))* *(fdc+3));
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = *(lambda+i+j*(ntpp1+2*ext))*dt* *(tmp_tpp+i+j*(ntpp1+2*ext))/dx;
      *(ptzzbvx+i+(j-ntpp2)*(ntpp1+2*ext)) = (*(bfull+ntpp2+2*ext-j-1) * *(ptzzbvx+i+(j-ntpp2)*(ntpp1+2*ext))
      + *(afull+ntpp2+2*ext-j-1)* *(tmp_tpp+i+j*(ntpp1+2*ext)));
      *(tzz+i+j*(ntpp1+2*ext)) = *(tzz+i+j*(ntpp1+2*ext)) + *(ptzzbvx+i+(j-ntpp2)*(ntpp1+2*ext)) + *(tmp_tpp+i+j*(ntpp1+2*ext));
    }
  }

  // txxbvz
  #pragma omp parallel for collapse(2)
  for(j=0; j<ntpp2+2*ext; j++)
  {
    for(i=2; i<ext; i++)
    {
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = (*(vz+i-2+j*(nvz1+2*ext))* *(fdc) + *(vz+i-1+j*(nvz1+2*ext))* *(fdc+1)
       + *(vz+i+j*(nvz1+2*ext))* *(fdc+2) + *(vz+i+1+j*(nvz1+2*ext))* *(fdc+3));
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = *(lambda+i+j*(ntpp1+2*ext))*dt* *(tmp_tpp+i+j*(ntpp1+2*ext))/dz;
      *(ptxxbvz+i+j*2*ext) = *(bfull+i) * *(ptxxbvz+i+j*2*ext) + *(afull+i)* *(tmp_tpp+i+j*(ntpp1+2*ext));
      *(txx+i+j*(ntpp1+2*ext)) = *(txx+i+j*(ntpp1+2*ext)) + *(tmp_tpp+i+j*(ntpp1+2*ext)) + *(ptxxbvz+i+j*2*ext);
    }
  }
  #pragma omp parallel for collapse(2)
  for(j=0; j<ntpp2+2*ext; j++)
  {
    for(i=ext; i<ext+ntpp1; i++)
    {
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = (*(vz+i-2+j*(nvz1+2*ext))* *(fdc) + *(vz+i-1+j*(nvz1+2*ext))* *(fdc+1)
      + *(vz+i+j*(nvz1+2*ext))* *(fdc+2) + *(vz+i+1+j*(nvz1+2*ext))* *(fdc+3));
      *(txx+i+j*(ntpp1+2*ext)) = *(txx+i+j*(ntpp1+2*ext)) + *(lambda+i+j*(ntpp1+2*ext))*dt* *(tmp_tpp+i+j*(ntpp1+2*ext))/dz;
    }
  }
  #pragma omp parallel for collapse(2)
  for(j=0; j<ntpp2+2*ext; j++)
  {
    for(i=ext+ntpp1; i<ntpp1+2*ext-2; i++)
    {
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = (*(vz+i-2+j*(nvz1+2*ext))* *(fdc) + *(vz+i-1+j*(nvz1+2*ext))* *(fdc+1)
      + *(vz+i+j*(nvz1+2*ext))* *(fdc+2) + *(vz+i+1+j*(nvz1+2*ext))* *(fdc+3));
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = *(lambda+i+j*(ntpp1+2*ext))*dt* *(tmp_tpp+i+j*(ntpp1+2*ext))/dz;
      *(ptxxbvz+i-ntpp1+j*2*ext) = (*(bfull+ntpp1+2*ext-i-1) * *(ptxxbvz+i-ntpp1+j*2*ext) + *(afull+ntpp1+2*ext-i-1)* *(tmp_tpp+i+j*(ntpp1+2*ext)));
      *(txx+i+j*(ntpp1+2*ext)) = *(txx+i+j*(ntpp1+2*ext)) + *(tmp_tpp+i+j*(ntpp1+2*ext)) + *(ptxxbvz+i-ntpp1+j*2*ext);
    }
  }

  // txxbvx
  #pragma omp parallel for collapse(2)
  for(i=0; i<ntpp1+2*ext; i++)
  {
    for(j=2; j<ext; j++)
    {
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = (*(vx+i+(j-2)*(nvx1+2*ext))* *(fdc) + *(vx+i+(j-1)*(nvx1+2*ext))* *(fdc+1)
      + *(vx+i+j*(nvx1+2*ext))* *(fdc+2) + *(vx+i+(j+1)*(nvx1+2*ext))* *(fdc+3));
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = *(lamu+i+j*(ntpp1+2*ext))* *(tmp_tpp+i+j*(ntpp1+2*ext))/dx;
      *(ptxxbvx+i+j*(ntpp1+2*ext)) = *(bfull+j) * *(ptxxbvx+i+j*(ntpp1+2*ext)) + *(afull+j)* *(tmp_tpp+i+j*(ntpp1+2*ext));
      *(txx+i+j*(ntpp1+2*ext)) = *(txx+i+j*(ntpp1+2*ext)) + *(tmp_tpp+i+j*(ntpp1+2*ext)) + *(ptxxbvx+i+j*(ntpp1+2*ext));
    }
  }
  #pragma omp parallel for collapse(2)
  for(i=0; i<ntpp1+2*ext; i++)
  {
    for(j=ext; j<ext+ntpp2; j++)
    {
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = (*(vx+i+(j-2)*(nvx1+2*ext))* *(fdc) + *(vx+i+(j-1)*(nvx1+2*ext))* *(fdc+1)
      + *(vx+i+j*(nvx1+2*ext))* *(fdc+2) + *(vx+i+(j+1)*(nvx1+2*ext))* *(fdc+3));
      *(txx+i+j*(ntpp1+2*ext)) = *(txx+i+j*(ntpp1+2*ext)) + *(lamu+i+j*(ntpp1+2*ext))* *(tmp_tpp+i+j*(ntpp1+2*ext))/dx;
    }
  }
  #pragma omp parallel for collapse(2)
  for(i=0; i<ntpp1+2*ext; i++)
  {
    for(j=ext+ntpp2; j<ntpp2+2*ext-2; j++)
    {
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = (*(vx+i+(j-2)*(nvx1+2*ext))* *(fdc) + *(vx+i+(j-1)*(nvx1+2*ext))* *(fdc+1)
      + *(vx+i+j*(nvx1+2*ext))* *(fdc+2) + *(vx+i+(j+1)*(nvx1+2*ext))* *(fdc+3));
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = *(lamu+i+j*(ntpp1+2*ext))* *(tmp_tpp+i+j*(ntpp1+2*ext))/dx;
      *(ptxxbvx+i+(j-ntpp2)*(ntpp1+2*ext)) = (*(bfull+ntpp2+2*ext-j-1) * *(ptxxbvx+i+(j-ntpp2)*(ntpp1+2*ext))
      + *(afull+ntpp2+2*ext-j-1)* *(tmp_tpp+i+j*(ntpp1+2*ext)));
      *(txx+i+j*(ntpp1+2*ext)) = *(txx+i+j*(ntpp1+2*ext)) + *(tmp_tpp+i+j*(ntpp1+2*ext)) + *(ptxxbvx+i+(j-ntpp2)*(ntpp1+2*ext));
    }
  }

  // txzbvx
  #pragma omp parallel for collapse(2)
  for(j=0; j<ntxz2+2*ext; j++)
  {
    for(i=1; i<ext; i++)
    {
      *(tmp_txz+i+j*(ntxz1+2*ext)) = (*(vx+i-1+j*(nvx1+2*ext))* *(fdc) + *(vx+i+j*(nvx1+2*ext))* *(fdc+1)
      + *(vx+i+1+j*(nvx1+2*ext))* *(fdc+2) + *(vx+i+2+j*(nvx1+2*ext))* *(fdc+3));
      *(tmp_txz+i+j*(ntxz1+2*ext)) = *(mu_txz+i+j*(ntxz1+2*ext))* *(tmp_txz+i+j*(ntxz1+2*ext))/dz;
      *(ptxzbvx+i+j*2*ext) = *(bhalf+i) * *(ptxzbvx+i+j*2*ext) + *(ahalf+i)* *(tmp_txz+i+j*(ntxz1+2*ext));
      *(txz+i+j*(ntxz1+2*ext)) = *(txz+i+j*(ntxz1+2*ext)) + *(tmp_txz+i+j*(ntxz1+2*ext)) + *(ptxzbvx+i+j*2*ext);
    }
  }
  #pragma omp parallel for collapse(2)
  for(j=0; j<ntxz2+2*ext; j++)
  {
    for(i=ext; i<ext+ntxz1; i++)
    {
      *(tmp_txz+i+j*(ntxz1+2*ext)) = (*(vx+i-1+j*(nvx1+2*ext))* *(fdc) + *(vx+i+j*(nvx1+2*ext))* *(fdc+1)
      + *(vx+i+1+j*(nvx1+2*ext))* *(fdc+2) + *(vx+i+2+j*(nvx1+2*ext))* *(fdc+3));
      *(txz+i+j*(ntxz1+2*ext)) = *(txz+i+j*(ntxz1+2*ext)) + *(mu_txz+i+j*(ntxz1+2*ext))* *(tmp_txz+i+j*(ntxz1+2*ext))/dz;
    }
  }
  #pragma omp parallel for collapse(2)
  for(j=0; j<ntxz2+2*ext; j++)
  {
    for(i=ext+ntxz1; i<ntxz1+2*ext-1; i++)
    {
      *(tmp_txz+i+j*(ntxz1+2*ext)) = (*(vx+i-1+j*(nvx1+2*ext))* *(fdc) + *(vx+i+j*(nvx1+2*ext))* *(fdc+1)
      + *(vx+i+1+j*(nvx1+2*ext))* *(fdc+2) + *(vx+i+2+j*(nvx1+2*ext))* *(fdc+3));
      *(tmp_txz+i+j*(ntxz1+2*ext)) = *(mu_txz+i+j*(ntxz1+2*ext))* *(tmp_txz+i+j*(ntxz1+2*ext))/dz;
      *(ptxzbvx+i-ntxz1+j*2*ext) = *(bhalf+ntxz1+2*ext-i-1) * *(ptxzbvx+i-ntxz1+j*2*ext) + *(ahalf+ntxz1+2*ext-i-1)* *(tmp_txz+i+j*(ntxz1+2*ext));
      *(txz+i+j*(ntxz1+2*ext)) = *(txz+i+j*(ntxz1+2*ext)) + *(tmp_txz+i+j*(ntxz1+2*ext)) + *(ptxzbvx+i-ntxz1+j*2*ext);
    }
  }

  // txzbvz
  #pragma omp parallel for collapse(2)
  for(i=0; i<ntxz1+2*ext; i++)
  {
    for(j=1; j<ext; j++)
    {
      *(tmp_txz+i+j*(ntxz1+2*ext)) = (*(vz+i+(j-1)*(nvz1+2*ext))* *(fdc) + *(vz+i+j*(nvz1+2*ext))* *(fdc+1)
      + *(vz+i+(j+1)*(nvz1+2*ext))* *(fdc+2) + *(vz+i+(j+2)*(nvz1+2*ext))* *(fdc+3));
      *(tmp_txz+i+j*(ntxz1+2*ext)) = *(mu_txz+i+j*(ntxz1+2*ext))* *(tmp_txz+i+j*(ntxz1+2*ext))/dx;
      *(ptxzbvz+i+j*(ntxz1+2*ext)) = *(bhalf+j) * *(ptxzbvz+i+j*(ntxz1+2*ext)) + *(ahalf+j)* *(tmp_txz+i+j*(ntxz1+2*ext));
      *(txz+i+j*(ntxz1+2*ext)) = *(txz+i+j*(ntxz1+2*ext)) + *(tmp_txz+i+j*(ntxz1+2*ext)) + *(ptxzbvz+i+j*(ntxz1+2*ext));
    }
  }
  #pragma omp parallel for collapse(2)
  for(i=0; i<ntxz1+2*ext; i++)
  {
    for(j=ext; j<ext+ntxz2; j++)
    {
      *(tmp_txz+i+j*(ntxz1+2*ext)) = (*(vz+i+(j-1)*(nvz1+2*ext))* *(fdc) + *(vz+i+j*(nvz1+2*ext))* *(fdc+1)
      + *(vz+i+(j+1)*(nvz1+2*ext))* *(fdc+2) + *(vz+i+(j+2)*(nvz1+2*ext))* *(fdc+3));
      *(txz+i+j*(ntxz1+2*ext)) = *(txz+i+j*(ntxz1+2*ext)) + *(mu_txz+i+j*(ntxz1+2*ext))* *(tmp_txz+i+j*(ntxz1+2*ext))/dx;
    }
  }
  #pragma omp parallel for collapse(2)
  for(i=0; i<ntxz1+2*ext; i++)
  {
    for(j=ext+ntxz2; j<ntxz2+2*ext-1; j++)
    {
      *(tmp_txz+i+j*(ntxz1+2*ext)) = (*(vz+i+(j-1)*(nvz1+2*ext))* *(fdc) + *(vz+i+j*(nvz1+2*ext))* *(fdc+1)
      + *(vz+i+(j+1)*(nvz1+2*ext))* *(fdc+2) + *(vz+i+(j+2)*(nvz1+2*ext))* *(fdc+3));
      *(tmp_txz+i+j*(ntxz1+2*ext)) = *(mu_txz+i+j*(ntxz1+2*ext))* *(tmp_txz+i+j*(ntxz1+2*ext))/dx;
      *(ptxzbvz+i+(j-ntxz2)*(ntxz1+2*ext)) = *(bhalf+ntxz2+2*ext-j-1) * *(ptxzbvz+i+(j-ntxz2)*(ntxz1+2*ext)) + *(ahalf+ntxz2+2*ext-j-1)* *(tmp_txz+i+j*(ntxz1+2*ext));
      *(txz+i+j*(ntxz1+2*ext)) = *(txz+i+j*(ntxz1+2*ext)) + *(tmp_txz+i+j*(ntxz1+2*ext)) + *(ptxzbvz+i+(j-ntxz2)*(ntxz1+2*ext));
    }
  }

}
