#include<stdio.h>

// vx
void acextrap2domp(double *vx, double *tmp_vx,  int nvx1, int nvx2, double *rho_vx, double *pvxbtxx,
  double *vz, double *tmp_vz,  int nvz1, int nvz2, double *rho_vz, double *pvzbtzz,
  double *txx, double *tzz, double *tmp_tpp, int ntpp1, int ntpp2, double *lambda,
  double *ptxxbvx, double *ptxxbvz, double *ptzzbvx, double *ptzzbvz,
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

  // tzzbvz
  #pragma omp parallel for collapse(2)
  for(j=0; j<ntpp2+2*ext; j++)
  {
    for(i=2; i<ext; i++)
    {
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = (*(vz+i-2+j*(nvz1+2*ext))* *(fdc) + *(vz+i-1+j*(nvz1+2*ext))* *(fdc+1)
       + *(vz+i+j*(nvz1+2*ext))* *(fdc+2) + *(vz+i+1+j*(nvz1+2*ext))* *(fdc+3));
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = *(lambda+i+j*(ntpp1+2*ext))*dt* *(tmp_tpp+i+j*(ntpp1+2*ext))/dz;
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
      *(tzz+i+j*(ntpp1+2*ext)) = *(tzz+i+j*(ntpp1+2*ext)) + *(lambda+i+j*(ntpp1+2*ext))*dt* *(tmp_tpp+i+j*(ntpp1+2*ext))/dz;
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
      *(tmp_tpp+i+j*(ntpp1+2*ext)) = *(lambda+i+j*(ntpp1+2*ext))*dt* *(tmp_tpp+i+j*(ntpp1+2*ext))/dx;
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
      *(txx+i+j*(ntpp1+2*ext)) = *(txx+i+j*(ntpp1+2*ext)) + *(lambda+i+j*(ntpp1+2*ext))*dt* *(tmp_tpp+i+j*(ntpp1+2*ext))/dx;
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
      *(ptxxbvx+i+(j-ntpp2)*(ntpp1+2*ext)) = (*(bfull+ntpp2+2*ext-j-1) * *(ptxxbvx+i+(j-ntpp2)*(ntpp1+2*ext))
      + *(afull+ntpp2+2*ext-j-1)* *(tmp_tpp+i+j*(ntpp1+2*ext)));
      *(txx+i+j*(ntpp1+2*ext)) = *(txx+i+j*(ntpp1+2*ext)) + *(tmp_tpp+i+j*(ntpp1+2*ext)) + *(ptxxbvx+i+(j-ntpp2)*(ntpp1+2*ext));
    }
  }
}
