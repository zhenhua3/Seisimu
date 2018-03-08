#include<stdio.h>
#include"fd_body_bound.c"

extern "C"{
  __global__ void kernel_v(double *dev_vx, double *dev_rho,
  double *dev_txx, double *dev_txy, double *dev_txz,
  double dx, double dy, double dz,
  int nx_vx, int ny_vx, int nz_vx,
  int nx_tpp, int ny_tpp, int nz_tpp,
  int nx_txy, int ny_txy, int nz_txy,
  int nx_txz, int ny_txz, int nz_txz,
  double *dev_fdc, int chunk)
  {
    int tX = threadIdx.x + blockIdx.x*blockDim.x;
    int tY = threadIdx.y + blockIdx.y*blockDim.y;
    int tZ = threadIdx.z + blockIdx.z*blockDim.z;

    // vx
    int tid_vx = tZ + tX*nz_vx + tY*nx_vx*nz_vx;

    int tid_txx_0 = tZ + (tX+0)*nz_tpp + tY*nx_txx*nz_tpp;
    int tid_txx_1 = tZ + (tX+1)*nz_tpp + tY*nx_tpp*nz_tpp;
    int tid_txx_2 = tZ + (tX+2)*nz_tpp + tY*nx_tpp*nz_tpp;
    int tid_txx_3 = tZ + (tX+3)*nz_tpp + tY*nx_tpp*nz_tpp;

    int tid_txy_0 = tZ + (tX+2)*nz_txy + (tY+0)*nx_txy*nz_txy;
    int tid_txy_1 = tZ + (tX+2)*nz_txy + (tY+1)*nx_txy*nz_txy;
    int tid_txy_2 = tZ + (tX+2)*nz_txy + (tY+2)*nx_txy*nz_txy;
    int tid_txy_3 = tZ + (tX+2)*nz_txy + (tY+3)*nx_txy*nz_txy;

    int tid_txz_0 = (tZ+0) + (tX+2)*nz_txz + tY*nx_txz*nz_txz;
    int tid_txz_1 = (tZ+1) + (tX+2)*nz_txz + tY*nx_txz*nz_txz;
    int tid_txz_2 = (tZ+2) + (tX+2)*nz_txz + tY*nx_txz*nz_txz;
    int tid_txz_3 = (tZ+3) + (tX+2)*nz_txz + tY*nx_txz*nz_txz;

    int tid_rho_vx_0 = tZ + (tX+0)*nz_tpp + tY*nx_txx*nz_tpp;
    int tid_rho_vx_1 = tZ + (tX+1)*nz_tpp + tY*nx_txx*nz_tpp;

    double tmp_vx;
    if(tX < nx_vx && tY < chunk && tZ < nz_vx){
      tmp_vx = (*(dev_txx+tid_txx_0)* *(dev_fdc)
      + *(dev_txx+tid_txx_1)* *(dev_fdc+1)
      + *(dev_txx+tid_txx_2)* *(dev_fdc+2)
      + *(dev_txx+tid_txx_3)* *(dev_fdc+3)) / dx;

      tmp_vx = tmp_vx + (*(dev_txy+tid_txy_0)* *(dev_fdc)
      + *(dev_txy+tid_txy_1)* *(dev_fdc+1)
      + *(dev_txy+tid_txy_2)* *(dev_fdc+2)
      + *(dev_txy+tid_txy_3)* *(dev_fdc+3)) / dy;

      tmp_vx = tmp_vx + (*(dev_txz+tid_txz_0)* *(dev_fdc)
      + *(dev_txz+tid_txz_1)* *(dev_fdc+1)
      + *(dev_txz+tid_txz_2)* *(dev_fdc+2)
      + *(dev_txz+tid_txz_3)* *(dev_fdc+3)) / dz;

      tmp_rho_vx = 2/(*(dev_rho + tid_rho_vx_0)
      + *(dev_rho + tid_rho_vx_1));

      *(dev_vx+tid_vx) = *(dev_vx+tid_vx) + tmp_rho_vx * tmp_vx;
    }

    // vy
    int tid_vy = tZ + tX*nz_vy + tY*nx_vy*nz_vy;
    int tid_txy_0 = tZ + (tX+0)*nz_txy + (tY+2)*nx_txy*nz_txy;
    int tid_txy_1 = tZ + (tX+1)*nz_txy + (tY+2)*nx_txy*nz_txy;
    int tid_txy_2 = tZ + (tX+2)*nz_txy + (tY+2)*nx_txy*nz_txy;
    int tid_txy_3 = tZ + (tX+3)*nz_txy + (tY+2)*nx_txy*nz_txy;

    int tid_tyy_0 = tZ + tX*nz_tyy + (tY+0)*nx_tyy*nz_tyy;
    int tid_tyy_1 = tZ + tX*nz_tyy + (tY+1)*nx_tyy*nz_tyy;
    int tid_tyy_2 = tZ + tX*nz_tyy + (tY+2)*nx_tyy*nz_tyy;
    int tid_tyy_3 = tZ + tX*nz_tyy + (tY+3)*nx_tyy*nz_tyy;

    int tid_tyz_0 = (tZ+0) + tX*nz_tyz + (tY+2)*nx_tyz*nz_tyz;
    int tid_tyz_1 = (tZ+1) + tX*nz_tyz + (tY+2)*nx_tyz*nz_tyz;
    int tid_tyz_2 = (tZ+2) + tX*nz_tyz + (tY+2)*nx_tyz*nz_tyz;
    int tid_tyz_3 = (tZ+3) + tX*nz_tyz + (tY+2)*nx_tyz*nz_tyz;

    int tid_rho_vy_0 = tZ + tX*nz_tyy + (tY+0)*nx_tyy*nz_tyy;
    int tid_rho_vy_1 = tZ + tX*nz_tyy + (tY+1)*nx_tyy*nz_tyy;

    double tmp_vy;
    if(tX < nx_vy && tY < chunk && tZ < nz_vy){
      tmp_vy = (*(dev_txy+tid_txy_0)* *(dev_fdc)
      + *(dev_txy+tid_txy_1)* *(dev_fdc+1)
      + *(dev_txy+tid_txy_2)* *(dev_fdc+2)
      + *(dev_txy+tid_txy_3)* *(dev_fdc+3)) / dx;

      tmp_vy = tmp_vy + (*(dev_tyy+tid_tyy_0)* *(dev_fdc)
      + *(dev_tyy+tid_tyy_1)* *(dev_fdc+1)
      + *(dev_tyy+tid_tyy_2)* *(dev_fdc+2)
      + *(dev_tyy+tid_tyy_3)* *(dev_fdc+3)) / dy;

      tmp_vy = tmp_vy + (*(dev_tyz+tid_tyz_0)* *(dev_fdc)
      + *(dev_tyz+tid_tyz_1)* *(dev_fdc+1)
      + *(dev_tyz+tid_tyz_2)* *(dev_fdc+2)
      + *(dev_tyz+tid_tyz_3)* *(dev_fdc+3)) / dz;

      tmp_rho_vy = 2/(*(dev_rho + tid_rho_vy_0)
      + *(dev_rho + tid_rho_vy_1));

      *(dev_vy+tid_vy) = *(dev_vy+tid_vy) + tmp_rho_vy * tmp_vy;
    }

    // vz
    int tid_vz = tZ + tX*nz_vz + tY*nx_vz*nz_vz;
    int tid_txz_0 = tZ+2 + (tX+0)*nz_txz + tY*nx_txz*nz_txz;
    int tid_txz_1 = tZ+2 + (tX+1)*nz_txz + tY*nx_txz*nz_txz;
    int tid_txz_2 = tZ+2 + (tX+2)*nz_txz + tY*nx_txz*nz_txz;
    int tid_txz_3 = tZ+2 + (tX+3)*nz_txz + tY*nx_txz*nz_txz;

    int tid_tyz_0 = tZ+2 + tX*nz_tyz + (tY+0)*nx_tyz*nz_tyz;
    int tid_tyz_1 = tZ+2 + tX*nz_tyz + (tY+1)*nx_tyz*nz_tyz;
    int tid_tyz_2 = tZ+2 + tX*nz_tyz + (tY+2)*nx_tyz*nz_tyz;
    int tid_tyz_3 = tZ+2 + tX*nz_tyz + (tY+3)*nx_tyz*nz_tyz;

    int tid_tzz_0 = (tZ+0) + tX*nz_tzz + tY*nx_tzz*nz_tzz;
    int tid_tzz_1 = (tZ+1) + tX*nz_tzz + tY*nx_tzz*nz_tzz;
    int tid_tzz_2 = (tZ+2) + tX*nz_tzz + tY*nx_tzz*nz_tzz;
    int tid_tzz_3 = (tZ+3) + tX*nz_tzz + tY*nx_tzz*nz_tzz;

    int tid_rho_vz_0 = (tZ+0) + tX*nz_tzz + tY*nx_tzz*nz_tzz;
    int tid_rho_vz_1 = (tZ+1) + tX*nz_tzz + tY*nx_tzz*nz_tzz;

    double tmp_vz;
    if(tX < nx_vz && tY < chunk && tZ < nz_vz){
      tmp_vz = (*(dev_txz+tid_txz_0)* *(dev_fdc)
      + *(dev_txz+tid_txz_1)* *(dev_fdc+1)
      + *(dev_txz+tid_txz_2)* *(dev_fdc+2)
      + *(dev_txz+tid_txz_3)* *(dev_fdc+3)) / dx;

      tmp_vz = tmp_vz + (*(dev_tyz+tid_tyz_0)* *(dev_fdc)
      + *(dev_tyz+tid_tyz_1)* *(dev_fdc+1)
      + *(dev_tyz+tid_tyz_2)* *(dev_fdc+2)
      + *(dev_tyz+tid_tyz_3)* *(dev_fdc+3)) / dy;

      tmp_vz = tmp_vz + (*(dev_tzz+tid_tzz_0)* *(dev_fdc)
      + *(dev_tzz+tid_tzz_1)* *(dev_fdc+1)
      + *(dev_tzz+tid_tzz_2)* *(dev_fdc+2)
      + *(dev_tzz+tid_tzz_3)* *(dev_fdc+3)) / dz;

      tmp_rho_vz = 2/(*(dev_rho + tid_rho_vz_0)
      + *(dev_rho + tid_rho_vz_1));

      *(dev_vz+tid_vz) = *(dev_vz+tid_vz) + tmp_rho_vz * tmp_vz;
    }
  }
  //double *fdc/
  //double *intvl: dz, dx, dt
  //double *modelsize : BDnDZ, BDnHX, nDZ, nHX
  el3d_cump(double *vx, double *tmp_vx, int nx_vx, int ny_vx, int nz_vx, int BD_nx_vx, int BD_ny_vx, int BD_vz_vx,
  double *pvxbtxx, double *pvxbtxy, double *pvxbtxz,
  double *vy, double *tmp_vy, int nx_vy, int ny_vy, int nz_vy, int BD_nx_vy, int BD_ny_vy, int BD_vz_vy,
  double *pvybtxy, double *pvybtyy, double *bvybtyz,
  double *vz, double *tmp_vz, int nx_vz, int ny_vz, int nz_vz, int BD_nx_vz, int BD_ny_vz, int BD_vz_vz,
  double *pvzbtxz, double *pvzbtyz, double *bvzbtzz,
  double *txx, double *tmp_txx,
  double *ptxxbvx, double *ptxxbvy, double *ptxxbvz,
  double *tyy, double *tmp_tyy,
  double *ptyybvx, double *ptyybvy, double *ptyybvz,
  double *tzz, double *tmp_tzz,
  double *ptzzbvx, double *ptzzbvy, double *ptzzbvz,
  int nx_tpp, int ny_tpp, int nz_tpp, int BD_nx_tpp, int BD_ny_tpp, int BD_vz_tpp,
  double *txy, double *tmp_txy, int nx_txy, int ny_txy, int nz_txy, int BD_nx_txy, int BD_ny_txy, int BD_nz_txy,
  double *ptxybvx, double *ptxybvy,
  double *tyz, double *tmp_tyz, int nx_tyz, int ny_tyz, int nz_tyz, int BD_nx_tyz, int BD_ny_tyz, int BD_nz_tyz,
  double *ptyzbvy, double *btyzbvz,
  double *txz, double *tmp_txz, int nx_txz, int ny_txz, int nz_txz, int BD_nx_txz, int BD_ny_txz, int BD_nz_txz,
  double *ptxzbvx, double *btxzbvz,
  double *rho, double dt, double dx, double dy, double dz,
  double *bhalf, double *ahalf, double *bfull, double *afull,
  long int *chk_group_full, long int *chk_group_half, long int Num_group, long int Max_group_dim,
  long int *chk_stream_full, long int *chk_stream_half, long int Max_num_stream, long int Max_stream_dim,
  long int *threadim, long int *blockdim)
  {
    //************************************************//
    //**************** GPU setting *******************//
    //************************************************//
    int gn,sn; // gn : group number; sn : stream number
    cudaStream_t stream[Num_stream]; // create streams for GPU
    for(sn=0 ; sn < Num_stream ; sn++){
             cudaStreamCreate(&stream[sn]); // create concurrent streams
    }
    cudaMalloc((void**)&dev_fdc, 4*sizeof(double)); // copy fdc to device
    cudaMemcpy(dev_fdc,fdc,4*sizeof(double),cudaMemcpyHostToDevice);
    // number of threads and blocks used
    dim3 threads(*threadim,*(threadim+1),*(threadim+2));
    dim3 blocks(*blockdim,*(blockdim+1),*(blockdim+2));
    // pinned host memory for faster transfer between host and device mem
    double *host_vx, *host_vy, *host_vz;
    double *host_txx; *host_tyy, *host_tzz;
    double *host_txy, *host_tyz, *host_txz;
    double *host_rho, *host_lambda, *host_mu;

    cudaHostAlloc((void**)&host_vx, (Max_group_dim+3)*(nx_tpp+3)*(nz_tpp+3)*sizeof(double), cudaHostAllocDefault);
    cudaHostAlloc((void**)&host_vy, (Max_group_dim+3)*(nx_tpp+3)*(nz_tpp+3)*sizeof(double), cudaHostAllocDefault);
    cudaHostAlloc((void**)&host_vz, (Max_group_dim+3)*(nx_tpp+3)*(nz_tpp+3)*sizeof(double), cudaHostAllocDefault);
    cudaHostAlloc((void**)&host_txx, Max_group_dim*(nx_tpp+3)*nz_tpp*sizeof(double), cudaHostAllocDefault);
    cudaHostAlloc((void**)&host_tyy, (Max_group_dim+3)*nx_tpp*nz_tpp*sizeof(double), cudaHostAllocDefault);
    cudaHostAlloc((void**)&host_tzz, Max_group_dim*nx_tpp*(nz_tpp+3)*sizeof(double), cudaHostAllocDefault);
    cudaHostAlloc((void**)&host_txy, (Max_group_dim+3)*(nx_tpp+3)*nz_txy*sizeof(double), cudaHostAllocDefault);
    cudaHostAlloc((void**)&host_tyz, (Max_group_dim+3)*nx_tpp*(nz_tpp+3)*sizeof(double), cudaHostAllocDefault);
    cudaHostAlloc((void**)&host_txz, Max_group_dim*(nx_tpp+3)*(nz_tpp+3)*sizeof(double), cudaHostAllocDefault);
    cudaHostAlloc((void**)&host_rho, (Max_group_dim+1)*nx_tpp*nz_tpp*sizeof(double), cudaHostAllocDefault);
    cudaHostAlloc((void**)&host_lambda, Max_group_dim*nx_tpp*nz_tpp*sizeof(double), cudaHostAllocDefault);
    cudaHostAlloc((void**)&host_mu, (Max_group_dim+1)*nx_tpp*nz_tpp*sizeof(double), cudaHostAllocDefault);

    // device memory for each stream
    double *dev_vx[Max_num_stream], *dev_vy[Max_num_stream], *dev_vz[Max_num_stream];
    double *dev_txx[Max_num_stream], *dev_tyy[Max_num_stream], *dev_tzz[Max_num_stream];
    double *dev_txy[Max_num_stream], *dev_tyz[Max_num_stream], *dev_txz[Max_num_stream];
    double *dev_rho[Max_num_stream], *dev_lambda[Max_num_stream], *dev_mu[Max_num_stream];

    for(i=0 ; i<Max_num_stream ; i++){
      cudaMalloc((void**)&dev_vx[i], Max_stream_dim*(nx_tpp+3)*(nz_tpp+3)*sizeof(double));
      cudaMalloc((void**)&dev_vy, (Max_stream_dim+3)*(nx_tpp+3)*(nz_tpp+3)*sizeof(double), cudaHostAllocDefault);
      cudaMalloc((void**)&dev_vz, (Max_stream_dim+3)*(nx_tpp+3)*(nz_tpp+3)*sizeof(double), cudaHostAllocDefault);
      cudaMalloc((void**)&dev_txx, Max_stream_dim*(nx_tpp+3)*nz_tpp*sizeof(double), cudaHostAllocDefault);
      cudaMalloc((void**)&dev_tyy, (Max_stream_dim+3)*nx_tpp*nz_tpp*sizeof(double), cudaHostAllocDefault);
      cudaMalloc((void**)&dev_tzz, Max_stream_dim*nx_tpp*(nz_tpp+3)*sizeof(double), cudaHostAllocDefault);
      cudaMalloc((void**)&dev_txy, (Max_stream_dim+3)*(nx_tpp+3)*nz_txy*sizeof(double), cudaHostAllocDefault);
      cudaMalloc((void**)&dev_tyz, (Max_stream_dim+3)*nx_tpp*(nz_tpp+3)*sizeof(double), cudaHostAllocDefault);
      cudaMalloc((void**)&dev_txz, Max_stream_dim*(nx_tpp+3)*(nz_tpp+3)*sizeof(double), cudaHostAllocDefault);
      cudaMalloc((void**)&dev_rho, (Max_stream_dim+1)*nx_tpp*nz_tpp*sizeof(double), cudaHostAllocDefault);
      cudaMalloc((void**)&dev_lambda, Max_stream_dim*nx_tpp*nz_tpp*sizeof(double), cudaHostAllocDefault);
      cudaMalloc((void**)&dev_mu, (Max_stream_dim+1)*nx_tpp*nz_tpp*sizeof(double), cudaHostAllocDefault);
    }

    //*************************************************//
    // PML boundary is calculated using openmp since   //
    // it is less efficient to create different        //
    // scenario for boundary and inner body part.      //
    //*************************************************//

    //*************** CPU part ***************//
    int i,j,k;
    #pragma omp parallel
    {
      //******************** vxbtxx ********************//
      #pragma omp for collapse(3) nowait
      for(i=0; i<BD_nz_vx; i++){
        for(j=1; j<ext; j++){
          for(k=0; k<BD_ny_vx; k++){
            bound_x(*vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              *txx, BD_nz_tpp, BD_nx_tpp, j-1,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dx, *BDwf, *BDcoeff_b, *BDcoeff_a, ext, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=0; i<ext; i++){
        for(j=ext; j<ext+nx_vx; j++){
          for(k=0; k<BD_ny_vx; k++){
            body_x(*vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
            *txx, BD_nz_tpp, BD_nx_tpp, j-1,
            0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            dx, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext+nz_vx; i<BD_nz_vx; i++){
        for(j=ext; j<ext+nx_vx; j++){
          for(k=0; k<BD_ny_vx; k++){
            body_x(*vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              *txx, BD_nz_tpp, BD_nx_tpp, j-1,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dx, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vx; i++){
        for(j=ext; j<ext+nx_vx; j++){
          for(k=0; k<ext; k++){
            body_x(*vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              *txx, BD_nz_tpp, BD_nx_tpp, j-1,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dx, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vx; i++){
        for(j=ext; j<ext+nx_vx; j++){
          for(k=ext+ny_vx; k<2*ext+ny_vx; k++){
            body_x(*vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              *txx, BD_nz_tpp, BD_nx_tpp, j-1,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dx, *fdc);
          }
        }
      }

      //******************* vxbtxy ********************//
      #pragma omp for collapse(3) nowait
      for(i=0; i<BD_nz_vx; i++){
        for(j=0; j<BD_nx_vx; j++){
          for(k=2; k<ext; k++){
            bound_y(*vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
                *txy, BD_nz_txy, BD_nx_txy, k-2,
                0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
                dy, *BDwf, *BDcoeff_b, *BDcoeff_a, ext, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=0; i<ext; i++){
        for(j=0; j<BD_nx_vx; j++){
          for(k=ext; k<ext+ny_vx; k++){
            body_y(*vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              *txy, BD_nz_txy, BD_nx_txy, k-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dy, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext+nz_vx; i<BD_nz_vx; i++){
        for(j=0; j<BD_nx_vx; j++){
          for(k=ext; k<ext+ny_vx; k++){
            body_y(*vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              *txy, BD_nz_txy, BD_nx_txy, k-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dy, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vx; i++){
        for(j=0; j<ext; j++){
          for(k=ext; k<ext+ny_vx; k++){
            body_y(*vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              *txy, BD_nz_txy, BD_nx_txy, k-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dy, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vx; i++){
        for(j=ext+nx_vx; j<BD_nx_vx; j++){
          for(k=ext; k<ext+ny_vx; k++){
            body_y(*vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              *txy, BD_nz_txy, BD_nx_txy, k-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+（j+1）* BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dy, *fdc);
          }
        }
      }

      //******************* vxbtxz ********************//
      #pragma omp for collapse(3) nowait
      for(i=2; i<ext; i++){
        for(j=0; j<BD_nx_vx; j++){
          for(k=0; k<BD_ny_vx; k++){
            bound_z(*vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              *txz, BD_nz_txz, BD_nx_txz, i-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dz, *BDwf, *BDcoeff_b, *BDcoeff_a, ext, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vx; i++){
        for(j=0; j<BD_nx_vx; j++){
          for(k=0; k<ext; k++){
            body_z(*vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              *txz, BD_nz_txz, BD_nx_txz, i-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+(1+j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dz, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vx; i++){
        for(j=0; j<BD_nx_vx; j++){
          for(k=ext+ny_vx; k<BD_ny_vx; k++){
            body_z(*vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              *txz, BD_nz_txz, BD_nx_txz, i-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+(1+j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dz, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vx; i++){
        for(j=0; j<ext; j++){
          for(k=ext; k<ext+ny_vx; k++){
            body_y(*vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
              *txz, BD_nz_txz, BD_nx_txz, j-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+(1+j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dz, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vx; i++){
        for(j=ext+nx_vx; j<BD_nx_vx; j++){
          for(k=ext; k<ext+ny_vx; k++){
            body_y(*vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
            *txz, BD_nz_txz, BD_nx_txz, j-2,
            0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+(1+j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            dz, *fdc);
          }
        }
      }

      //******************** vybtxy ********************//
      #pragma omp for collapse(3) nowait
      for(i=0; i<BD_nz_vy; i++){
        for(j=2; j<ext; j++){
          for(k=0; k<BD_ny_vy; k++){
            bound_x(*vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              *txy, BD_nz_txy, BD_nx_txy, j-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              dx, *BDwf, *BDcoeff_b, *BDcoeff_a, ext, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=0; i<ext; i++){
        for(j=ext; j<ext+nx_vy; j++){
          for(k=0; k<BD_ny_vy; k++){
            body_x(*vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              *txy, BD_nz_txy, BD_nx_txy, j-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              dx, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext+nz_vy; i<BD_nz_vy; i++){
        for(j=ext; j<ext+nx_vy; j++){
          for(k=0; k<BD_ny_vy; k++){
            body_x(*vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              *txy, BD_nz_txy, BD_nx_txy, j-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+j*BD_nz_tpp+(k+1)BD_nz_tpp*BD_nx_tpp)),
              dx, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vy; i++){
        for(j=ext; j<ext+nx_vy; j++){
          for(k=0; k<ext; k++){
            body_x(*vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              *txy, BD_nz_txy, BD_nx_txy, j-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              dx, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vy; i++){
        for(j=ext; j<ext+nx_vy; j++){
          for(k=ext+ny_vy; k<2*ext+ny_vy; k++){
            body_x(*vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              *txy, BD_nz_txy, BD_nx_txy, j-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              dx, *fdc);
          }
        }
      }

      //******************* vybtyy ********************//
      #pragma omp for collapse(3) nowait
      for(i=0; i<BD_nz_vy; i++){
        for(j=0; j<BD_nx_vy; j++){
          for(k=1; k<ext; k++){
            bound_y(*vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              *tyy, BD_nz_tpp, BD_nx_tpp, k-1,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              dy, *BDwf, *BDcoeff_b, *BDcoeff_a, ext, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=0; i<ext; i++){
        for(j=0; j<BD_nx_vy; j++){
          for(k=ext; k<ext+ny_vy; k++){
            body_y(*vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              *tyy, BD_nz_tpp, BD_nx_tpp, k-1,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              dy, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext+nz_vy; i<BD_nz_vy; i++){
        for(j=0; j<BD_nx_vy; j++){
          for(k=ext; k<ext+ny_vy; k++){
            body_y(*vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              *tyy, BD_nz_tpp, BD_nx_tpp, k-1,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              dy, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vy; i++){
        for(j=0; j<ext; j++){
          for(k=ext; k<ext+ny_vy; k++){
            body_y(*vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              *tyy, BD_nz_tpp, BD_nx_tpp, k-1,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              dy, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vy; i++){
        for(j=ext+nx_vy; j<BD_nx_vy; j++){
          for(k=ext; k<ext+ny_vy; k++){
            body_y(*vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              *tyy, BD_nz_tpp, BD_nx_tpp, k-1,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+j*BD_nz_tpp+（k+1）*BD_nz_tpp*BD_nx_tpp)),
              dy, *fdc);
          }
        }
      }

      //******************* vybtyz ********************//
      #pragma omp for collapse(3) nowait
      for(i=2; i<ext; i++){
        for(j=0; j<BD_nx_vy; j++){
          for(k=0; k<BD_ny_vy; k++){
            bound_z(*vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              *tyz, BD_nz_tyz, BD_nx_tyz, i-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              dz, *BDwf, *BDcoeff_b, *BDcoeff_a, ext, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vy; i++){
        for(j=0; j<BD_nx_vy; j++){
          for(k=0; k<ext; k++){
            body_z(*vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              *tyz, BD_nz_tyz, BD_nx_tyz, i-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              dz, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vy; i++){
        for(j=0; j<BD_nx_vy; j++){
          for(k=ext+ny_vy; k<BD_ny_vy; k++){
            body_z(*vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              *tyz, BD_nz_tyz, BD_nx_tyz, i-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              dz, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vy; i++){
        for(j=0; j<ext; j++){
          for(k=ext; k<ext+ny_vy; k++){
            body_y(*vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              *tyz, BD_nz_tyz, BD_nx_tyz, j-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              dz, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vy; i++){
        for(j=ext+nx_vy; j<BD_nx_vy; j++){
          for(k=ext; k<ext+ny_vy; k++){
            body_y(*vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
              *tyz, BD_nz_tyz, BD_nx_tyz, j-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
              dz, *fdc);
          }
        }
      }

      //******************** vzbtxz ********************//
      #pragma omp for collapse(3) nowait
      for(i=0; i<BD_nz_vz; i++){
        for(j=2; j<ext; j++){
          for(k=0; k<BD_ny_vz; k++){
            bound_x(*vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              *txz, BD_nz_txz, BD_nx_txz, j-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dx, *BDwf, *BDcoeff_b, *BDcoeff_a, ext, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=0; i<ext; i++){
        for(j=ext; j<ext+nx_vz; j++){
          for(k=0; k<BD_ny_vz; k++){
            body_x(*vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              *txz, BD_nz_txz, BD_nx_txz, j-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dx, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext+nz_vz; i<BD_nz_vz; i++){
        for(j=ext; j<ext+nx_vz; j++){
          for(k=0; k<BD_ny_vz; k++){
            body_x(*vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              *txz, BD_nz_txz, BD_nx_txz, j-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dx, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vz; i++){
        for(j=ext; j<ext+nx_vz; j++){
          for(k=0; k<ext; k++){
            body_x(*vy, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              *txz, BD_nz_txz, BD_nx_txz, j-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dx, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vz; i++){
        for(j=ext; j<ext+nx_vz; j++){
          for(k=ext+ny_vz; k<2*ext+ny_vz; k++){
            body_x(*vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              *txz, BD_nz_txz, BD_nx_txz, j-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dx, *fdc);
          }
        }
      }

      //******************* vzbtyz ********************//
      #pragma omp for collapse(3) nowait
      for(i=0; i<BD_nz_vz; i++){
        for(j=0; j<BD_nx_vz; j++){
          for(k=2; k<ext; k++){
            bound_y(*vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              *tyz, BD_nz_tyz, BD_nx_tyz, k-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dy, *BDwf, *BDcoeff_b, *BDcoeff_a, ext, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=0; i<ext; i++){
        for(j=0; j<BD_nx_vz; j++){
          for(k=ext; k<ext+ny_vz; k++){
            body_y(*vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              *tyz, BD_nz_tyz, BD_nx_tyz, k-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dy, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext+nz_vz; i<BD_nz_vz; i++){
        for(j=0; j<BD_nx_vz; j++){
          for(k=ext; k<ext+ny_vz; k++){
            body_y(*vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              *tyz, BD_nz_tyz, BD_nx_tyz, k-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dy, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vz; i++){
        for(j=0; j<ext; j++){
          for(k=ext; k<ext+ny_vz; k++){
            body_y(*vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              *tyz, BD_nz_tyz, BD_nx_tyz, k-2,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dy, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vz; i++){
        for(j=ext+nx_vz; j<BD_nx_vz; j++){
          for(k=ext; k<ext+ny_vz; k++){
            body_y(*vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              *tyz, BD_nz_tyz, BD_nx_tyz, k-1,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dy, *fdc);
          }
        }
      }

      //******************* vzbtzz ********************//
      #pragma omp for collapse(3) nowait
      for(i=1; i<ext; i++){
        for(j=0; j<BD_nx_vz; j++){
          for(k=0; k<BD_ny_vz; k++){
            bound_z(*vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              *tzz, BD_nz_tpp, BD_nx_tpp, i-1,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dz, *BDwf, *BDcoeff_b, *BDcoeff_a, ext, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vz; i++){
        for(j=0; j<BD_nx_vz; j++){
          for(k=0; k<ext; k++){
            body_z(*vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              *tzz, BD_nz_tpp, BD_nx_tpp, i-1,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dz, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vz; i++){
        for(j=0; j<BD_nx_vz; j++){
          for(k=ext+ny_vz; k<BD_ny_vz; k++){
            body_z(*vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              *tzz, BD_nz_tzz, BD_nx_tzz, i-1,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dz, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vz; i++){
        for(j=0; j<ext; j++){
          for(k=ext; k<ext+ny_vz; k++){
            body_y(*vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              *tzz, BD_nz_tzz, BD_nx_tzz, j-1,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dz, *fdc);
          }
        }
      }
      #pragma omp for collapse(3) nowait
      for(i=ext; i<ext+nz_vz; i++){
        for(j=ext+nx_vz; j<BD_nx_vz; j++){
          for(k=ext; k<ext+ny_vz; k++){
            body_y(*vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
              *tzz, BD_nz_tpp, BD_nx_tpp, j-1,
              0.5*(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
              dz, *fdc);
          }
        }
      }
      ////////////////// GPU part /////////////////////
      // copy memory to pinned memory with a group
      // group number is larger than 1
      int offset;
      int start_point = 0;
      #pragma omp for nowait
      for(gn=0 ; gn<Num_group; gn++){
        start_point = start_point + *(chk_group_full+gn) + ext;
        offset_full = *(chk_group_full+gn+1);
        offset_half = *(chk_group_half+gn+1);
        // vx, vz, txx, tzz, txz
        for(k=0 ; k<offset_full ; k++){
          // copy vx to pinned memory
          for(j=ext;j<ext+nx_vx;j++){
            memcpy(host_vx+k*nz_vx*nx_vx+(j-ext)*nz_vx, vx+(start_point+k)*(BD_nz_vx*BD_nx_vx)+j*BD_nz_vx+ext, nz_vx*sizeof(double));
          }
          // copy vz to pinned memory
          for(j=ext;j<ext+nx_vz;j++){
            memcpy(host_vz+k*nz_vz*nx_vz+(j-ext)*nz_vx, vz+(start_point+k)*(BD_nz_vz*BD_nx_vz)+j*BD_nz_vz+ext, nz_vz*sizeof(double));
          }
          // copy txx to pinned memory
          for(j=ext-1;j<ext+nx_vx+2;j++){
            memcpy(host_txx+k*nz_tpp*(nx_vx+3)+j*nz_tpp, txx+(start_point+k)*(BD_nz_tpp*BD_nx_tpp)+j*BD_nz_tpp+ext, nz_tpp*sizeof(double));
          }
          // copy tzz to pinned memory
          for(j=ext;j<ext+nx_tpp;j++){
            memcpy(host_tzz+k*(nz_vz+3)*nx_tpp+j*(nz_vz+3), tzz+(start_point+k)*(BD_nz_tpp*BD_nx_tpp)+j*BD_nz_tpp+ext-1, (nz_vz+3)*sizeof(double));
          }
          // copy txz to pinned memory
          for(j=ext-2;j<ext+nx_txz+2,j++){
            memcpy(host_txz+k*(nz_txz+4)*(nx_txz+4), txz+(start_point+k)*(BD_nz_tpp*BD_nx_tpp)+j*BD_nz_tpp+ext-2, (nz_txz+4)*sizeof(double))
            }
          }
        // vy, tyy, txy, tyz
        for(k=0 ; k<offset_half ; k++){
          // copy vy to pinned memory
          for(j=ext;j<ext+nx_vy;j++){
            memcpy(host_vy+k*nz_vy*nx_vy, vy+(start_point+k)*(BD_nz_vy*BD_nx_vy)+j*BD_nz_vy+ext, nz_vy*sizeof(double));
          }
        }
        for(k=0 ; k<offset_half+3; k++){
          // copy tyy to pinned memory
          for(j=ext;j<ext+nx_tpp;j++){
            memcpy(host_tyy+k*nz_tpp*nx_tpp, tyy+(start_point+k-1)*(BD_nz_tpp*BD_nx_tpp)+j*BD_nz_tpp+ext, nz_tpp*sizeof(double));
          }
        }
        for(k=0 ; k<offset_full+3; k++){
          //copy txy to pinned memory
          for(j=ext-2;j<ext+nx_txy+2){
            memcpy(host_txy+k*nz_txy*(nx_txy+4)), txy+(start_point+k-2)*(BD_nz_txy*BD_nx_txy)+j*BD_nz_txy+ext, nz_txy*sizeof(double));
          }
          //copy tyz to pinned memory
          for(j=ext;j<ext+nx_tyz){
            memcpy(host_tyz+k*(nz_tyz+4)*nx_tyz, tyz+(start_point+k-2)*(BD_nz_tyz*BD_nx_tyz)+j*BD_nz_tyz+ext-2, (nz_tyz+4)*sizeof(double));
          }
        }
        for(k=0; k<offset_half+1;k++){
          for(j=ext;j<ext+nx_tpp;j++){
            memcpy(host_rho+k*nz_tpp*nx_tpp,rho+(start_point+k)*(BD_nz_tpp*BD_nx_tpp)+j*BD_nz_tpp+ext, nz_tpp*sizeof(double));
          }
        }

        /// copy host memory data to device memory()
        start_point = 0;
        for(sn=0 ; sn<*(chk_stream_full+gn*(Max_num_stream+2)); sn++){
          start_point = start_point + *(chk_stream_full+gn*(Max_num_stream+2)+sn+1);
          offset = *(chk_stream_full+gn*(Max_num_stream+2)+sn+2);
          cudaMemcpyAsync(dev_vx[sn], host_vx+start_point*nz_vx*nx_vx, offset*nz_vx*nx_vx*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
          cudaMemcpyAsync(dev_vz[sn], host_vz+start_point*nz_vz*nx_vz, offset*nz_vz*nx_vz*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
          cudaMemcpyAsync(dev_txx[sn], host_txx+start_point*nz_tpp*(nx_vx+3), offset*nz_vx*(nx_vx+3)*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
          cudaMemcpyAsync(dev_tzz[sn], host_tzz+start_point*(nz_vz+3)*nx_tpp, offset*(nz_vz+3)*nx_vz*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);)
          cudaMemcpyAsync(dev_txz[sn], host_txz+start_point*(nz_vx+3)*(nx_vz+3), offset*(nz_vx+3)*(nx_vz+3)*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
          cudaMemcpyAsync(dev_txy[sn], host_txy+start_point*nz_txy*(nx_vy+3), (offset+3)*nz_vx*(nx_vx+3)*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
          cudaMemcpyAsync(dev_tyz[sn], host_tyz+start_point*(nz_vy+3)*nx_tyz, (offset+3)*(nz_vz+3)*nx_tyz*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        }
        for(sn=0 ; sn<*(chk_stream_half+gn*(Max_num_stream+2)); sn++){
          start_point = start_point + *(chk_stream_half+gn*(Max_num_stream+2)+sn+1);
          offset = *(chk_stream_half+gn*(Max_num_stream+2)+sn+2);
          cudaMemcpyAsync(dev_vy[sn], host_vy+start_point*nz_vy*nx_vy, offset*nz_vy*nx_vy*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
          cudaMemcpyAsync(dev_tyy[sn], host_tyy+start_point*nz_tpp*nx_tpp, (offset+3)*nz_tpp*nx_tpp*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        }

        for(sn=0 ; sn<*(chk_stream_full+gn*(Max_num_stream+2)); sn++){
          kernel_v<<<blocks,threads,0,stream[sn]>>>(dev_vx[stream_n], dev_rho_vx[stream_n], dev_txx[stream_n], dev_fdc, nvx1, ntpp1, *(sgmt+3), dx);
        }

        for(stream_n=0 ; stream_n<*(sgmt+2); stream_n++){
          cudaMemcpyAsync(host_vx+stream_n**(sgmt+3)*nvx1, dev_vx[stream_n], *(sgmt+3)*nvx1*sizeof(double), cudaMemcpyDeviceToHost, stream[stream_n]);

        }
        cudaMemcpyAsync(host_vx+(*(sgmt+2)-1)**(sgmt+3)*nvx1, dev_vx[*(sgmt+2)-1], *(sgmt+4)*nvx1*sizeof(double), cudaMemcpyDeviceToHost, stream[*(sgmt+2)-1]);


        for(stream_n=0 ; stream_n<*(sgmt+2); stream_n++){
                    cudaStreamSynchronize(stream[stream_n]);
        }

        for(k=0 ; k<*(sgmt+1) ; k++){
                    memcpy(vx+(group_n**(sgmt+1)+ext+k)*(nvx1+2*ext)+ext, host_vx+k*nvx1, nvx1*sizeof(double));
        }
      }
}
                // end group
               group_n = *sgmt-1;
              for(k=0 ; k<*(sgmt+4) ; k++){
                               memcpy(host_vx+k*nvx1,vx+(group_n**(sgmt+1)+ext+k)*(nvx1+2*ext)+ext,nvx1*sizeof(double));
                               memcpy(host_rho_vx+k*nvx1,rho_vx+(group_n**(sgmt+1)+ext+k)*(nvx1+2*ext)+ext,nvx1*sizeof(double));
                               memcpy(host_txx+k*ntpp1,txx+(group_n**(sgmt+1)+ext-1+k)*(ntpp1+2*ext)+ext,ntpp1*sizeof(double));
             }
            for(k=*(sgmt+4) ; k<*(sgmt+4)+3 ; k++){
                             memcpy(host_txx+k*ntpp1,txx+(group_n**(sgmt+1)+ext-1+k)*(ntpp1+2*ext)+ext,ntpp1*sizeof(double));
           }
           /// copy host memory data to device memory
           for(stream_n=0 ; stream_n<(*(sgmt+5)-1) ; stream_n++){
                           cudaMemcpyAsync(dev_vx[stream_n], host_vx+stream_n**(sgmt+6)*nvx1,*(sgmt+6)*nvx1*sizeof(double), cudaMemcpyHostToDevice, stream[stream_n]);
                           cudaMemcpyAsync(dev_rho_vx[stream_n], host_rho_vx+stream_n**(sgmt+6)*nvx1,*(sgmt+6)*nvx1*sizeof(double), cudaMemcpyHostToDevice, stream[stream_n]);
                           cudaMemcpyAsync(dev_txx[stream_n], host_txx+stream_n**(sgmt+6)*ntpp1,(*(sgmt+6)+3)*ntpp1*sizeof(double), cudaMemcpyHostToDevice, stream[stream_n]);
         }
        cudaMemcpyAsync(dev_vx[*(sgmt+5)-1], host_vx+(*(sgmt+5)-1)**(sgmt+6)*nvx1,*(sgmt+7)*nvx1*sizeof(double), cudaMemcpyHostToDevice, stream[*(sgmt+5)-1]);
        cudaMemcpyAsync(dev_rho_vx[*(sgmt+5)-1], host_rho_vx+(*(sgmt+5)-1)**(sgmt+6)*nvx1,*(sgmt+7)*nvx1*sizeof(double), cudaMemcpyHostToDevice, stream[*(sgmt+5)-1]);
        cudaMemcpyAsync(dev_txx[*(sgmt+5)-1], host_txx+(*(sgmt+5)-1)**(sgmt+6)*ntpp1,(*(sgmt+7)+3)*ntpp1*sizeof(double), cudaMemcpyHostToDevice, stream[*(sgmt+5)-1]);

        for(stream_n=0 ; stream_n<*(sgmt+5) ; stream_n++){
                         kernel_vx<<<blocks,threads,0,stream[stream_n]>>>(dev_vx[stream_n], dev_rho_vx[stream_n], dev_txx[stream_n], dev_fdc, nvx1, ntpp1, *(sgmt+6), dx);
        }
        for(stream_n=0 ; stream_n<(*(sgmt+5)-1) ; stream_n++){
        cudaMemcpyAsync(host_vx+stream_n**(sgmt+6)*nvx1,dev_vx[stream_n], *(sgmt+6)*nvx1*sizeof(double), cudaMemcpyDeviceToHost, stream[stream_n]);
        // cudaMemcpyAsync(host_rho_vx+stream_n**(sgmt+7)*nvx1, dev_rho_vx[stream_n], *(sgmt+7)*nvx1*sizeof(double), cudaMemcpyDeviceToHost, stream[stream_n]);
        // cudaMemcpyAsync(host_txx+stream_n**(sgmt+7)*ntpp1, dev_txx[stream_n], (*(sgmt+7)+3)*ntpp1*sizeof(double), cudaMemcpyDeviceToHost, stream[stream_n]);
                 }
          cudaMemcpyAsync(host_vx+(*(sgmt+5)-1)**(sgmt+6)*nvx1, dev_vx[*(sgmt+5)-1], *(sgmt+7)*nvx1*sizeof(double), cudaMemcpyDeviceToHost, stream[*(sgmt+5)-1]);
      // cudaMemcpyAsync(host_rho_vx+(*(sgmt+6)-1)**(sgmt+3)*nvx1, dev_rho_vx[*(sgmt+6)-1], *(sgmt+8)*nvx1*sizeof(double), cudaMemcpyDeviceToHost, stream[*(sgmt+6)-1]);
      // cudaMemcpyAsync(host_txx+(*(sgmt+6)-1)**(sgmt+3)*ntpp1, dev_txx[*(sgmt+6)-1], (*(sgmt+8)+3)*ntpp1*sizeof(double), cudaMemcpyDeviceToHost, stream[*(sgmt+6)-1]);

      for(stream_n=0 ; stream_n<*(sgmt+5); stream_n++){
        cudaStreamSynchronize(stream[stream_n]);
      }

      for(k=0 ; k<*(sgmt+4) ; k++){
        memcpy(vx+(group_n**(sgmt+1)+ext+k)*(nvx1+2*ext)+ext,host_vx+k*nvx1,nvx1*sizeof(double));
      }
    for(i=0 ; i<*(sgmt+2) ; i++){
      cudaFree(dev_vx[i]);
      cudaFree(dev_rho_vx[i]);
      cudaFree(dev_txx[i]);
      cudaStreamDestroy(stream[i]);
    }
    cudaFreeHost(host_vx);
    cudaFreeHost(host_rho_vx);
    cudaFreeHost(host_txx);
  }
}
