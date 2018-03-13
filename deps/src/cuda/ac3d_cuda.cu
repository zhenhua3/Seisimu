
#include<stdio.h>
#include"fd_cuda_lib.c"

extern "C"{
  __global__ void kernel_v(double *dev_vx, int BD_nx_vx, int BD_nz_vx,
    double *dev_vy, int BD_nx_vy, int BD_nz_vy,
    double *dev_vz, int BD_nx_vz, int BD_nz_vz,
    double *dev_tpp, int BD_nx_tpp, int BD_nz_tpp,
    double *dev_rho, double *dev_fdc,
    double dx, double dy, double dz, double dt,
    int chunk_full, int chunk_half)
  {
    int tX = threadIdx.x + blockIdx.x*blockDim.x;
    int tY = threadIdx.y + blockIdx.y*blockDim.y;
    int tZ = threadIdx.z + blockIdx.z*blockDim.z;

    int tid_vx, tid_vy, tid_vz;
    int tid_tpp_0, tid_tpp_1, tid_tpp_2, tid_tpp_3;
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

    //***************** vy ********************//
    tid_vy = tZ + tX*BD_nz_vy + tY*BD_nx_vy*BD_nz_vy;
    tid_rho_0 = tZ + tX*BD_nz_tpp + (tY+0)*BD_nx_tpp*BD_nz_tpp;
    tid_rho_1 = tZ + tX*BD_nz_tpp + (tY+1)*BD_nx_tpp*BD_nz_tpp;

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

    //***************** vz ********************//
    tid_vz = tZ + tX*BD_nz_vz + tY*BD_nx_vz*BD_nz_vz;
    tid_rho_0 = (tZ+0) + tX*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;
    tid_rho_1 = (tZ+1) + tX*BD_nz_tpp + tY*BD_nx_tpp*BD_nz_tpp;

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
    double *dev_tpp, int BD_nx_tpp, int BD_nz_tpp,
    double *dev_lambda, double *dev_fdc,
    double dx, double dy, double dz, double dt,
    int chunk_full, int chunk_half)
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

        *(dev_tpp+tid_tpp) = *(dev_tpp+tid_tpp)
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

        *(dev_tpp+tid_tpp) = *(dev_tpp+tid_tpp)
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

        *(dev_tpp+tid_tpp) = *(dev_tpp+tid_tpp)
        + *(dev_lambda+tid_lambda_tpp) * tmp_vz;
      }
    }
  //double *fdc/
  //double *intvl: dz, dx, dt
  //double *modelsize : BDnDZ, BDnHX, nDZ, nHX
  void ac3d_cump(double *vx, int BD_nx_vx, int BD_ny_vx, int BD_nz_vx,
  double *pvxbtpp,
  double *vy, int BD_nx_vy, int BD_ny_vy, int BD_nz_vy,
  double *pvybtpp,
  double *vz, int BD_nx_vz, int BD_ny_vz, int BD_nz_vz,
  double *pvzbtpp,
  double *tpp, double *ptppbvx, double *ptppbvy, double *ptppbvz,
  int BD_nx_tpp, int BD_ny_tpp, int BD_nz_tpp,
  double *rho, double *lambda, double *fdc,
  double dt, double dx, double dy, double dz, int ext,
  double *bhalf, double *ahalf, double *bfull, double *afull,
  long int *chk_group_full, long int *chk_group_half, long int Num_group, long int Max_group_dim,
  long int *chk_stream_full, long int *chk_stream_half, long int Max_num_stream, long int Max_stream_dim,
  long int *threadim, long int *blockdim)
  {
    int i,j,k,t;
    //************************************************//
    //**************** GPU setting *******************//
    //************************************************//
    int gn,sn; // gn : group number; sn : stream number
    cudaStream_t stream[Max_num_stream]; // create streams for GPU
    double *dev_fdc;
    cudaMalloc((void**)&dev_fdc, 4*sizeof(double)); // copy fdc to device
    cudaMemcpy(dev_fdc,fdc,4*sizeof(double),cudaMemcpyHostToDevice);

    for(sn=0 ; sn < Max_num_stream ; sn++){
      cudaStreamCreate(&stream[sn]); // create concurrent streams
    }

    // number of threads and blocks used
    dim3 threads(*threadim,*(threadim+1),*(threadim+2));
    dim3 blocks(*blockdim,*(blockdim+1),*(blockdim+2));
    // pinned host memory for faster transfer between host and device mem
    double *host_vx, *host_vy, *host_vz;
    double *host_tpp,
    double *host_rho, *host_lambda;
    cudaHostAlloc((void**)&host_vx, (Max_group_dim+3)*BD_nx_vx*BD_nz_vx*sizeof(double), cudaHostAllocDefault);
    cudaHostAlloc((void**)&host_vy, (Max_group_dim+3)*BD_nx_vy*BD_nz_vy*sizeof(double), cudaHostAllocDefault);
    cudaHostAlloc((void**)&host_vz, (Max_group_dim+3)*BD_nx_vz*BD_nz_vz*sizeof(double), cudaHostAllocDefault);
    cudaHostAlloc((void**)&host_tpp, Max_group_dim*BD_nx_tpp*BD_nz_tpp*sizeof(double), cudaHostAllocDefault);
    cudaHostAlloc((void**)&host_rho, (Max_group_dim+1)*BD_nx_tpp*BD_nz_tpp*sizeof(double), cudaHostAllocDefault);
    cudaHostAlloc((void**)&host_lambda, Max_group_dim*BD_nx_tpp*BD_nz_tpp*sizeof(double), cudaHostAllocDefault);

    // device memory for each stream
    double *dev_vx[Max_num_stream], *dev_vy[Max_num_stream], *dev_vz[Max_num_stream];
    double *dev_tpp[Max_num_stream];
    double *dev_rho[Max_num_stream], *dev_lambda[Max_num_stream];

    for(i=0 ; i<Max_num_stream ; i++){
      cudaMalloc((void**)&dev_vx[i], (Max_stream_dim+3)*BD_nx_vx*BD_nz_vx*sizeof(double));
      cudaMalloc((void**)&dev_vy[i], (Max_stream_dim+3)*BD_nx_vy*BD_nz_vy*sizeof(double));
      cudaMalloc((void**)&dev_vz[i], (Max_stream_dim+3)*BD_nx_vz*BD_nz_vz*sizeof(double));
      cudaMalloc((void**)&dev_txx[i], Max_stream_dim*BD_nx_tpp*BD_nz_tpp*sizeof(double));
      cudaMalloc((void**)&dev_rho[i], (Max_stream_dim+1)*BD_nx_tpp*BD_nz_tpp*sizeof(double));
      cudaMalloc((void**)&dev_lambda[i], Max_stream_dim*BD_nx_tpp*BD_nz_tpp*sizeof(double));
    }

    int start_group_full = 0;
    int start_group_half = 0;
    int offset_group_full;
    int offset_group_half;
    int start_stream_full[Max_num_stream];
    int start_stream_half[Max_num_stream];
    int offset_stream_full[Max_num_stream];
    int offset_stream_half[Max_num_stream];
    for(gn=0 ; gn<Num_group; gn++)
    {
      start_stream_full[0] = 0;
      start_stream_half[0] = 0;
      for(sn=0;sn<Max_num_stream;sn++)
      {
        offset_stream_full[sn] = *(chk_stream_full+gn*(Max_num_stream+2)+sn+2);
        offset_stream_half[sn] = *(chk_stream_half+gn*(Max_num_stream+2)+sn+2);
      }
      for(sn=1;sn<Max_num_stream;sn++)
      {
        start_stream_full[sn] = start_stream_full[sn-1] + offset_stream_full[sn-1];
        start_stream_half[sn] = start_stream_half[sn-1] + offset_stream_half[sn-1];
      }
    }

    //************ time iteration ************//
    for(t=0;t<1;t++)
    {
      // *****************************************************************//
      // ********************** GPU particle velocities ******************//
      // *****************************************************************//
      start_group_full = 2;
      start_group_half = 1;
    for(gn=0 ; gn<Num_group; gn++)
    {
      start_group_full = start_group_full + *(chk_group_full+gn);
      start_group_half = start_group_half + *(chk_group_half+gn);
      offset_group_full = *(chk_group_full+gn+1);
      offset_group_half = *(chk_group_half+gn+1);
      // vx, vz, txx, tzz, txz
      memcpy(host_vx, vx+start_group_full*BD_nz_vx*BD_nx_vx, offset_group_full*BD_nz_vx*BD_nx_vx*sizeof(double));
      memcpy(host_vy, vy+start_group_half*BD_nz_vy*BD_nx_vy, offset_group_half*BD_nz_vy*BD_nx_vy*sizeof(double));
      memcpy(host_vz, vz+start_group_full*BD_nz_vz*BD_nx_vz, offset_group_full*BD_nz_vz*BD_nx_vz*sizeof(double));
      memcpy(host_tpp, tpp+start_group_full*BD_nz_tpp*BD_nx_tpp, offset_group_full*BD_nz_tpp*BD_nx_tpp*sizeof(double));
      memcpy(host_rho, rho+start_group_half*BD_nz_tpp*BD_nx_tpp, (offset_group_half+1)*BD_nz_tpp*BD_nx_tpp*sizeof(double));

      // copy host memory data to device memory()
      for(sn=0 ; sn<*(chk_stream_full+gn*(Max_num_stream+2)); sn++)
      {
        cudaMemcpyAsync(dev_vx[sn], host_vx+start_stream_full[sn]*BD_nz_vx*BD_nx_vx, offset_stream_full[sn]*BD_nz_vx*BD_nx_vx*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_vz[sn], host_vz+start_stream_full[sn]*BD_nz_vz*BD_nx_vz, offset_stream_full[sn]*BD_nz_vz*BD_nx_vz*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_txx[sn], host_txx+start_stream_full[sn]*BD_nz_tpp*BD_nx_tpp, offset_stream_full[sn]*BD_nz_tpp*BD_nx_tpp*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_tzz[sn], host_tzz+start_stream_full[sn]*BD_nz_tpp*BD_nx_tpp, offset_stream_full[sn]*BD_nz_tpp*BD_nx_tpp*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_txz[sn], host_txz+start_stream_full[sn]*BD_nz_txz*BD_nx_txz, offset_stream_full[sn]*BD_nz_txz*BD_nx_txz*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_txy[sn], host_txy+start_stream_full[sn]*BD_nz_txy*BD_nx_txy, (offset_stream_full[sn]+3)*BD_nz_txy*BD_nx_txy*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_tyz[sn], host_tyz+start_stream_full[sn]*BD_nz_tyz*BD_nx_tyz, (offset_stream_full[sn]+3)*BD_nz_tyz*BD_nx_tyz*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_vy[sn], host_vy+start_stream_half[sn]*BD_nz_vy*BD_nx_vy, offset_stream_half[sn]*BD_nz_vy*BD_nx_vy*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_tyy[sn], host_tyy+start_stream_half[sn]*BD_nz_tpp*BD_nx_tpp, (offset_stream_half[sn]+3)*BD_nz_tpp*BD_nx_tpp*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_rho[sn], host_rho+start_stream_half[sn]*BD_nz_tpp*BD_nx_tpp, (offset_stream_half[sn]+1)*BD_nz_tpp*BD_nx_tpp*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
      }
      // if(err!= cudaSuccess){printf("%s\n","error");}
      for(sn=0 ; sn<*(chk_stream_full+gn*(Max_num_stream+2)); sn++)
      {
        kernel_v<<<blocks,threads,0,stream[sn]>>>
        (dev_vx[sn], BD_nx_vx, BD_nz_vx,
          dev_vy[sn], BD_nx_vy, BD_nz_vy,
          dev_vz[sn], BD_nx_vz, BD_nz_vz,
          dev_txx[sn], dev_tyy[sn], dev_tzz[sn], BD_nx_tpp, BD_nz_tpp,
          dev_txy[sn], BD_nx_txy, BD_nz_txy,
          dev_tyz[sn], BD_nx_tyz, BD_nz_tyz,
          dev_txz[sn], BD_nx_txz, BD_nz_txz,
          dev_rho[sn], dev_fdc, dx, dy, dz, dt,
          offset_stream_full[sn],
          offset_stream_half[sn]);
      }
      for(sn=0 ; sn<*(chk_stream_full+gn*(Max_num_stream+2)); sn++)
      {
        cudaMemcpyAsync(host_vx+start_stream_full[sn]*BD_nz_vx*BD_nx_vx, dev_vx[sn], offset_stream_full[sn]*BD_nz_vx*BD_nx_vx*sizeof(double), cudaMemcpyDeviceToHost, stream[sn]);
        cudaMemcpyAsync(host_vz+start_stream_full[sn]*BD_nz_vz*BD_nx_vz, dev_vz[sn], offset_stream_full[sn]*BD_nz_vz*BD_nx_vz*sizeof(double), cudaMemcpyDeviceToHost, stream[sn]);
        cudaMemcpyAsync(host_vy+start_stream_half[sn]*BD_nz_vy*BD_nx_vy, dev_vy[sn], offset_stream_half[sn]*BD_nz_vy*BD_nx_vy*sizeof(double), cudaMemcpyDeviceToHost, stream[sn]);
        cudaMemcpyAsync(host_txx+start_stream_full[sn]*BD_nz_tpp*BD_nx_tpp, dev_txx[sn], offset_stream_full[sn]*BD_nz_tpp*BD_nx_tpp*sizeof(double), cudaMemcpyDeviceToHost, stream[sn]);
      }
      cudaDeviceSynchronize();

      memcpy(vx+start_group_full*BD_nz_vx*BD_nx_vx, host_vx, offset_group_full*BD_nz_vx*BD_nx_vx*sizeof(double));
      memcpy(vy+start_group_half*BD_nz_vy*BD_nx_vy, host_vy, offset_group_half*BD_nz_vy*BD_nx_vy*sizeof(double));
      memcpy(vz+start_group_full*BD_nz_vz*BD_nx_vz, host_vz, offset_group_full*BD_nz_vz*BD_nx_vz*sizeof(double));
      memcpy(txx+start_group_full*BD_nz_tpp*BD_nx_tpp, host_txx, offset_group_full*BD_nz_tpp*BD_nx_tpp*sizeof(double));
    }
    // *************************************************//
    // ******* openmp particle velocity boundary *******//
    // *************************************************//
    // vxbtxx
    for(k=0; k<BD_ny_vx; k++)
    {
      for(i=0; i<BD_nz_vx; i++)
      {
        for(j=1; j<ext; j++)
        {
          bound_x(vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
            txx, BD_nz_tpp, BD_nx_tpp, 1,
          2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
          2.0/(*(rho+i+(BD_nx_vx-1-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(BD_nx_vx-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
          dx, dt, pvxbtxx, bhalf, ahalf, ext, fdc);
        }
      }
    }
    // vybtxy
    for(k=0; k<BD_ny_vy; k++)
    {
      for(i=0; i<BD_nz_vy; i++)
      {
        for(j=2; j<ext; j++)
        {
          bound_x(vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
            txy, BD_nz_txy, BD_nx_txy, 2,
          2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
          2.0/(*(rho+i+(BD_nx_vy-1-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(BD_nx_vy-1-j)*BD_nz_tpp+(BD_ny_vy-k)*BD_nz_tpp*BD_nx_tpp)),
          dx, dt, pvybtxy, bfull, afull, ext, fdc);
        }
      }
    }
    // vzbtxz
    for(k=0; k<BD_ny_vz; k++)
    {
      for(i=0; i<BD_nz_vz; i++)
      {
        for(j=2; j<ext; j++)
        {
          bound_x(vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
            txz, BD_nz_txz, BD_nx_txz, 2,
          2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
          2.0/(*(rho+i+(BD_nx_vz-1-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(i+1)+(BD_nx_vz-1-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
          dx,dt, pvzbtxz, bfull, afull, ext, fdc);
        }
      }
    }
  //********************* V_Y *********************//
    // vxbtxy
    for(j=0; j<BD_nx_vx; j++)
    {
      for(i=0; i<BD_nz_vx; i++)
      {
        for(k=2; k<ext; k++)
        {
          bound_y(vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
            txy, BD_nz_txy, BD_nx_txy, 2,
            2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            2.0/(*(rho+i+j*BD_nz_tpp+(BD_ny_vx-1-k)*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(j+1)*BD_nz_tpp+(BD_ny_vx-1-k)*BD_nz_tpp*BD_nx_tpp)),
            dy,dt, pvxbtxy, bfull, afull, ext, fdc);
        }
      }
    }
    // vybtyy
    for(j=0; j<BD_nx_vy; j++)
    {
      for(i=0; i<BD_nz_vy; i++)
      {
        for(k=1; k<ext; k++)
        {
          bound_y(vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
            tyy, BD_nz_tpp, BD_nx_tpp, 1,
            2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
            2.0/(*(rho+i+j*BD_nz_tpp+(BD_ny_vy-1-k)*BD_nz_tpp*BD_nx_tpp)+*(rho+i+j*BD_nz_tpp+(BD_ny_vy-k)*BD_nz_tpp*BD_nx_tpp)),
            dy,dt, pvybtyy, bhalf, ahalf, ext, fdc);
        }
      }
    }
    // vzbtyz
    for(j=0; j<BD_nx_vz; j++)
    {
      for(i=0; i<BD_nz_vz; i++)
      {
        for(k=2; k<ext; k++)
        {
          bound_y(vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
            tyz, BD_nz_tyz, BD_nx_tyz, 2,
            2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            2.0/(*(rho+i+j*BD_nz_tpp+(BD_ny_vz-1-k)*BD_nz_tpp*BD_nx_tpp)+*(rho+(i+1)+j*BD_nz_tpp+(BD_ny_vz-1-k)*BD_nz_tpp*BD_nx_tpp)),
            dy,dt, pvzbtyz, bfull, afull, ext, fdc);
        }
      }
    }

  //********************* V_Z *********************//
    // vxbtxz
    for(j=0; j<BD_nx_vx; j++)
    {
      for(k=0; k<BD_ny_vx; k++)
      {
        for(i=2; i<ext; i++)
        {
          unlimited_bound_z(vx, BD_nz_vx, BD_nx_vx, BD_ny_vx, i, j, k,
            txz, BD_nz_txz, BD_nx_txz, 2,
          2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
          2.0/(*(rho+(BD_nz_vx-1-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(BD_nz_vx-1-i)+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
          dz,dt, pvxbtxz, bfull, afull, ext, fdc);
        }
      }
    }
    // vybtyz
    for(j=0; j<BD_nx_vy; j++)
    {
      for(k=0; k<BD_ny_vy; k++)
      {
        for(i=2; i<ext; i++)
        {
          unlimited_bound_z(vy, BD_nz_vy, BD_nx_vy, BD_ny_vy, i, j, k,
            tyz, BD_nz_tyz, BD_nx_tyz, 2,
            2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
            2.0/(*(rho+(BD_nz_vy-1-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(BD_nz_vy-1-i)+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)),
            dz,dt, pvybtyz, bfull, afull, ext, fdc);
        }
      }
    }
    // vzbtzz
    for(j=0; j<BD_nx_vz; j++)
    {
      for(k=0; k<BD_ny_vz; k++)
      {
        for(i=1;i<ext;i++)
        {
          unlimited_bound_z(vz, BD_nz_vz, BD_nx_vz, BD_ny_vz, i, j, k,
            tzz, BD_nz_tpp, BD_nx_tpp, 1,
            2.0/(*(rho+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            2.0/(*(rho+(BD_nz_vz-1-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)+*(rho+(BD_nz_vz-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)),
            dz,dt, pvzbtzz, bhalf, ahalf, ext, fdc);
        }
      }
    }

    //*****************************************************************//
    //**************************** GPU stress *************************//
    //*****************************************************************//
    start_group_full = 2;
    start_group_half = 1;
    for(gn=0 ; gn<Num_group; gn++)
    {
      start_group_full = start_group_full + *(chk_group_full+gn);
      start_group_half = start_group_half + *(chk_group_half+gn);
      offset_group_full = *(chk_group_full+gn+1);
      offset_group_half = *(chk_group_half+gn+1);
      // vx, vz, txx, tzz, txz
        memcpy(host_txx, txx+start_group_full*BD_nz_tpp*BD_nx_tpp, offset_group_full*BD_nz_tpp*BD_nx_tpp*sizeof(double));
        memcpy(host_tzz, tzz+start_group_full*BD_nz_tpp*BD_nx_tpp, offset_group_full*BD_nz_tpp*BD_nx_tpp*sizeof(double));
        memcpy(host_tyy, tyy+start_group_full*BD_nz_tpp*BD_nx_tpp, offset_group_half*BD_nz_tpp*BD_nx_tpp*sizeof(double));
        memcpy(host_txy, txy+start_group_half*BD_nz_txy*BD_nx_txy, offset_group_half*BD_nz_txy*BD_nx_txy*sizeof(double));
        memcpy(host_tyz, tyz+start_group_half*BD_nz_tyz*BD_nx_tyz, offset_group_half*BD_nz_tyz*BD_nx_tyz*sizeof(double));
        memcpy(host_txz, txz+start_group_full*BD_nz_txz*BD_nx_txz, offset_group_full*BD_nz_txz*BD_nx_txz*sizeof(double));
        memcpy(host_vx, vx+(start_group_half-1)*BD_nz_vx*BD_nx_vx, (offset_group_half+3)*BD_nz_vx*BD_nx_vx*sizeof(double));
        memcpy(host_vy, vy+(start_group_full-2)*BD_nz_vy*BD_nx_vy, (offset_group_full+3)*BD_nz_vy*BD_nx_vy*sizeof(double));
        memcpy(host_vz, vz+(start_group_half-1)*BD_nz_vz*BD_nx_vz, (offset_group_half+3)*BD_nz_vz*BD_nx_vz*sizeof(double));
        memcpy(host_lambda, lambda+start_group_full*BD_nz_tpp*BD_nx_tpp, offset_group_full*BD_nz_tpp*BD_nx_tpp*sizeof(double));
        memcpy(host_mu, mu+start_group_half*BD_nz_tpp*BD_nx_tpp, (offset_group_half+1)*BD_nz_tpp*BD_nx_tpp*sizeof(double));

      // copy host memory data to device memory()
      for(sn=0 ; sn<*(chk_stream_full+gn*(Max_num_stream+2)); sn++)
      {
        cudaMemcpyAsync(dev_txx[sn], host_txx+start_stream_full[sn]*BD_nz_tpp*BD_nx_tpp, offset_stream_full[sn]*BD_nz_tpp*BD_nx_tpp*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_tzz[sn], host_tzz+start_stream_full[sn]*BD_nz_tpp*BD_nx_tpp, offset_stream_full[sn]*BD_nz_tpp*BD_nx_tpp*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_tyy[sn], host_tyy+start_stream_full[sn]*BD_nz_tpp*BD_nx_tpp, offset_stream_full[sn]*BD_nz_tpp*BD_nx_tpp*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_txz[sn], host_txz+start_stream_full[sn]*BD_nz_txz*BD_nx_txz, offset_stream_full[sn]*BD_nz_txz*BD_nx_txz*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_txy[sn], host_txy+start_stream_half[sn]*BD_nz_txy*BD_nx_txy, offset_stream_half[sn]*BD_nz_txy*BD_nx_txy*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_tyz[sn], host_tyz+start_stream_half[sn]*BD_nz_tyz*BD_nx_tyz, offset_stream_half[sn]*BD_nz_tyz*BD_nx_tyz*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_vx[sn], host_vx+start_stream_half[sn]*BD_nz_vx*BD_nx_vx, (offset_stream_half[sn]+3)*BD_nz_vx*BD_nx_vx*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_vz[sn], host_vz+start_stream_half[sn]*BD_nz_vz*BD_nx_vz, (offset_stream_half[sn]+3)*BD_nz_vz*BD_nx_vz*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_vy[sn], host_vy+start_stream_full[sn]*BD_nz_vy*BD_nx_vy, (offset_stream_full[sn]+3)*BD_nz_vy*BD_nx_vy*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_lambda[sn], host_lambda+start_stream_full[sn]*BD_nz_tpp*BD_nx_tpp, offset_stream_full[sn]*BD_nz_tpp*BD_nx_tpp*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
        cudaMemcpyAsync(dev_mu[sn], host_mu+start_stream_half[sn]*BD_nz_tpp*BD_nx_tpp, (offset_stream_half[sn]+1)*BD_nz_tpp*BD_nx_tpp*sizeof(double), cudaMemcpyHostToDevice, stream[sn]);
      }
      for(sn=0 ; sn<*(chk_stream_full+gn*(Max_num_stream+2)); sn++)
      {
        kernel_tau<<<blocks,threads,0,stream[sn]>>>
        (dev_vx[sn], BD_nx_vx, BD_nz_vx,
          dev_vy[sn], BD_nx_vy, BD_nz_vy,
          dev_vz[sn], BD_nx_vz, BD_nz_vz,
          dev_txx[sn], dev_tyy[sn], dev_tzz[sn], BD_nx_tpp, BD_nz_tpp,
          dev_txy[sn], BD_nx_txy, BD_nz_txy,
          dev_tyz[sn], BD_nx_tyz, BD_nz_tyz,
          dev_txz[sn], BD_nx_txz, BD_nz_txz,
          dev_lambda[sn], dev_mu[sn], dev_fdc, dx, dy, dz, dt,
          offset_stream_full[sn],
          offset_stream_half[sn]);
      }
      for(sn=0 ; sn<*(chk_stream_full+gn*(Max_num_stream+2)); sn++)
      {
        cudaMemcpyAsync(host_txx+start_stream_full[sn]*BD_nz_tpp*BD_nx_tpp, dev_txx[sn], offset_stream_full[sn]*BD_nz_tpp*BD_nx_tpp*sizeof(double), cudaMemcpyDeviceToHost, stream[sn]);
        cudaMemcpyAsync(host_tyy+start_stream_full[sn]*BD_nz_tpp*BD_nx_tpp, dev_tyy[sn], offset_stream_full[sn]*BD_nz_tpp*BD_nx_tpp*sizeof(double), cudaMemcpyDeviceToHost, stream[sn]);
        cudaMemcpyAsync(host_tzz+start_stream_full[sn]*BD_nz_tpp*BD_nx_tpp, dev_tzz[sn], offset_stream_full[sn]*BD_nz_tpp*BD_nx_tpp*sizeof(double), cudaMemcpyDeviceToHost, stream[sn]);
        cudaMemcpyAsync(host_txz+start_stream_full[sn]*BD_nz_txz*BD_nx_txz, dev_txz[sn], offset_stream_full[sn]*BD_nz_txz*BD_nx_txz*sizeof(double), cudaMemcpyDeviceToHost, stream[sn]);
        cudaMemcpyAsync(host_txy+start_stream_half[sn]*BD_nz_txy*BD_nx_txy, dev_txy[sn], offset_stream_half[sn]*BD_nz_txy*BD_nx_txy*sizeof(double), cudaMemcpyDeviceToHost, stream[sn]);
        cudaMemcpyAsync(host_tyz+start_stream_half[sn]*BD_nz_tyz*BD_nx_tyz, dev_tyz[sn], offset_stream_half[sn]*BD_nz_tyz*BD_nx_tyz*sizeof(double), cudaMemcpyDeviceToHost, stream[sn]);
      }

      cudaDeviceSynchronize();

      memcpy(txx+start_group_full*BD_nz_tpp*BD_nx_tpp, host_txx, offset_group_full*BD_nz_tpp*BD_nx_tpp*sizeof(double));
      memcpy(tzz+start_group_full*BD_nz_tpp*BD_nx_tpp, host_tzz, offset_group_full*BD_nz_tpp*BD_nx_tpp*sizeof(double));
      memcpy(tyy+start_group_full*BD_nz_tpp*BD_nx_tpp, host_tyy, offset_group_full*BD_nz_tpp*BD_nx_tpp*sizeof(double));
      memcpy(txy+start_group_half*BD_nz_txy*BD_nx_txy, host_txy, offset_group_half*BD_nz_txy*BD_nx_txy*sizeof(double));
      memcpy(tyz+start_group_half*BD_nz_tyz*BD_nx_tyz, host_tyz, offset_group_half*BD_nz_tyz*BD_nx_tyz*sizeof(double));
      memcpy(txz+start_group_full*BD_nz_txz*BD_nx_txz, host_txz, offset_group_full*BD_nz_txz*BD_nx_txz*sizeof(double));
    }

  //*************************************************//
  //******* openmp stress boundary *******//
  //*************************************************//
    // #pragma omp for collapse(3)
      for(k=0; k<BD_ny_tpp; k++)
      {
        for(i=0; i<BD_nz_tpp; i++)
        {
          for(j=2; j<ext; j++)
          {
            bound_tpp_x(txx, tyy, tzz, BD_nz_tpp, BD_nx_tpp, BD_ny_tpp,
            i, j, k, vx, BD_nz_vx, BD_nx_vx, 2, lambda, mu, dx, dt,
            ptxxbvx, ptyybvx, ptzzbvx, bhalf, ahalf, ext, fdc);
          }
        }
      }
      // #pragma omp for collapse(3)
      for(k=0; k<BD_ny_txy; k++)
      {
        for(i=0; i<BD_nz_txy; i++)
        {
          for(j=1; j<ext; j++)
          {
            bound_x(txy, BD_nz_txy, BD_nx_txy, BD_ny_txy, i, j, k,
            vy, BD_nz_vy, BD_nx_vy, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+i+(j+1)*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp) + *(mu+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp))/4,
            (*(mu+i+(BD_nx_txy-1-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(BD_nx_txy-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+i+(BD_nx_txy-j)*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(BD_nx_txy-1-j)*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp))/4,
            dx, dt, ptxybvy, bfull, afull, ext, fdc);
          }
        }
      }
  //     // #pragma omp for collapse(3)
      for(k=0; k<BD_ny_txz; k++)
      {
        for(i=0; i<BD_nz_txz; i++)
        {
          for(j=1; j<ext; j++)
          {
            bound_x(txz, BD_nz_txz, BD_nx_txz, BD_ny_txz, i, j, k,
            vz, BD_nz_vz, BD_nx_vz, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(i+1)+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            (*(mu+i+(BD_nx_txz-1-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(BD_nx_txz-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(i+1)+(BD_nx_txz-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+(i+1)+(BD_nx_txz-1-j)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            dx, dt, ptxzbvz, bfull, afull, ext, fdc);
          }
        }
      }
  //
    //************ T_Y ************//
      // #pragma omp for collapse(3)
      for(i=0; i<BD_nz_tpp; i++)
      {
        for(j=0; j<BD_nx_tpp; j++)
        {
          for(k=2; k<ext; k++)
          {
            bound_tpp_y(txx, tyy, tzz, BD_nz_tpp, BD_nx_tpp, BD_ny_tpp,
            i, j, k, vy, BD_nz_vy, BD_nx_vy, 2, lambda, mu, dy, dt,
            ptxxbvy, ptyybvy, ptzzbvy, bhalf, ahalf, ext, fdc);
          }
        }
      }
  //     // #pragma omp for collapse(3)
      for(i=0; i<BD_nz_txy; i++)
      {
        for(j=0; j<BD_nx_txy; j++)
        {
          for(k=1; k<ext; k++)
          {
            bound_y(txy, BD_nz_txy, BD_nx_txy, BD_ny_txy, i, j, k,
            vx, BD_nz_vx, BD_nx_vx, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+i+(j+1)*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp) + *(mu+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp))/4,
            (*(mu+i+j*BD_nz_tpp+(BD_ny_txy-k-1)*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(j+1)*BD_nz_tpp+(BD_ny_txy-k-1)*BD_nz_tpp*BD_nx_tpp)
            + *(mu+i+(j+1)*BD_nz_tpp+(BD_ny_txy-k)*BD_nz_tpp*BD_nx_tpp) + *(mu+i+j*BD_nz_tpp+(BD_ny_txy-k)*BD_nz_tpp*BD_nx_tpp))/4,
            dy, dt, ptxybvx, bfull, afull, ext, fdc);
          }
        }
      }
      for(i=0; i<BD_nz_tyz; i++)
      {
        for(j=0; j<BD_nx_tyz; j++)
        {
          for(k=1; k<ext; k++)
          {
            bound_y(tyz, BD_nz_tyz, BD_nx_tyz, BD_ny_tyz, i, j, k,
            vz, BD_nz_vz, BD_nx_vz, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(i+1)+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp) + *(mu+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            (*(mu+i+j*BD_nz_tpp+(BD_ny_tyz-1-k)*BD_nz_tpp*BD_nx_tpp) + *(mu+i+j*BD_nz_tpp+(BD_ny_tyz-k)*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(i+1)+j*BD_nz_tpp+(BD_ny_tyz-k)*BD_nz_tpp*BD_nx_tpp) + *(mu+(i+1)+j*BD_nz_tpp+(BD_ny_tyz-1-k)*BD_nz_tpp*BD_nx_tpp))/4,
            dy, dt, ptyzbvz, bfull, afull, ext, fdc);
          }
        }
      }
    //************ T_Z ************//
      // #pragma omp for collapse(3)
      for(j=0; j<BD_nx_tpp; j++)
      {
        for(k=0; k<BD_ny_tpp; k++)
        {
          for(i=2; i<ext; i++)
          {
            unlimited_bound_tpp_z(txx, tyy, tzz, BD_nz_tpp, BD_nx_tpp, BD_ny_tpp,
            i, j, k, vz, BD_nz_vz, BD_nx_vz, 2, lambda, mu, dz, dt,
            ptxxbvz, ptyybvz, ptzzbvz, bhalf, ahalf, ext, fdc);
          }
        }
      }
  //     // #pragma omp for collapse(3)
      for(j=0; j<BD_nx_txz; j++)
      {
        for(k=0; k<BD_ny_txz; k++)
        {
          for(i=1; i<ext; i++)
          {
            unlimited_bound_z(txz, BD_nz_txz, BD_nx_txz, BD_ny_txz, i, j, k,
            vx, BD_nz_vx, BD_nx_vx, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(i+1)+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            (*(mu+(BD_nz_txz-1-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+(BD_nz_txz-1-i)+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(BD_nz_txz-i)+(j+1)*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+(BD_nz_txz-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            dz, dt, ptxzbvx, bfull, afull, ext, fdc);
          }
        }
      }
  //     // #pragma omp for collapse(3)
      for(j=0; j<BD_nx_tyz; j++)
      {
        for(k=0; k<BD_ny_tyz; k++)
        {
          for(i=0; i<ext; i++)
          {
            unlimited_bound_z(tyz, BD_nz_tyz, BD_nx_tyz, BD_ny_tyz, i, j, k,
            vy, BD_nz_vy, BD_nx_vy, 1,
            (*(mu+i+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+i+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(i+1)+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp) + *(mu+(i+1)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            (*(mu+(BD_nz_tyz-1-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp) + *(mu+(BD_nz_tyz-1-i)+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp)
            + *(mu+(BD_nz_tyz-i)+j*BD_nz_tpp+(k+1)*BD_nz_tpp*BD_nx_tpp) + *(mu+(BD_nz_tyz-i)+j*BD_nz_tpp+k*BD_nz_tpp*BD_nx_tpp))/4,
            dz, dt, ptyzbvy, bfull, afull, ext, fdc);
          }
        }
      }
  }
  for(i=0 ; i<Max_num_stream ; i++)
    {
      cudaFree(dev_vx[i]);
      cudaFree(dev_vy[i]);
      cudaFree(dev_vz[i]);
      cudaFree(dev_txx[i]);
      cudaFree(dev_tyy[i]);
      cudaFree(dev_tzz[i]);
      cudaFree(dev_txy[i]);
      cudaFree(dev_tyz[i]);
      cudaFree(dev_txz[i]);
      cudaFree(dev_rho[i]);
      cudaFree(dev_lambda[i]);
      cudaFree(dev_mu[i]);
      cudaStreamDestroy(stream[i]);
    }
    cudaFreeHost(host_vx);
    cudaFreeHost(host_vy);
    cudaFreeHost(host_vz);
    cudaFreeHost(host_txx);
    cudaFreeHost(host_tyy);
    cudaFreeHost(host_tzz);
    cudaFreeHost(host_txy);
    cudaFreeHost(host_tyz);
    cudaFreeHost(host_txz);
    cudaFreeHost(host_rho);
    cudaFreeHost(host_lambda);
    cudaFreeHost(host_mu);
  }
}
