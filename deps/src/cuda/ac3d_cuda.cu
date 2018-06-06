extern "C"{
  #include<stdio.h>
  #include<unistd.h>
  #include<sys/stat.h>
  #include<sys/mman.h>
  #include<fcntl.h>
  #include<string.h>
  #include<stdlib.h>
  #include"lib/lib.h"

  //double *fdc/
  //double *intvl: dz, dx, dt
  //double *modelsize : BDnDZ, BDnHX, nDZ, nHX
  void ac3d_cuda(
    double *vx, int BD_nx_vx, int BD_ny_vx, int BD_nz_vx, double *pvxbtpp,
    double *vy, int BD_nx_vy, int BD_ny_vy, int BD_nz_vy, double *pvybtpp,
    double *vz, int BD_nx_vz, int BD_ny_vz, int BD_nz_vz, double *pvzbtpp,
    double *tpp, double *ptppbvx, double *ptppbvy, double *ptppbvz,
    int BD_nx_tpp, int BD_ny_tpp, int BD_nz_tpp,
    double *rho, double *lambda, double *fdc,
    int nT, double dt, double dx, double dy, double dz, int ext,
    double *bhalf, double *ahalf, double *bfull, double *afull,
    char *snp_path, long long dim, int intvl,
    long long *threadim, long long *blockdim,
    int AssignedStreamNum, int TotalStreamNum, int RegStreamDim,
    int *vx_PV_start, int *vy_PV_start, int *vz_PV_start,
    int *tpp_PV_start, int *rho_PV_start,
    int *vx_PV_offset, int *vy_PV_offset, int *vz_PV_offset,
    int *tpp_PV_offset, int *rho_PV_offset,
    int *vx_SS_start, int *vy_SS_start, int *vz_SS_start,
    int *tpp_SS_start, int *lambda_SS_start,
    int *vx_SS_offset, int *vy_SS_offset, int *vz_SS_offset,
    int *tpp_SS_offset, int *lambda_SS_offset)
  {

    int t,nstream;
    cudaStream_t stream[AssignedStreamNum];
    for(nstream=0 ; nstream < AssignedStreamNum ; nstream++){
      cudaStreamCreate(&stream[nstream]); // create concurrent streams
    }

    // double *snp_ptr = mmap_snapshot(snp_path, dim);

    // copy finite difference coefficient to constant memory
    cudaMemcpyToSymbol(dev_fdc,fdc,4*sizeof(double));

    // number of threads and blocks used
    dim3 threads(*threadim,*(threadim+1),*(threadim+2));
    dim3 blocks(*blockdim,*(blockdim+1),*(blockdim+2));

    // pinned host memory for faster transfer between host and device mem
    acHostRegister(vx, BD_nx_vx, BD_ny_vx, BD_nz_vx,
    vy, BD_nx_vy, BD_ny_vy, BD_nz_vy,
    vz, BD_nx_vz, BD_ny_vz, BD_nz_vz,
    tpp, BD_nx_tpp, BD_ny_tpp, BD_nz_tpp,
    rho, lambda);


    // Assign device space for each stream
    double *dev_vx[AssignedStreamNum], *dev_vy[AssignedStreamNum], *dev_vz[AssignedStreamNum];
    double *dev_tpp[AssignedStreamNum];
    double *dev_rho[AssignedStreamNum], *dev_lambda[AssignedStreamNum];

    acDeviceMalloc(dev_vx, BD_nx_vx, BD_nz_vx,
      dev_vy, BD_nx_vy, BD_nz_vy,
      dev_vz, BD_nx_vz, BD_nz_vz,
      dev_tpp, BD_nx_tpp, BD_nz_tpp,
      dev_rho, dev_lambda,
      AssignedStreamNum, RegStreamDim);



    //************ time iteration ************//
    for(t=0;t<nT;t++)
    {
      // Particle velocities part
      // Copy data from host memory to device memory
      acMemcpyHToDforParticleVel(
        dev_vx, vx, BD_nx_vx, BD_nz_vx,
        dev_vy, vy, BD_nx_vy, BD_nz_vy,
        dev_vz, vz, BD_nx_vz, BD_nz_vz,
        dev_tpp, tpp, BD_nx_tpp, BD_nz_tpp,
        dev_rho, rho,
        TotalStreamNum, AssignedStreamNum, stream,
        vx_PV_start, vy_PV_start, vz_PV_start,
        tpp_PV_start, rho_PV_start,
        vx_PV_offset, vy_PV_offset, vz_PV_offset,
        tpp_PV_offset, rho_PV_offset);

      // Kernekl Execution for particle velocities
      acKernelExecforParticleVel(
        dev_vx, vx, BD_nx_vx, BD_nz_vx,
        dev_vy, vy, BD_nx_vy, BD_nz_vy,
        dev_vz, vz, BD_nx_vz, BD_nz_vz,
        dev_tpp, tpp, BD_nx_tpp, BD_nz_tpp,
        dev_rho, rho,
        dx, dy, dz, dt,
        blocks, threads,
        TotalStreamNum, AssignedStreamNum, stream,
        vx_PV_offset, vy_PV_offset, vz_PV_offset,
        vx_PV_start, vy_PV_start, vz_PV_start,
        tpp_PV_start, rho_PV_start);

      // Copy data from device meory to host memory
      acMemcpyDToHforParticleVel(
        dev_vx, vx, BD_nx_vx, BD_nz_vx,
        dev_vy, vy, BD_nx_vy, BD_nz_vy,
        dev_vz, vz, BD_nx_vz, BD_nz_vz,
        TotalStreamNum, AssignedStreamNum, stream,
        vx_PV_start, vy_PV_start, vz_PV_start,
        vx_PV_offset, vy_PV_offset, vz_PV_offset);

      // finish computing particle velocities before computing stress
      cudaDeviceSynchronize();

      // Stress part
      // Copy data from host memory to device memory
      acMemcpyHToDforStress(
        dev_vx, vx, BD_nx_vx, BD_nz_vx,
        dev_vy, vy, BD_nx_vy, BD_nz_vy,
        dev_vz, vz, BD_nx_vz, BD_nz_vz,
        dev_tpp, tpp, BD_nx_tpp, BD_nz_tpp,
        dev_lambda, lambda,
        TotalStreamNum, AssignedStreamNum, stream,
        vx_SS_start, vy_SS_start, vz_SS_start,
        tpp_SS_start, lambda_SS_start,
        vx_SS_offset, vy_SS_offset, vz_SS_offset,
        tpp_SS_offset, lambda_SS_offset);

      // Kernekl Execution for stress
      acKernelExecforStress(
        dev_vx, vx, BD_nx_vx, BD_nz_vx,
        dev_vy, vy, BD_nx_vy, BD_nz_vy,
        dev_vz, vz, BD_nx_vz, BD_nz_vz,
        dev_tpp, tpp, BD_nx_tpp, BD_nz_tpp,
        dev_lambda, lambda,
        dx, dy, dz, dt,
        blocks, threads,
        TotalStreamNum, AssignedStreamNum, stream,
        tpp_SS_offset);

      // Copy data from device meory to host memory
      acMemcpyDToHforStress(
        dev_tpp, tpp, BD_nx_tpp, BD_nz_tpp,
        TotalStreamNum, AssignedStreamNum, stream,
        tpp_SS_start, tpp_SS_offset);

      // finish computing stress before going to the next time step
      cudaDeviceSynchronize();

      // output snapshot
      // if(t%intvl==0){
      //   memcpy(snp_ptr+nsnp*BD_nx_vz*BD_ny_vz*BD_nz_vz,vz,BD_nx_vz*BD_ny_vz*BD_nz_vz*sizeof(double));
      //   nsnp++;
      // }
    }

    acHostUnRegister(vx, vy, vz, tpp, rho, lambda);
    acDeviceFree(dev_vx, dev_vy, dev_vz, dev_tpp,
      dev_lambda, dev_rho, AssignedStreamNum, stream);
    // munmap_snapshot(snp_ptr, snp_path);
  }
}
