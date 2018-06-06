__host__ void acKernelExecforParticleVel(
  double *dev_vx[], double *vx, int BD_nx_vx, int BD_nz_vx,
  double *dev_vy[], double *vy, int BD_nx_vy, int BD_nz_vy,
  double *dev_vz[], double *vz, int BD_nx_vz, int BD_nz_vz,
  double *dev_tpp[], double *tpp, int BD_nx_tpp, int BD_nz_tpp,
  double *dev_rho[], double *rho,
  double dx, double dy, double dz, double dt,
  dim3 blocks, dim3 threads,
  int TotalStreamNum, int AssignedStreamNum, cudaStream_t stream[],
  int *vx_PV_offset, int *vy_PV_offset, int *vz_PV_offset,
  int *vx_PV_start, int *vy_PV_start, int *vz_PV_start,
  int *tpp_PV_start, int *rho_PV_start)
  {
    int nstream;
    for(nstream=0; nstream < TotalStreamNum; nstream++){
      ackernel_v<<<blocks,threads,0,stream[nstream%AssignedStreamNum]>>>
        (dev_vx[nstream%AssignedStreamNum], BD_nx_vx, BD_nz_vx,
         dev_vy[nstream%AssignedStreamNum], BD_nx_vy, BD_nz_vy,
         dev_vz[nstream%AssignedStreamNum], BD_nx_vz, BD_nz_vz,
         dev_tpp[nstream%AssignedStreamNum], BD_nx_tpp, BD_nz_tpp,
         dev_rho[nstream%AssignedStreamNum], dx, dy, dz, dt,
         *(vx_PV_offset+nstream),*(vy_PV_offset+nstream),*(vz_PV_offset+nstream),
       *(vx_PV_start+nstream),*(vy_PV_start+nstream),*(vz_PV_start+nstream),
     *(tpp_PV_start+nstream),*(rho_PV_start+nstream));
    }
  }


__host__ void acKernelExecforStress(
  double *dev_vx[], double *vx, int BD_nx_vx, int BD_nz_vx,
  double *dev_vy[], double *vy, int BD_nx_vy, int BD_nz_vy,
  double *dev_vz[], double *vz, int BD_nx_vz, int BD_nz_vz,
  double *dev_tpp[], double *tpp, int BD_nx_tpp, int BD_nz_tpp,
  double *dev_lambda[], double *lambda,
  double dx, double dy, double dz, double dt,
  dim3 blocks, dim3 threads,
  int TotalStreamNum, int AssignedStreamNum, cudaStream_t stream[],
  int *tpp_SS_offset)
  {
    int nstream;
    for(nstream=0; nstream < TotalStreamNum; nstream++){
      ackernel_tau<<<blocks,threads,0,stream[nstream%AssignedStreamNum]>>>
        (dev_vx[nstream%AssignedStreamNum], BD_nx_vx, BD_nz_vx,
          dev_vy[nstream%AssignedStreamNum], BD_nx_vy, BD_nz_vy,
          dev_vz[nstream%AssignedStreamNum], BD_nx_vz, BD_nz_vz,
          dev_tpp[nstream%AssignedStreamNum], BD_nx_tpp, BD_nz_tpp,
          dev_lambda[nstream%AssignedStreamNum], dx, dy, dz, dt, *(tpp_SS_offset+nstream));
    }
  }
