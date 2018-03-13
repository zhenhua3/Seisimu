#include<stdio.h>
 //********************************************//
 //** bound and body are two functions that  **//
 //** are used in CPU parallel               **//
 //** bound : PML boundary part              **//
 //** body : real medium part                **//
 //********************************************//
double fd3d(double *in_wf, int idx0_in_wf, int idx1_in_wf, int idx2_in_wf, int idx3_in_wf, double *fdc)
 {
  double tmp;
  tmp = *(in_wf+idx0_in_wf)* *(fdc+0) + *(in_wf+idx1_in_wf)* *(fdc+1) + *(in_wf+idx2_in_wf)* *(fdc+2) + *(in_wf+idx3_in_wf)* *(fdc+3);
  return tmp;
 }

  double body_x_3d(double *out_wf, int BD_nz_out, int BD_nx_out, int BD_ny_out,
                  int idx_i_out, int idx_j_out, int idx_k_out,
                  double *in_wf, int BD_nz_in, int BD_nx_in,
                  int idx_j_in_start, double para, double dx,
                  double dt, double *fdc)
 {
  double tmp;
  int idx_out_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf;
  idx_out_wf = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
  idx0_in_wf = idx_i_out + (idx_j_out - idx_j_in_start) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
  idx1_in_wf = idx_i_out + (idx_j_out - idx_j_in_start + 1) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
  idx2_in_wf = idx_i_out + (idx_j_out - idx_j_in_start + 2) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
  idx3_in_wf = idx_i_out + (idx_j_out - idx_j_in_start + 3) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;

  tmp = fd3d(in_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf, fdc);
  tmp = tmp / dx * para * dt;
  return tmp;
 }


   double body_y_3d(double *out_wf, int BD_nz_out, int BD_nx_out, int BD_ny_out,
                   int idx_i_out, int idx_j_out, int idx_k_out,
                   double *in_wf, int BD_nz_in, int BD_nx_in,
                   int idx_k_in_start, double para, double dy,
                   double dt, double *fdc)
  {
   double tmp;
   int idx_out_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf;
   idx_out_wf = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
   idx0_in_wf = idx_i_out + idx_j_out * BD_nz_in + (idx_k_out - idx_k_in_start) * BD_nz_in * BD_nx_in;
   idx1_in_wf = idx_i_out + idx_j_out * BD_nz_in + (idx_k_out - idx_k_in_start + 1) * BD_nz_in * BD_nx_in;
   idx2_in_wf = idx_i_out + idx_j_out * BD_nz_in + (idx_k_out - idx_k_in_start + 2) * BD_nz_in * BD_nx_in;
   idx3_in_wf = idx_i_out + idx_j_out * BD_nz_in + (idx_k_out - idx_k_in_start + 3) * BD_nz_in * BD_nx_in;

   tmp = fd3d(in_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf, fdc);
   tmp = para * tmp / dy * dt;
   return tmp;
  }


  double body_z_3d(double *out_wf, int BD_nz_out, int BD_nx_out, int BD_ny_out,
                  int idx_i_out, int idx_j_out, int idx_k_out,
                  double *in_wf, int BD_nz_in, int BD_nx_in,
                  int idx_i_in_start, double para, double dz,
                  double dt, double *fdc)
 {
  double tmp;
  int idx_out_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf;
  idx_out_wf = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
  idx0_in_wf = (idx_i_out - idx_i_in_start) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
  idx1_in_wf = (idx_i_out - idx_i_in_start + 1) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
  idx2_in_wf = (idx_i_out - idx_i_in_start + 2) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
  idx3_in_wf = (idx_i_out - idx_i_in_start + 3) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;

  tmp = fd3d(in_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf, fdc);
  tmp = para * tmp / dz * dt;
  return tmp;
 }

  void bound_x_3d(double *out_wf, int BD_nz_out, int BD_nx_out, int BD_ny_out,
                  int idx_i_out, int idx_j_out, int idx_k_out,
                  double *in_wf, int BD_nz_in, int BD_nx_in,
                  int idx_j_in_start,
                  double para_head, double para_toe, double dx,
                  double dt,
                  double *BDwf, double *BDcoeff_b, double *BDcoeff_a, int ext,
                  double *fdc)
  {
   int idx_BDcoeff = idx_j_out;
   double tmp_head, tmp_toe;
   int idx_out_wf_head, idx_out_wf_toe, idx_BDwf_head, idx_BDwf_toe;
   tmp_head = body_x_3d(out_wf, BD_nz_out, BD_nx_out, BD_ny_out,
          idx_i_out, idx_j_out, idx_k_out,
          in_wf, BD_nz_in, BD_nx_in, idx_j_in_start, para_head, dx, dt, fdc);

   idx_out_wf_head = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
   idx_BDwf_head = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * 2 * ext;
   *(BDwf+idx_BDwf_head) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf+idx_BDwf_head) + *(BDcoeff_a+idx_BDcoeff) * tmp_head;
   *(out_wf+idx_out_wf_head) = *(out_wf+idx_out_wf_head) + *(BDwf+idx_BDwf_head);

   tmp_toe = body_x_3d(out_wf, BD_nz_out, BD_nx_out,BD_ny_out,
          idx_i_out, BD_nx_out - idx_j_out - 1, idx_k_out,
          in_wf, BD_nz_in, BD_nx_in, idx_j_in_start, para_toe, dx, dt, fdc);

   idx_out_wf_toe = idx_i_out + (BD_nx_out - 1 - idx_j_out) * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
   idx_BDwf_toe = idx_i_out + (2 * ext - 1 - idx_j_out) * BD_nz_out + idx_k_out * BD_nz_out * 2 * ext;
   *(BDwf+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe;
   *(out_wf+idx_out_wf_toe) = *(out_wf+idx_out_wf_toe) + *(BDwf+idx_BDwf_toe);
  }

   void bound_y_3d(double *out_wf, int BD_nz_out, int BD_nx_out, int BD_ny_out,
                   int idx_i_out, int idx_j_out, int idx_k_out,
                   double *in_wf, int BD_nz_in, int BD_nx_in,
                   int idx_k_in_start,
                   double para_head, double para_toe, double dy,
                   double dt,
                   double *BDwf, double *BDcoeff_b, double *BDcoeff_a, int ext,
                   double *fdc)
   {
    int idx_BDcoeff = idx_k_out;
    double tmp_head, tmp_toe;
    int idx_out_wf_head, idx_out_wf_toe, idx_BDwf_head, idx_BDwf_toe;
    tmp_head = body_y_3d(out_wf, BD_nz_out, BD_nx_out, BD_ny_out,
           idx_i_out, idx_j_out, idx_k_out,
           in_wf, BD_nz_in, BD_nx_in, idx_k_in_start, para_head, dy, dt, fdc);

    idx_out_wf_head = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
    idx_BDwf_head = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
    *(BDwf+idx_BDwf_head) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf+idx_BDwf_head) + *(BDcoeff_a+idx_BDcoeff) * tmp_head;
    *(out_wf+idx_out_wf_head) = *(out_wf+idx_out_wf_head) + *(BDwf+idx_BDwf_head);

    tmp_toe = body_y_3d(out_wf, BD_nz_out, BD_nx_out, BD_ny_out,
           idx_i_out, idx_j_out, BD_ny_out - idx_k_out - 1,
           in_wf, BD_nz_in, BD_nx_in, idx_k_in_start, para_toe, dy, dt, fdc);

    idx_out_wf_toe = idx_i_out + idx_j_out * BD_nz_out + (BD_ny_out - 1 - idx_k_out) * BD_nz_out * BD_nx_out;
    idx_BDwf_toe = idx_i_out + idx_j_out * BD_nz_out + (2 * ext - 1 - idx_k_out) * BD_nz_out * BD_nx_out;
    *(BDwf+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe;
    *(out_wf+idx_out_wf_toe) = *(out_wf+idx_out_wf_toe) + *(BDwf+idx_BDwf_toe);
   }

    void unlimited_bound_z_3d(double *out_wf, int BD_nz_out, int BD_nx_out, int BD_ny_out,
                    int idx_i_out, int idx_j_out, int idx_k_out,
                    double *in_wf, int BD_nz_in, int BD_nx_in,
                    int idx_i_in_start,
                    double para_head, double para_toe, double dz,
                    double dt,
                    double *BDwf, double *BDcoeff_b, double *BDcoeff_a, int ext,
                    double *fdc)
    {
     int idx_BDcoeff = idx_i_out;
     double tmp_head, tmp_toe;
     int idx_out_wf_head, idx_out_wf_toe, idx_BDwf_head, idx_BDwf_toe;
     tmp_head = body_z_3d(out_wf, BD_nz_out, BD_nx_out, BD_ny_out,
            idx_i_out, idx_j_out, idx_k_out,
            in_wf, BD_nz_in, BD_nx_in, idx_i_in_start, para_head, dz, dt, fdc);

     idx_out_wf_head = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
     idx_BDwf_head = idx_i_out + idx_j_out * 2 * ext + idx_k_out * 2 * ext * BD_nx_out;
     *(BDwf+idx_BDwf_head) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf+idx_BDwf_head) + *(BDcoeff_a+idx_BDcoeff) * tmp_head ;
     *(out_wf+idx_out_wf_head) = *(out_wf+idx_out_wf_head) + *(BDwf+idx_BDwf_head);

     tmp_toe = body_z_3d(out_wf, BD_nz_out, BD_nx_out, BD_ny_out,
            BD_nz_out - idx_i_out - 1, idx_j_out, idx_k_out,
            in_wf, BD_nz_in, BD_nx_in, idx_i_in_start, para_toe, dz, dt, fdc);

     idx_out_wf_toe = (BD_nz_out - 1 - idx_i_out) + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
     idx_BDwf_toe = (2 * ext - 1 - idx_i_out) + idx_j_out * 2 * ext + idx_k_out * 2 * ext * BD_nx_out;
     *(BDwf+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe;
     *(out_wf+idx_out_wf_toe) = *(out_wf+idx_out_wf_toe) + *(BDwf+idx_BDwf_toe);
    }

     void free_bound_z_3d(double *out_wf, int BD_nz_out, int BD_nx_out, int BD_ny_out,
                    int idx_i_out, int idx_j_out, int idx_k_out,
                    double *in_wf, int BD_nz_in, int BD_nx_in,
                    int idx_i_in_start, double para_toe, double dz,
                    double dt,
                    double *BDwf, double *BDcoeff_b, double *BDcoeff_a, int ext,
                    double *fdc)
    {
     int idx_BDcoeff = idx_i_out;
     double tmp_toe;
     int idx_out_wf_toe, idx_BDwf_toe;

     tmp_toe = body_z_3d(out_wf, BD_nz_out, BD_nx_out, BD_ny_out,
            BD_nz_out - idx_i_out - 1, idx_j_out, idx_k_out,
            in_wf, BD_nz_in, BD_nx_in, idx_i_in_start, para_toe, dz, dt, fdc);

     idx_out_wf_toe = (BD_nz_out - 1 - idx_i_out) + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
     idx_BDwf_toe = (ext - 1 - idx_i_out) + idx_j_out * ext + idx_k_out * ext * BD_nx_out;
     *(BDwf+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe;
     *(out_wf+idx_out_wf_toe) = *(out_wf+idx_out_wf_toe) + *(BDwf+idx_BDwf_toe);
    }

 void el_bound_tpp_x_3d(double *out_txx, double *out_tyy, double *out_tzz,
                  int BD_nz_out, int BD_nx_out, int BD_ny_out,
                  int idx_i_out, int idx_j_out, int idx_k_out,
                  double *in_wf, int BD_nz_in, int BD_nx_in,
                  int idx_j_in_start, double *lambda, double *mu, double dx,
                  double dt,
                  double *BDwf_txx, double *BDwf_tyy, double *BDwf_tzz,
                  double *BDcoeff_b, double *BDcoeff_a, int ext,
                  double *fdc)
  {
   int idx_BDcoeff = idx_j_out;
   int idx_BDwf_head, idx_BDwf_toe;
   double tmp, tmp_head_txx, tmp_head_tyy, tmp_head_tzz,tmp_toe_txx, tmp_toe_tyy, tmp_toe_tzz;
   int idx_out_head_tpp, idx_out_toe_tpp, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf;
   idx_out_head_tpp = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
   idx0_in_wf = idx_i_out + (idx_j_out - idx_j_in_start) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
   idx1_in_wf = idx_i_out + (idx_j_out - idx_j_in_start + 1) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
   idx2_in_wf = idx_i_out + (idx_j_out - idx_j_in_start + 2) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
   idx3_in_wf = idx_i_out + (idx_j_out - idx_j_in_start + 3) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;

   tmp = fd3d(in_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf, fdc);
   tmp_head_txx = tmp / dx * (*(lambda+idx_out_head_tpp) + 2**(mu+idx_out_head_tpp)) * dt;
   tmp_head_tyy = tmp / dx * *(lambda+idx_out_head_tpp) * dt;
   tmp_head_tzz = tmp_head_tyy;

   idx_BDwf_head = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * 2 * ext;
   *(BDwf_txx+idx_BDwf_head) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_txx+idx_BDwf_head) + *(BDcoeff_a+idx_BDcoeff) * tmp_head_txx;
   *(BDwf_tyy+idx_BDwf_head) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tyy+idx_BDwf_head) + *(BDcoeff_a+idx_BDcoeff) * tmp_head_tyy;
   *(BDwf_tzz+idx_BDwf_head) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tzz+idx_BDwf_head) + *(BDcoeff_a+idx_BDcoeff) * tmp_head_tzz;

   *(out_txx+idx_out_head_tpp) = *(out_txx+idx_out_head_tpp) + *(BDwf_txx+idx_BDwf_head);
   *(out_tyy+idx_out_head_tpp) = *(out_tyy+idx_out_head_tpp) + *(BDwf_tyy+idx_BDwf_head);
   *(out_tzz+idx_out_head_tpp) = *(out_tzz+idx_out_head_tpp) + *(BDwf_tzz+idx_BDwf_head);

   idx_out_toe_tpp = idx_i_out + (BD_nx_out - 1 - idx_j_out) * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
   idx0_in_wf = idx_i_out + (BD_nx_out - 1 - idx_j_out - idx_j_in_start) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
   idx1_in_wf = idx_i_out + (BD_nx_out - 1 - idx_j_out - idx_j_in_start + 1) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
   idx2_in_wf = idx_i_out + (BD_nx_out - 1 - idx_j_out - idx_j_in_start + 2) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
   idx3_in_wf = idx_i_out + (BD_nx_out - 1 - idx_j_out - idx_j_in_start + 3) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;

   tmp = fd3d(in_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf, fdc);
   tmp_toe_txx = tmp / dx * (*(lambda+idx_out_toe_tpp) + 2**(mu+idx_out_toe_tpp)) * dt;
   tmp_toe_tyy = tmp / dx * *(lambda+idx_out_toe_tpp) * dt;
   tmp_toe_tzz = tmp_toe_tyy;

   idx_BDwf_toe = idx_i_out + (2 * ext - 1 - idx_j_out) * BD_nz_out + idx_k_out * BD_nz_out * 2 * ext;
   *(BDwf_txx+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_txx+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe_txx;
   *(BDwf_tyy+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tyy+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe_tyy;
   *(BDwf_tzz+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tzz+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe_tzz;

   *(out_txx+idx_out_toe_tpp) = *(out_txx+idx_out_toe_tpp) + *(BDwf_txx+idx_BDwf_toe);
   *(out_tyy+idx_out_toe_tpp) = *(out_tyy+idx_out_toe_tpp) + *(BDwf_tyy+idx_BDwf_toe);
   *(out_tzz+idx_out_toe_tpp) = *(out_tzz+idx_out_toe_tpp) + *(BDwf_tzz+idx_BDwf_toe);
  }

  void el_bound_tpp_y_3d(double *out_txx, double *out_tyy, double *out_tzz,
                   int BD_nz_out, int BD_nx_out, int BD_ny_out,
                   int idx_i_out, int idx_j_out, int idx_k_out,
                   double *in_wf, int BD_nz_in, int BD_nx_in,
                   int idx_k_in_start, double *lambda, double *mu, double dy,
                   double dt,
                   double *BDwf_txx, double *BDwf_tyy, double *BDwf_tzz,
                   double *BDcoeff_b, double *BDcoeff_a, int ext,
                   double *fdc)
   {
    int idx_BDcoeff = idx_k_out;
    int idx_BDwf_head, idx_BDwf_toe;
    double tmp, tmp_head_txx, tmp_head_tyy, tmp_head_tzz, tmp_toe_txx, tmp_toe_tyy, tmp_toe_tzz;
    int idx_out_head_tpp, idx_out_toe_tpp, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf;
    idx_out_head_tpp = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
    idx0_in_wf = idx_i_out + idx_j_out * BD_nz_in + (idx_k_out  - idx_k_in_start)* BD_nz_in * BD_nx_in;
    idx1_in_wf = idx_i_out + idx_j_out * BD_nz_in + (idx_k_out - idx_k_in_start + 1) * BD_nz_in * BD_nx_in;
    idx2_in_wf = idx_i_out + idx_j_out * BD_nz_in + (idx_k_out - idx_k_in_start + 2) * BD_nz_in * BD_nx_in;
    idx3_in_wf = idx_i_out + idx_j_out * BD_nz_in + (idx_k_out - idx_k_in_start + 3) * BD_nz_in * BD_nx_in;

    tmp = fd3d(in_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf, fdc);
    tmp_head_txx = tmp / dy * *(lambda+idx_out_head_tpp) * dt;
    tmp_head_tyy = tmp / dy * (*(lambda+idx_out_head_tpp) + 2**(mu+idx_out_head_tpp)) * dt;
    tmp_head_tzz = tmp_head_txx;

    idx_BDwf_head = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
    *(BDwf_txx + idx_BDwf_head) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_txx+idx_BDwf_head) + *(BDcoeff_a+idx_BDcoeff) * tmp_head_txx;
    *(BDwf_tyy + idx_BDwf_head) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tyy+idx_BDwf_head) + *(BDcoeff_a+idx_BDcoeff) * tmp_head_tyy;
    *(BDwf_tzz + idx_BDwf_head) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tzz+idx_BDwf_head) + *(BDcoeff_a+idx_BDcoeff) * tmp_head_tzz;

    *(out_txx+idx_out_head_tpp) = *(out_txx+idx_out_head_tpp) + *(BDwf_txx+idx_BDwf_head);
    *(out_tyy+idx_out_head_tpp) = *(out_tyy+idx_out_head_tpp) + *(BDwf_tyy+idx_BDwf_head);
    *(out_tzz+idx_out_head_tpp) = *(out_tzz+idx_out_head_tpp) + *(BDwf_tzz+idx_BDwf_head);

    idx_out_toe_tpp = idx_i_out + idx_j_out * BD_nz_out + (BD_ny_out - 1 - idx_k_out) * BD_nz_out * BD_nx_out;
    idx0_in_wf = idx_i_out + idx_j_out * BD_nz_in + (BD_ny_out - 1 - idx_k_out - idx_k_in_start) * BD_nz_in * BD_nx_in;
    idx1_in_wf = idx_i_out + idx_j_out * BD_nz_in + (BD_ny_out - 1 - idx_k_out - idx_k_in_start + 1) * BD_nz_in * BD_nx_in;
    idx2_in_wf = idx_i_out + idx_j_out * BD_nz_in + (BD_ny_out - 1 - idx_k_out - idx_k_in_start + 2) * BD_nz_in * BD_nx_in;
    idx3_in_wf = idx_i_out + idx_j_out * BD_nz_in + (BD_ny_out - 1 - idx_k_out - idx_k_in_start + 3) * BD_nz_in * BD_nx_in;

    tmp = fd3d(in_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf, fdc);
    tmp_toe_txx = tmp / dy * *(lambda+idx_out_toe_tpp) * dt;
    tmp_toe_tyy = tmp / dy * (*(lambda+idx_out_toe_tpp) + 2**(mu+idx_out_toe_tpp)) * dt;
    tmp_toe_tzz = tmp_toe_txx;

    idx_BDwf_toe = idx_i_out + idx_j_out * BD_nz_out + (2 * ext - 1 - idx_k_out) * BD_nz_out * BD_nx_out;
    *(BDwf_txx+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_txx+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe_txx;
    *(BDwf_tyy+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tyy+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe_tyy;
    *(BDwf_tzz+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tzz+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe_tzz;

    *(out_txx+idx_out_toe_tpp) = *(out_txx+idx_out_toe_tpp) + *(BDwf_txx+idx_BDwf_toe);
    *(out_tyy+idx_out_toe_tpp) = *(out_tyy+idx_out_toe_tpp) + *(BDwf_tyy+idx_BDwf_toe);
    *(out_tzz+idx_out_toe_tpp) = *(out_tzz+idx_out_toe_tpp) + *(BDwf_tzz+idx_BDwf_toe);
   }

   void el_unlimited_bound_tpp_z_3d(double *out_txx, double *out_tyy, double *out_tzz,
                    int BD_nz_out, int BD_nx_out, int BD_ny_out,
                    int idx_i_out, int idx_j_out, int idx_k_out,
                    double *in_wf, int BD_nz_in, int BD_nx_in,
                    int idx_i_in_start, double *lambda, double *mu, double dz,
                    double dt,
                    double *BDwf_txx, double *BDwf_tyy, double *BDwf_tzz,
                    double *BDcoeff_b, double *BDcoeff_a, int ext,
                    double *fdc)
    {
     int idx_BDcoeff = idx_i_out;
     int idx_BDwf_head, idx_BDwf_toe;
     double tmp, tmp_head_txx, tmp_head_tyy, tmp_head_tzz, tmp_toe_txx, tmp_toe_tyy, tmp_toe_tzz;
     int idx_out_head_tpp, idx_out_toe_tpp, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf;
     idx_out_head_tpp = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
     idx0_in_wf = (idx_i_out - idx_i_in_start) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
     idx1_in_wf = (idx_i_out - idx_i_in_start + 1) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
     idx2_in_wf = (idx_i_out - idx_i_in_start + 2) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
     idx3_in_wf = (idx_i_out - idx_i_in_start + 3) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;

     tmp = fd3d(in_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf, fdc);
     tmp_head_txx = tmp / dz * *(lambda + idx_out_head_tpp) * dt;
     tmp_head_tyy = tmp_head_txx;
     tmp_head_tzz = tmp / dz * (*(lambda + idx_out_head_tpp) + 2* *(mu + idx_out_head_tpp)) * dt;

     idx_BDwf_head = idx_i_out + idx_j_out * 2 * ext + idx_k_out * 2 * ext * BD_nx_out;
     *(BDwf_txx + idx_BDwf_head) = *(BDcoeff_b + idx_BDcoeff)* *(BDwf_txx + idx_BDwf_head) + *(BDcoeff_a + idx_BDcoeff) * tmp_head_txx;
     *(BDwf_tyy + idx_BDwf_head) = *(BDcoeff_b + idx_BDcoeff)* *(BDwf_tyy + idx_BDwf_head) + *(BDcoeff_a + idx_BDcoeff) * tmp_head_tyy;
     *(BDwf_tzz + idx_BDwf_head) = *(BDcoeff_b + idx_BDcoeff)* *(BDwf_tzz + idx_BDwf_head) + *(BDcoeff_a + idx_BDcoeff) * tmp_head_tzz;

     *(out_txx + idx_out_head_tpp) = *(out_txx + idx_out_head_tpp) + *(BDwf_txx + idx_BDwf_head);
     *(out_tyy + idx_out_head_tpp) = *(out_tyy + idx_out_head_tpp) + *(BDwf_tyy + idx_BDwf_head);
     *(out_tzz + idx_out_head_tpp) = *(out_tzz + idx_out_head_tpp) + *(BDwf_tzz + idx_BDwf_head);

     idx_out_toe_tpp = (BD_nz_out - 1 - idx_i_out) + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
     idx0_in_wf = (BD_nz_out - 1 - idx_i_out - idx_i_in_start) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
     idx1_in_wf = (BD_nz_out - 1 - idx_i_out - idx_i_in_start + 1) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
     idx2_in_wf = (BD_nz_out - 1 - idx_i_out - idx_i_in_start + 2) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
     idx3_in_wf = (BD_nz_out - 1 - idx_i_out - idx_i_in_start + 3) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;

     tmp = fd3d(in_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf, fdc);
     tmp_toe_txx = tmp / dz * *(lambda+idx_out_toe_tpp) * dt;
     tmp_toe_tyy = tmp_head_txx;
     tmp_toe_tzz = tmp / dz * (*(lambda+idx_out_toe_tpp) + 2**(mu+idx_out_toe_tpp)) * dt;

     idx_BDwf_toe = (2 * ext - 1 - idx_i_out) + idx_j_out * 2 * ext + idx_k_out * 2 * ext * BD_nx_out;
     *(BDwf_txx+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_txx+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe_txx;
     *(BDwf_tyy+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tyy+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe_tyy;
     *(BDwf_tzz+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tzz+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe_tzz;

     *(out_txx+idx_out_toe_tpp) = *(out_txx+idx_out_toe_tpp) + *(BDwf_txx+idx_BDwf_toe);
     *(out_tyy+idx_out_toe_tpp) = *(out_tyy+idx_out_toe_tpp) + *(BDwf_tyy+idx_BDwf_toe);
     *(out_tzz+idx_out_toe_tpp) = *(out_tzz+idx_out_toe_tpp) + *(BDwf_tzz+idx_BDwf_toe);
    }

    void el_free_bound_tpp_z_3d(double *out_txx, double *out_tyy, double *out_tzz,
                     int BD_nz_out, int BD_nx_out, int BD_ny_out,
                     int idx_i_out, int idx_j_out, int idx_k_out,
                     double *in_wf, int BD_nz_in, int BD_nx_in,
                     int idx_i_in_start, double *lambda, double *mu, double dz,
                     double dt,
                     double *BDwf_txx, double *BDwf_tyy, double *BDwf_tzz,
                     double *BDcoeff_b, double *BDcoeff_a, int ext,
                     double *fdc)
     {
      int idx_BDcoeff = idx_i_out;
      int idx_BDwf_toe;
      double tmp, tmp_toe_txx, tmp_toe_tyy, tmp_toe_tzz;
      int idx_out_toe_tpp, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf;

      idx_out_toe_tpp = (BD_nz_out - 1 - idx_i_out) + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
      idx0_in_wf = (BD_nz_out - 1 - idx_i_out - idx_i_in_start) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
      idx1_in_wf = (BD_nz_out - 1 - idx_i_out - idx_i_in_start + 1) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
      idx2_in_wf = (BD_nz_out - 1 - idx_i_out - idx_i_in_start + 2) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
      idx3_in_wf = (BD_nz_out - 1 - idx_i_out - idx_i_in_start + 3) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;

      tmp = fd3d(in_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf, fdc);
      tmp_toe_txx = tmp / dz * *(lambda+idx_out_toe_tpp) * dt;
      tmp_toe_tyy = tmp_toe_txx;
      tmp_toe_tzz = tmp / dz * (*(lambda+idx_out_toe_tpp) + 2**(mu+idx_out_toe_tpp)) * dt;

      idx_BDwf_toe = (ext - 1 - idx_i_out) + idx_j_out * ext + idx_k_out * ext * BD_nx_out;
      *(BDwf_txx+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_txx+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe_txx;
      *(BDwf_tyy+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tyy+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe_tyy;
      *(BDwf_tzz+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tzz+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe_tzz;

      *(out_txx+idx_out_toe_tpp) = *(out_txx+idx_out_toe_tpp) + *(BDwf_txx+idx_BDwf_toe);
      *(out_tyy+idx_out_toe_tpp) = *(out_tyy+idx_out_toe_tpp) + *(BDwf_tyy+idx_BDwf_toe);
      *(out_tzz+idx_out_toe_tpp) = *(out_tzz+idx_out_toe_tpp) + *(BDwf_tzz+idx_BDwf_toe);
     }

            /*********************/
            /* Acoustic boundary */
            /*********************/
     void ac_bound_tpp_x_3d(double *out_tpp,
                      int BD_nz_out, int BD_nx_out, int BD_ny_out,
                      int idx_i_out, int idx_j_out, int idx_k_out,
                      double *in_wf, int BD_nz_in, int BD_nx_in,
                      int idx_j_in_start, double *lambda, double dx,
                      double dt,
                      double *BDwf_tpp,
                      double *BDcoeff_b, double *BDcoeff_a, int ext,
                      double *fdc)
      {
       int idx_BDcoeff = idx_j_out;
       int idx_BDwf_head, idx_BDwf_toe;
       double tmp, tmp_head_tpp,tmp_toe_tpp;
       int idx_out_head_tpp, idx_out_toe_tpp, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf;
       idx_out_head_tpp = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
       idx0_in_wf = idx_i_out + (idx_j_out - idx_j_in_start) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
       idx1_in_wf = idx_i_out + (idx_j_out - idx_j_in_start + 1) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
       idx2_in_wf = idx_i_out + (idx_j_out - idx_j_in_start + 2) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
       idx3_in_wf = idx_i_out + (idx_j_out - idx_j_in_start + 3) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;

       tmp = fd3d(in_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf, fdc);
       tmp_head_tpp = tmp / dx * *(lambda+idx_out_head_tpp) * dt;

       idx_BDwf_head = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * 2 * ext;
       *(BDwf_tpp+idx_BDwf_head) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tpp+idx_BDwf_head) + *(BDcoeff_a+idx_BDcoeff) * tmp_head_tpp;

       *(out_tpp+idx_out_head_tpp) = *(out_tpp+idx_out_head_tpp)+ *(BDwf_tpp+idx_BDwf_head);

       idx_out_toe_tpp = idx_i_out + (BD_nx_out - 1 - idx_j_out) * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
       idx0_in_wf = idx_i_out + (BD_nx_out - 1 - idx_j_out - idx_j_in_start) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
       idx1_in_wf = idx_i_out + (BD_nx_out - 1 - idx_j_out - idx_j_in_start + 1) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
       idx2_in_wf = idx_i_out + (BD_nx_out - 1 - idx_j_out - idx_j_in_start + 2) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
       idx3_in_wf = idx_i_out + (BD_nx_out - 1 - idx_j_out - idx_j_in_start + 3) * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;

       tmp = fd3d(in_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf, fdc);
       tmp_toe_tpp = tmp / dx * *(lambda+idx_out_toe_tpp) * dt;

       idx_BDwf_toe = idx_i_out + (2 * ext - 1 - idx_j_out) * BD_nz_out + idx_k_out * BD_nz_out * 2 * ext;
       *(BDwf_tpp+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tpp+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe_tpp;

       *(out_tpp+idx_out_toe_tpp) = *(out_tpp+idx_out_toe_tpp) + *(BDwf_tpp+idx_BDwf_toe);
      }

      void ac_bound_tpp_y_3d(double *out_tpp,
                       int BD_nz_out, int BD_nx_out, int BD_ny_out,
                       int idx_i_out, int idx_j_out, int idx_k_out,
                       double *in_wf, int BD_nz_in, int BD_nx_in,
                       int idx_k_in_start, double *lambda, double dy,
                       double dt,
                       double *BDwf_tpp,
                       double *BDcoeff_b, double *BDcoeff_a, int ext,
                       double *fdc)
       {
        int idx_BDcoeff = idx_k_out;
        int idx_BDwf_head, idx_BDwf_toe;
        double tmp, tmp_head_tpp,tmp_toe_tpp;
        int idx_out_head_tpp, idx_out_toe_tpp, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf;
        idx_out_head_tpp = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
        idx0_in_wf = idx_i_out + idx_j_out * BD_nz_in + (idx_k_out  - idx_k_in_start)* BD_nz_in * BD_nx_in;
        idx1_in_wf = idx_i_out + idx_j_out * BD_nz_in + (idx_k_out - idx_k_in_start + 1) * BD_nz_in * BD_nx_in;
        idx2_in_wf = idx_i_out + idx_j_out * BD_nz_in + (idx_k_out - idx_k_in_start + 2) * BD_nz_in * BD_nx_in;
        idx3_in_wf = idx_i_out + idx_j_out * BD_nz_in + (idx_k_out - idx_k_in_start + 3) * BD_nz_in * BD_nx_in;

        tmp = fd3d(in_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf, fdc);
        tmp_head_tpp = tmp / dy * *(lambda+idx_out_head_tpp) * dt;

        idx_BDwf_head = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
        *(BDwf_tpp + idx_BDwf_head) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tpp+idx_BDwf_head) + *(BDcoeff_a+idx_BDcoeff) * tmp_head_tpp;

        *(out_tpp+idx_out_head_tpp) = *(out_tpp+idx_out_head_tpp) + *(BDwf_tpp+idx_BDwf_head);

        idx_out_toe_tpp = idx_i_out + idx_j_out * BD_nz_out + (BD_ny_out - 1 - idx_k_out) * BD_nz_out * BD_nx_out;
        idx0_in_wf = idx_i_out + idx_j_out * BD_nz_in + (BD_ny_out - 1 - idx_k_out - idx_k_in_start) * BD_nz_in * BD_nx_in;
        idx1_in_wf = idx_i_out + idx_j_out * BD_nz_in + (BD_ny_out - 1 - idx_k_out - idx_k_in_start + 1) * BD_nz_in * BD_nx_in;
        idx2_in_wf = idx_i_out + idx_j_out * BD_nz_in + (BD_ny_out - 1 - idx_k_out - idx_k_in_start + 2) * BD_nz_in * BD_nx_in;
        idx3_in_wf = idx_i_out + idx_j_out * BD_nz_in + (BD_ny_out - 1 - idx_k_out - idx_k_in_start + 3) * BD_nz_in * BD_nx_in;

        tmp = fd3d(in_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf, fdc);
        tmp_toe_tpp = tmp / dy * *(lambda+idx_out_toe_tpp) * dt;

        idx_BDwf_toe = idx_i_out + idx_j_out * BD_nz_out + (2 * ext - 1 - idx_k_out) * BD_nz_out * BD_nx_out;
        *(BDwf_tpp+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tpp+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe_tpp;

        *(out_tpp+idx_out_toe_tpp) = *(out_tpp+idx_out_toe_tpp) + *(BDwf_tpp+idx_BDwf_toe);
       }

       void ac_unlimited_bound_tpp_z_3d(double *out_tpp,
                        int BD_nz_out, int BD_nx_out, int BD_ny_out,
                        int idx_i_out, int idx_j_out, int idx_k_out,
                        double *in_wf, int BD_nz_in, int BD_nx_in,
                        int idx_i_in_start, double *lambda, double dz,
                        double dt,
                        double *BDwf_tpp,
                        double *BDcoeff_b, double *BDcoeff_a, int ext,
                        double *fdc)
        {
         int idx_BDcoeff = idx_i_out;
         int idx_BDwf_head, idx_BDwf_toe;
         double tmp, tmp_head_tpp, tmp_toe_tpp;
         int idx_out_head_tpp, idx_out_toe_tpp, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf;
         idx_out_head_tpp = idx_i_out + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
         idx0_in_wf = (idx_i_out - idx_i_in_start) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
         idx1_in_wf = (idx_i_out - idx_i_in_start + 1) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
         idx2_in_wf = (idx_i_out - idx_i_in_start + 2) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
         idx3_in_wf = (idx_i_out - idx_i_in_start + 3) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;

         tmp = fd3d(in_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf, fdc);
         tmp_head_tpp = tmp / dz * *(lambda + idx_out_head_tpp) * dt;

         idx_BDwf_head = idx_i_out + idx_j_out * 2 * ext + idx_k_out * 2 * ext * BD_nx_out;
         *(BDwf_tpp + idx_BDwf_head) = *(BDcoeff_b + idx_BDcoeff)* *(BDwf_tpp + idx_BDwf_head) + *(BDcoeff_a + idx_BDcoeff) * tmp_head_tpp;

         *(out_tpp + idx_out_head_tpp) = *(out_tpp + idx_out_head_tpp) + *(BDwf_tpp + idx_BDwf_head);

         idx_out_toe_tpp = (BD_nz_out - 1 - idx_i_out) + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
         idx0_in_wf = (BD_nz_out - 1 - idx_i_out - idx_i_in_start) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
         idx1_in_wf = (BD_nz_out - 1 - idx_i_out - idx_i_in_start + 1) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
         idx2_in_wf = (BD_nz_out - 1 - idx_i_out - idx_i_in_start + 2) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
         idx3_in_wf = (BD_nz_out - 1 - idx_i_out - idx_i_in_start + 3) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;

         tmp = fd3d(in_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf, fdc);
         tmp_toe_tpp = tmp / dz * *(lambda+idx_out_toe_tpp) * dt;

         idx_BDwf_toe = (2 * ext - 1 - idx_i_out) + idx_j_out * 2 * ext + idx_k_out * 2 * ext * BD_nx_out;
         *(BDwf_tpp+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tpp+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe_tpp;

         *(out_tpp+idx_out_toe_tpp) = *(out_tpp+idx_out_toe_tpp) + *(BDwf_tpp+idx_BDwf_toe);
        }

        void ac_free_bound_tpp_z_3d(double *out_tpp,
                         int BD_nz_out, int BD_nx_out, int BD_ny_out,
                         int idx_i_out, int idx_j_out, int idx_k_out,
                         double *in_wf, int BD_nz_in, int BD_nx_in,
                         int idx_i_in_start, double *lambda, double dz,
                         double dt,
                         double *BDwf_tpp,
                         double *BDcoeff_b, double *BDcoeff_a, int ext,
                         double *fdc)
         {
          int idx_BDcoeff = idx_i_out;
          int idx_BDwf_toe;
          double tmp, tmp_toe_tpp;
          int idx_out_toe_tpp, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf;

          idx_out_toe_tpp = (BD_nz_out - 1 - idx_i_out) + idx_j_out * BD_nz_out + idx_k_out * BD_nz_out * BD_nx_out;
          idx0_in_wf = (BD_nz_out - 1 - idx_i_out - idx_i_in_start) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
          idx1_in_wf = (BD_nz_out - 1 - idx_i_out - idx_i_in_start + 1) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
          idx2_in_wf = (BD_nz_out - 1 - idx_i_out - idx_i_in_start + 2) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;
          idx3_in_wf = (BD_nz_out - 1 - idx_i_out - idx_i_in_start + 3) + idx_j_out * BD_nz_in + idx_k_out * BD_nz_in * BD_nx_in;

          tmp = fd3d(in_wf, idx0_in_wf, idx1_in_wf, idx2_in_wf, idx3_in_wf, fdc);
          tmp_toe_tpp = tmp / dz * *(lambda+idx_out_toe_tpp) * dt;

          idx_BDwf_toe = (ext - 1 - idx_i_out) + idx_j_out * ext + idx_k_out * ext * BD_nx_out;
          *(BDwf_tpp+idx_BDwf_toe) = *(BDcoeff_b+idx_BDcoeff)* *(BDwf_tpp+idx_BDwf_toe) + *(BDcoeff_a+idx_BDcoeff) * tmp_toe_tpp;

          *(out_tpp+idx_out_toe_tpp) = *(out_tpp+idx_out_toe_tpp) + *(BDwf_tpp+idx_BDwf_toe);
         }
