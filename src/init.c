/*
 *  Part of R package npsp
 *  Copyright (C) 2012-2017 R. Fernandez-Casal
 */
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
/* .Fortran calls */
extern void F77_NAME(bin_den)(int *nd, int *nbin, double *x, int *ny, 
    double *bin_min, double *bin_max, double *bin_w, int *itype);
extern void F77_NAME(binning_r)(int *nd, int *nbin, double *x, int *ny, 
    double *y, double *bin_min, double *bin_max, double *bin_med, 
    double *bin_y, double *bin_w);
extern void F77_NAME(disc_sbv)(int *nx, double *x, int *dim, double *rango);
extern void F77_NAME(dnrm2_r)(int *nd, double *x, int *nx, double *z, int *nz, 
    double *dist);
extern void F77_NAME(dposv_r)(int *N, int *NRHS, double *A, double *B, int *INFO);
extern void F77_NAME(interp_data_grid)(int *nd, int *nbin, double *bin_min, 
    double *bin_max, int *ngrid, double *gy, double *x, int *ny, double *y);
extern void F77_NAME(lp_bin)(int *nd, int *nbin, int *ntbin, double *bin_min, 
    double *bin_max, double *bin_med, double *bin_y, double *bin_w, double *h, 
    double *lpe, int *degree, int *ideriv, double *deriv, int *ihat, 
    double *hatlp, int *ncv, double *rm, double *rss, int *nrl0);
extern void F77_NAME(lp_data_grid)(int *nd, int *nbin, int *ntbin, double *bin_min, 
    double *bin_max, double *bin_med, double *bin_y, double *h, double *lpe, 
    int *degree, int *ideriv, double *deriv, int *ihat, double *hatlp, 
    int *ncv, double *rm, double *rss, int *nrl0);
extern void F77_NAME(lp_raw)(int *nd, int *nbin, int *ntbin, double *x, int *ny, 
    double *y, double *bin_min, double *bin_max, double *bin_med, double *bin_y, 
    double *bin_w, double *h, double *lpe, int *degree, int *ideriv, 
    double *deriv, int *ihat, double *hatlp, int *ncv, double *rm, 
    double *rss, int *nrl0);
extern void F77_NAME(predict_locpol)(int *nd, int *nbin, double *bin_min, 
    double *bin_max, double *bin_med, double *bin_y, double *bin_w, int *ngrid, 
    double *lpe, int *ihat, double *hatlp, double *x, int *ny, double *lpy, 
    double *hatlpy);
extern void F77_NAME(svar_iso_bin)(int *nd, double *x, int *ny, double *y, 
    int *nlags, double *minlag, double *maxlag, int *itipo, double *bin_lag, 
    double *bin_med, double *bin_y, double *bin_w);
extern void F77_NAME(svar_iso_np)(int *nd, double *x, int *ny, double *y, 
    int *nlags, double *minlag, double *maxlag, double *bin_lag, 
    double *bin_med, double *bin_y, double *bin_w, double *h, double *lpe, 
    int *degree, int *ideriv, double *deriv, int *ihat, double *hatlp, 
    int *ndelcv, double *rm, double *rss, int *nrl0);

static const R_FortranMethodDef FortranEntries[] = {
  {"bin_den",          (DL_FUNC) &F77_NAME(bin_den),           8},
  {"binning_r",        (DL_FUNC) &F77_NAME(binning_r),        10},
  {"disc_sbv",         (DL_FUNC) &F77_NAME(disc_sbv),          4},
  {"dnrm2_r",          (DL_FUNC) &F77_NAME(dnrm2_r),           6},
  {"dposv_r",          (DL_FUNC) &F77_NAME(dposv_r),           5},
  {"interp_data_grid", (DL_FUNC) &F77_NAME(interp_data_grid),  9},
  {"lp_bin",           (DL_FUNC) &F77_NAME(lp_bin),           19},
  {"lp_data_grid",     (DL_FUNC) &F77_NAME(lp_data_grid),     18},
  {"lp_raw",           (DL_FUNC) &F77_NAME(lp_raw),           22},
  {"predict_locpol",   (DL_FUNC) &F77_NAME(predict_locpol),   15},
  {"svar_iso_bin",     (DL_FUNC) &F77_NAME(svar_iso_bin),     12},
  {"svar_iso_np",      (DL_FUNC) &F77_NAME(svar_iso_np),      22},
  {NULL, NULL, 0}
};

void R_init_npsp(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}