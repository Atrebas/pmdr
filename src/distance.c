#include <omp.h>
#include <R.h>
#include <Rinternals.h>


// The c_hello function below has been written by Drew Schmidt
// https://raw.githubusercontent.com/wrathematics/Romp/master/src/examples_c.c


SEXP c_hello()
{
  int tid, nthreads;
  
  #pragma omp parallel private(tid)
  {
    nthreads = omp_get_num_threads();
    tid = omp_get_thread_num();
    
    Rprintf("Hello from thread %d of %d\n", tid, nthreads);
  }
  
  return R_NilValue;
}



// ------------------------------------------------------
// The three functions below are used to compute manhattan distances
// with different looping strategies
// Note that the code handles INTEGER VALUES ONLY
// Numeric or missing values are not allowed

// These functions take an integer matrix as input
// and return an integer vector of distances for all the pairwise combinations
// corresponding the lower triangular matrix
// ------------------------------------------------------


// Method 1: Using two loops on i1 and i2
// i1 and i2 are the row indices in the input matrix
// computation is rowwise, results are colwise

SEXP c_dist_two_loops(SEXP x)
{
  const int I = nrows(x), J = ncols(x);
  const int K = I * (I - 1) / 2;
  SEXP res;
  PROTECT(res = allocVector(INTSXP, K));
  int *ptx = INTEGER(x), *ptres = INTEGER(res);
  int i1, i2, j, k; 

#pragma omp parallel for default(none) \
                         shared(ptx, ptres) private(i1, i2, j, k)

  for (i1 = 1; i1 < I; i1++) {
    for (i2 = 0; i2 < i1; i2++) {
      k = K - (I - i2) * ((I - i2) - 1) / 2 + i1 - i2 - 1;
      ptres[k] = 0;
      for (j = 0; j < J; j++) {
        ptres[k] += abs(ptx[i1 + I * j] - ptx[i2 + I * j]);
      }
    }
  }
  UNPROTECT(1);
  return res;
}



// Method 2: Using one loop on k
// k is a linear index referring to the pairwise combinations (0:(K-1))
// i1 and i2 are then computed from k with colwise numbering

SEXP c_dist_one_loop_colwise(SEXP x)
{
  const int I = nrows(x), J = ncols(x);
  const int K = I * (I - 1) / 2;
  SEXP res;
  PROTECT(res = allocVector(INTSXP, K));
  int *ptx = INTEGER(x), *ptres = INTEGER(res);
  int i1, i2, j, k, kp, p; 

#pragma omp parallel for default(none) \
                         shared(ptx, ptres) private(i1, i2, j, k, kp, p)

  for (k = 0; k < K; k++) {
    ptres[k] = 0;
    kp = K - k - 1;
    p  = floor((sqrt(1 + 8 * kp) - 1) / 2);
    i1 = I - (kp - p * (p + 1) / 2) - 1;
    i2 = I - 2 - p;
      for (j = 0; j < J; j++) {
        ptres[k] += abs(ptx[i1 + I * j] - ptx[i2 + I * j]);
      }
  }
  UNPROTECT(1);
  return res;
}



// Method 3: Using one loop on k
// k is a linear index referring to the pairwise combinations (0:(K-1))
// i1 and i2 are then computed from k with diagwise numbering

SEXP c_dist_one_loop_diagwise(SEXP x)
{
  const int I = nrows(x), J = ncols(x);
  const int K = I * (I - 1) / 2;
  SEXP res;
  PROTECT(res = allocVector(INTSXP, K));
  int *ptx = INTEGER(x), *ptres = INTEGER(res);
  int i1, i2, j, k, kp, p; 

#pragma omp parallel for default(none) \
                         shared(ptx, ptres) private(i1, i2, j, k, kp, p)

  for (k = 0; k < K; k++) {
    ptres[k] = 0;
    kp = K - k - 1;
    p  = floor((sqrt(1 + 8 * kp) - 1) / 2);
    i1 = I - (kp - p * (p + 1) / 2) - 1;
    i2 = i1 - (I - 1 - p);
      for (j = 0; j < J; j++) {
        ptres[k] += abs(ptx[i1 + I * j] - ptx[i2 + I * j]);
      }
  }
  UNPROTECT(1);
  return res;
}
