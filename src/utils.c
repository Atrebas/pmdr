#include <R.h>
#include <Rinternals.h>
#include <omp.h>


// The code below has been written by Drew Schmidt
// https://raw.githubusercontent.com/wrathematics/Romp/master/src/utils.c


SEXP R_num_procs()
{
  SEXP ret;
  PROTECT(ret = allocVector(INTSXP, 1));

  INTEGER(ret)[0] = omp_get_num_procs();

  UNPROTECT(1);
  return ret;
}


SEXP R_omp_num_threads(SEXP n)
{
  omp_set_num_threads(INTEGER(n)[0]);

  return R_NilValue;
}
