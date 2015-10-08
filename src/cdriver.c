#include <stdlib.h>
#include <stdio.h>
#include "ral_nlls.h"

/* A c driver for the ral_nlls program */
int main(void) {
  
  /* Problem data */
  const int n = 2;
  const int m = 67;
  
  /* Derived types */
  struct NLLS_control_type options;
  struct NLLS_inform_type status;

  int len_work_int = n;
  int *Work_int[len_work_int];
  /*  int *Work_int;
      %  *Work_int = (int*) malloc( len_work_int * sizeof(int) );*/
  int len_work_real = n;
  double *Work_real[len_work_real];
  /*  double *Work_real;
   *Work_real = (double*) malloc( len_work_real * sizeof(double) );*/

  /*  double *X;
   *X = (double*) malloc( n );*/
  double *X[n];
  (*X)[0] = 1.0;
  (*X)[1] = 2.0;

  NLLS_default_control(&options);
  
  RAL_NLLS(n, m, *X, *Work_int, len_work_int, 
	   *Work_real, len_work_real,
	   eval_F, eval_J,
	   &status, &options);

  free(X); free(Work_int); free(Work_real);

  printf("Seems to be working \n");
  
  return 0; /* success */
}

