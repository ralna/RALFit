#include <stdlib.h>
#include <stdio.h>
#include "ral_nlls.h"

/* A c driver for the ral_nlls program */
int main(void) {
  
  /* Problem data */
  const int n = 2;
  const int m = 67;
  
  /* Derived types */
  struct nlls_control_type options;
  struct nlls_inform_type status;

  printf("===============\n");
  printf("RAL NLLS driver\n");
  printf("~  C version  ~\n");
  printf("===============\n");

  /*  double *X;
   *X = (double*) malloc( n );*/
  double X[n];
  X[0] = 1.0;
  X[1] = 2.0;

  nlls_default_control(&options);
  
  options.print_level = 3;
  
  ral_nlls_int_func(n, m, X, 
	   &status, &options);

  int i;
  printf("\nX = \n");
  for(i=0; i < n; i++) {
    printf("  %5.4f \n",X[i]);
  }
  
  return 0; /* success */
}

