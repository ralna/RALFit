// test/nlls_c_test.c
//
// Test basic functionality of the c interface, plus any C-specific routines
#include "ral_nlls.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct params_type {
  double *t; // t_i
  double *y; // y__i
};

// use the test example r_i(x; t_i, y_i) = x_1 e^(x_2 * t_i) - y_i

int eval_r(int n, int m, void * params, double const* x, double* r) {
   double x1 = x[0];
   double x2 = x[1];
   double const* t = ((struct params_type*) params)->t;
   double const* y = ((struct params_type*) params)->y;

   for(int i=0; i<m; i++)
      r[i] = x1 * exp(x2*t[i]) - y[i];

   return 0; // Success
}

// Calculate:
// J_i1 = e^(x_2 * t_i)
// J_i2 = t_i x_1 e^(x_2 * t_i)
int eval_J(int n, int m, void * params, double const* x, double* J) {
   double x1 = x[0];
   double x2 = x[1];
   double const* t = ((struct params_type*) params)->t;

   for(int i=0; i<m; i++) {
      J[0*m+i] = exp(x2*t[i]);               // J_i1
      J[1*m+i] = t[i] * x1 * exp(x2*t[i]);   // J_i2
   }

   return 0; // Success
}


// Calculate:
// HF = sum_i r_i H_i
// Where H_i = [ 1                t_i e^(x_2 t_i)    ]
//             [ t_i e^(x_2 t_i)  t_i^2 e^(x_2 t_i)  ]
int eval_HF(int n, int m, void * params, double const* x, double const* r, double* HF) {
  //   double x1 = x[0];
   double x2 = x[1];
   double const* t = ((struct params_type*) params)->t;

   for(int i=0; i<n*n; i++) HF[i] = 0.0;
   for(int i=0; i<m; i++) {
      HF[    0] += r[i];                              // H_11
      HF[    1] += r[i] * t[i] * exp(x2*t[i]);        // H_21
      HF[1*n+1] += r[i] * t[i]*t[i] * exp(x2*t[i]);   // H_22
   }
   HF[1*n+0] = HF[1]; // H_12 by symmetry of Hessian

   return 0; // Success
}

int generic_test(int model, int method){
  // Data to be fitted
  int m = 5;
  char errstr[81];

  struct params_type params = {
    .t = (double []) { 1.0, 2.0, 4.0,  5.0,  8.0 },
    .y = (double []) { 3.0, 4.0, 6.0, 11.0, 20.0 }
  };
  
  // Initialize options values
  struct ral_nlls_options options;
  ral_nlls_default_options(&options);
  
  options.model = model;
  options.nlls_method = method;
  
  // Call fitting routine
  double x[2] = { 2.5, 0.25 }; // Initial guess
  struct ral_nlls_inform inform;
  nlls_solve(2, m, x, eval_r, eval_J, eval_HF, &params, &options, &inform,
	     NULL, NULL, NULL, NULL);
  if (model == 0) {
    printf("%s \n", inform.error_message);
    if (inform.status != -3){
      printf("ral_nlls() returned with error flag %d (expected -3)",inform.status);
      return -3;
    }
  }
  else{
    if(inform.status != 0) {
      printf("ral_nlls() returned with error flag %d\n", inform.status);
      return inform.status; // Error
    }
  }

  // If model is expected to pass,
  // call fitting routine with weights and bounds
  if (model > 0) { 
    x[0] =  2.5;
    x[1] =  0.25; // Reset Initial guess
    double weights[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double lower_bounds[2] = {0.0, 0.0};
    double upper_bounds[2] = {10.0, 10.0};
    
    nlls_solve(2, m, x, eval_r, eval_J, eval_HF, &params, &options, &inform,
	       weights, NULL, lower_bounds, upper_bounds);
    if(inform.status != 0) {
      printf("ral_nlls() returned with error flag %d\n", inform.status);
      return inform.status; // Error
    }
  }
    
  return 0; // Success!
}

int main(void){
  
  int no_errors = 0;
  int no_methods = 4;
  int status = 0;
  // passing tests....
  int model_array[4] = {1,2,3,0};
  for(int i=0; i<4;i++) { // loop over the methods
    for (int method=1;method<no_methods+1;method++){
      status = generic_test(model_array[i],method);
      if (status != 0) {
	status = 0;
	no_errors += 1;
	printf("Test failed, model = %d, method = %d\n",model_array[i],method);
      }
    }
  }
  
  if (no_errors > 0) {
    return no_errors;
  }
  printf("** C tests passed succcessfully! **\n");
  return 0; // success! 
}

