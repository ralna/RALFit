/* Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
 * Copyright (c) 2019, The Science and Technology Facilities Council (STFC)
 * All rights reserved.
 */

// examples/C/nlls_example.c
//
// Attempts to fit the model y_i = x_1 e^(x_2 t_i)
// For parameters x_1 and x_2, and input data (t_i, y_i)
#include "ral_nlls.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

struct params_type {
   double *t; // The m data points t_i
   double *y; // The m data points y_i
};

// Calculate r_i(x; t_i, y_i) = x_1 e^(x_2 * t_i) - y_i
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
// Where H_i = [ 0                t_i e^(x_2 t_i)        ]
//             [ t_i e^(x_2 t_i)  t_i^2 x_1 e^(x_2 t_i)  ]
int eval_HF(int n, int m, void * params, double const* x, double const* r, double* HF) {
   double x1 = x[0];
   double x2 = x[1];
   double const* t = ((struct params_type*) params)->t;

   for(int i=0; i<n*n; i++) HF[i] = 0.0;
   for(int i=0; i<m; i++) {
      HF[    0] += 0;                                      // H_11
      HF[    1] += r[i] * t[i] * exp(x2*t[i]);             // H_21
      HF[1*n+1] += r[i] * t[i]*t[i] * x1 * exp(x2*t[i]);   // H_22
   }
   HF[1*n+0] = HF[1]; // H_12 by symmetry of Hessian

   return 0; // Success
}

int main(void) {
   // Data to be fitted
   int m = 5;
   struct params_type params = {
      .t = (double []) { 1.0, 2.0, 4.0,  5.0,  8.0 },
      .y = (double []) { 3.0, 4.0, 6.0, 11.0, 20.0 }
   };

   // Initialize options values
   struct ral_nlls_options options;
   ral_nlls_default_options(&options);
   options.print_level = 2;
   options.check_derivatives = 1;
   options.print_options = true;
   options.derivative_test_tol = 2.0e-5;

   // Call fitting routine
   double x[2] = { 2.5, 0.25 }; // Initial guess
   struct ral_nlls_inform inform;
   printf("sending to nlls_solve\n");
   nlls_solve(2, m, x, eval_r, eval_J, eval_HF, &params, &options, &inform,
	      NULL, NULL, NULL, NULL);
   if(inform.status != 0) {
      printf("ral_nlls() returned with error flag %d\n", inform.status);
      return 1; // Error
   }

   // Print result
   printf ("Found a local optimum at x = %e %e\n", x[0], x[1]);
   printf ("Took %d iterations\n", inform.iter);
   printf ("     %d function evaluations\n", inform.f_eval);
   printf ("     %d gradient evaluations\n", inform.g_eval);
   printf ("     %d hessian evaluations\n", inform.h_eval);
}
