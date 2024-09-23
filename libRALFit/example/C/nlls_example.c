/* Copyright (c) 2019, The Science and Technology Facilities Council (STFC)
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
   ral_real *t; // The m data points t_i
   ral_real *y; // The m data points y_i
};

// Calculate r_i(x; t_i, y_i) = x_1 e^(x_2 * t_i) - y_i
ral_int eval_r(ral_int n, ral_int m, void * params, ral_real const* x, ral_real* r) {
   ral_real x1 = x[0];
   ral_real x2 = x[1];
   ral_real const* t = ((struct params_type*) params)->t;
   ral_real const* y = ((struct params_type*) params)->y;

   for(ral_int i=0; i<m; i++)
      r[i] = x1 * exp(x2*t[i]) - y[i];

   return 0; // Success
}

// Calculate:
// J_i1 = e^(x_2 * t_i)
// J_i2 = t_i x_1 e^(x_2 * t_i)
ral_int eval_J(ral_int n, ral_int m, void * params, ral_real const* x, ral_real* J) {
   ral_real x1 = x[0];
   ral_real x2 = x[1];
   ral_real const* t = ((struct params_type*) params)->t;

   for(ral_int i=0; i<m; i++) {
      J[0*m+i] = exp(x2*t[i]);               // J_i1
      J[1*m+i] = t[i] * x1 * exp(x2*t[i]);   // J_i2
   }

   return 0; // Success
}

// Calculate:
// HF = sum_i r_i H_i
// Where H_i = [ 0                t_i e^(x_2 t_i)        ]
//             [ t_i e^(x_2 t_i)  t_i^2 x_1 e^(x_2 t_i)  ]
ral_int eval_HF(ral_int n, ral_int m, void * params, ral_real const* x, ral_real const* r, ral_real* HF) {
   ral_real x1 = x[0];
   ral_real x2 = x[1];
   ral_real const* t = ((struct params_type*) params)->t;

   for(ral_int i=0; i<n*n; i++) HF[i] = 0.0;
   for(ral_int i=0; i<m; i++) {
      HF[    0] += 0;                                      // H_11
      HF[    1] += r[i] * t[i] * exp(x2*t[i]);             // H_21
      HF[1*n+1] += r[i] * t[i]*t[i] * x1 * exp(x2*t[i]);   // H_22
   }
   HF[1*n+0] = HF[1]; // H_12 by symmetry of Hessian

   return 0; // Success
}

int main(void) {
   // Data to be fitted
   ral_int m = 5;
   struct params_type params = {
      .t = (ral_real []) { 1.0, 2.0, 4.0,  5.0,  8.0 },
      .y = (ral_real []) { 3.0, 4.0, 6.0, 11.0, 20.0 }
   };

   // Initialize options values
   struct ral_nlls_options options;
   ral_nlls_default_options(&options);
   options.print_level = 2;
   options.check_derivatives = 1;
   options.print_options = true;
   options.derivative_test_tol = 5.0e-4;
   options.fd_step = 1.0e-4;

   // Call fitting routine
   ral_real x[2] = { 2.5, 0.25 }; // Initial guess
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
