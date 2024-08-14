/* Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
 * Copyright (c) 2019, The Numerical Algorithms Group Ltd (NAG)
 * All rights reserved.
 * Copyright (c) 2015, The Science and Technology Facilities Council (STFC)
 * All rights reserved.
 */

// examples/Fortran/nlls_example2.f90
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

   // initialize the workspace
   void * workspace;
   void * inner_workspace;

   // init workspace allocates and links together workspace with inner_workspace
   ral_nlls_init_workspace(&workspace, &inner_workspace);

   // Call fitting routine
   double x[2] = { 1.0, 1.0 }; // Initial guess
   double x_exp[2] = {6.914573e-3, 1.0};

   double lower_bounds[2] = { 0.0, 1.0};
   double upper_bounds[2] = { 1.0, 10.0};

   struct ral_nlls_inform inform;

   nlls_solve(2, m, x, eval_r, eval_J, NULL, &params,
	      &options, &inform, NULL, NULL, lower_bounds, upper_bounds);

   ral_nlls_free_workspace(&workspace);
   ral_nlls_free_workspace(&inner_workspace);

  int ok = 0;
  if (inform.status != 0) {
    printf("Status = %i [%s]\n", inform.status, inform.error_message);
  } else {
    // Print result
    double tol = 1.0e-6;
    char ok0 = (fabs(x[0]-x_exp[0]) <= tol) ? ' ' : 'X';
    char ok1 = (fabs(x[1]-x_exp[1]) <= tol) ? ' ' : 'X';
    ok = ok0 == ' ' && ok1 == ' ';
    printf("Found a local optimum at x = %e, %e (expected: %e %c, %e %c)\n",
           x[0], x[1], x_exp[0], ok0, x_exp[1], ok1);
    printf("Took %d iterations\n", inform.iter);
    printf("     %d function evaluations\n", inform.f_eval);
    printf("     %d gradient evaluations\n", inform.g_eval);
  }
  return ok ? 0 : 6;
}
