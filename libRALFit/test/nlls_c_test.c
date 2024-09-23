/* Copyright (c) 2015, The Science and Technology Facilities Council (STFC)
 * All rights reserved.
 * Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
 */
// test/nlls_c_test.c
//
// Test basic functionality of the c interface, plus any C-specific routines
#include "ral_nlls.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define VERBOSE 0

struct params_type {
  ral_real *t; // t_i
  ral_real *y; // y_i
};

// use the test example r_i(x; t_i, y_i) = x_1 e^(x_2 * t_i) - y_i

ral_int eval_r(ral_int n, ral_int m, void *params, ral_real const *x,
               ral_real *r) {
  ral_real x1 = x[0];
  ral_real x2 = x[1];
  ral_real const *t = ((struct params_type *)params)->t;
  ral_real const *y = ((struct params_type *)params)->y;

  for (ral_int i = 0; i < m; i++)
    r[i] = x1 * exp(x2 * t[i]) - y[i];

  return 0; // Success
}

// Calculate:
// J_i1 = e^(x_2 * t_i)
// J_i2 = t_i x_1 e^(x_2 * t_i)
ral_int eval_J(ral_int n, ral_int m, void *params, ral_real const *x,
               ral_real *J) {
  ral_real x1 = x[0];
  ral_real x2 = x[1];
  ral_real const *t = ((struct params_type *)params)->t;

  for (ral_int i = 0; i < m; i++) {
    J[0 * m + i] = exp(x2 * t[i]);             // J_i1
    J[1 * m + i] = t[i] * x1 * exp(x2 * t[i]); // J_i2
  }

  return 0; // Success
}

// Calculate:
// HF = sum_i r_i H_i
// Where H_i = [ 1                t_i e^(x_2 t_i)    ]
//             [ t_i e^(x_2 t_i)  t_i^2 e^(x_2 t_i)  ]
ral_int eval_HF(ral_int n, ral_int m, void *params, ral_real const *x,
                ral_real const *r, ral_real *HF) {
  //   ral_real x1 = x[0];
  ral_real x2 = x[1];
  ral_real const *t = ((struct params_type *)params)->t;

  for (ral_int i = 0; i < n * n; i++)
    HF[i] = 0.0;
  for (ral_int i = 0; i < m; i++) {
    HF[0] += r[i];                                        // H_11
    HF[1] += r[i] * t[i] * exp(x2 * t[i]);                // H_21
    HF[1 * n + 1] += r[i] * t[i] * t[i] * exp(x2 * t[i]); // H_22
  }
  HF[1 * n + 0] = HF[1]; // H_12 by symmetry of Hessian

  return 0; // Success
}

ral_int generic_test(ral_int model, ral_int method) {
  // Data to be fitted
  ral_int m = 5;

  struct params_type params = {.t = (ral_real[]){1.0, 2.0, 4.0, 5.0, 8.0},
                               .y = (ral_real[]){3.0, 4.0, 6.0, 11.0, 20.0}};

  // Initialize options values
  struct ral_nlls_options options;
  ral_nlls_default_options(&options);

  options.model = model;
  options.nlls_method = method;
  options.print_level = VERBOSE;
  options.print_options = VERBOSE > 0;
#ifdef SINGLE_PRECISION
  // Relax and tweak
  options.stop_s = 1.0e-6;
#endif

  // Call fitting routine
  ral_real x[2] = {2.5, 0.25}; // Initial guess
  struct ral_nlls_inform inform;
  nlls_solve(2, m, x, eval_r, eval_J, eval_HF, &params, &options, &inform, NULL,
             NULL, NULL, NULL);
  if (model == 0) {
    printf("%s \n", inform.error_message);
    if (inform.status != -3) {
      printf("nlls_solve() returned with error flag %d (expected -3)",
             inform.status);
      return -3;
    }
  } else {
    if (inform.status != 0) {
      printf("nlls_solve() returned with error flag %d\n", inform.status);
      return inform.status; // Error
    }
  }

  // If model is expected to pass,
  // call fitting routine with weights and bounds
  if (model > 0) {
    x[0] = 2.5;
    x[1] = 0.25; // Reset Initial guess
    ral_real weights[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
    ral_real lower_bounds[2] = {0.0, 0.0};
    ral_real upper_bounds[2] = {10.0, 10.0};

    nlls_solve(2, m, x, eval_r, eval_J, eval_HF, &params, &options, &inform,
               weights, NULL, lower_bounds, upper_bounds);
    if (inform.status != 0) {
      printf("nlls_solve() returned with flag %d\n", inform.status);
      return inform.status; // Error
    }
  }

  return 0; // Success!
}

ral_int main(void) {

  ral_int no_errors = 0;
  ral_int no_methods = 4;
  ral_int status = 0;
  ral_int cnt = 0;
  // passing tests....
  ral_int model_array[4] = {1, 2, 3, 0};
  for (ral_int i = 0; i < 4; i++) { // loop over the methods
    for (ral_int method = 1; method < no_methods + 1; method++) {
      ++cnt;
      printf("\n [Test: #%i.......] config(model=%i, method=%i)\n", cnt,
             model_array[i], method);
      status = generic_test(model_array[i], method);
      if (status != 0) {
        status = 0;
        no_errors += 1;
        printf(" [Test: #%i...FAIL] config(model=%i, method=%i)\n", cnt,
               model_array[i], method);
      } else {
        printf(" [Test: #%i...PASS] config(model=%i, method=%i)\n", cnt,
               model_array[i], method);
      }
    }
  }

  if (no_errors > 0) {
    printf("\n\n!! C tests FAILED: %i / %i\n", no_errors, cnt);
    return no_errors;
  }
  printf("\n\n** C tests passed succcessfully! **\n");
  return 0; // success!
}
