/* Copyright (c) 2015, The Science and Technology Facilities Council (STFC)
 * All rights reserved.
 * Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
 */
#include "ral_nlls.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>*/

/* define the usertype */
struct usertype {
  ral_real y_data[67];
  ral_real x_data[67];
};

void generate_data_example(ral_real *x_data, ral_real *y_data,
                           const ral_int m); // prototype
ral_int eval_F(const ral_int n, const ral_int m, void *params,
               const ral_real *X, ral_real *f);
ral_int eval_J(const ral_int n, const ral_int m, void *params,
               const ral_real *X, ral_real *f);
ral_int eval_HF(const ral_int n, const ral_int m, void *params,
                const ral_real *X, const ral_real *f, ral_real *hf);

/* A c driver for the ral_nlls program */
ral_int main(void) {

  /* Problem data */
  const ral_int n = 2;
  const ral_int m = 67;

  /* Derived types */
  struct ral_nlls_options options;
  struct ral_nlls_inform status;
  struct usertype params;

  printf("===============\n");
  printf("RAL NLLS driver\n");
  printf("~  C version  ~\n");
  printf("===============\n");

  /* Generate the data... */
  generate_data_example(params.x_data, params.y_data, m);

  ral_real x[n];
  x[0] = 1.0;
  x[1] = 2.0;
  ral_real x_exp[n];
  x_exp[0] = 0.32;
  x_exp[1] = 0.0275;

  ral_real tol = 1.0e-4;

  ral_nlls_default_options(&options);

  options.print_level = 3;
  options.print_options = true;
  options.check_derivatives = 1;

#ifdef SINGLE_PRECISION
  options.fd_step = 1.0e-3;
  options.derivative_test_tol = 3.0e-3;
#endif

  nlls_solve(n, m, x, eval_F, eval_J, NULL, &params, &options, &status, NULL,
             NULL, NULL, NULL);

  bool ok = status.status == 0;

  if (!ok)
    printf(" ** 1st call to solver did not exit successfully.\n");

  printf("\nSolution:\n");
  for (ral_int i = 0; i < n; i++) {
    ral_real err = fabs(x[i] - x_exp[i]);
    bool oki = err <= tol;
    printf("  %5.4f (%5.4f) %9.2e  %s\n", x[i], x_exp[i], err,
           (oki ? "Ok" : "FAIL!"));
    ok &= oki;
  }

  ral_int exit_status = !ok;

  options.print_options = false;
  options.check_derivatives = 1; // Keep at 1
  options.print_level = 0;
  x[0] = 1.0;
  x[1] = 2.0;
  nlls_solve(n, m, x, eval_F, NULL, eval_HF, &params, &options, &status, NULL,
             NULL, NULL, NULL);

  ok = status.status == 0;
  exit_status += !ok;

  if (!ok)
    printf(" ** 2nd call to solver (FD) did not exit successfully.\n");

  printf("\nSolution:\n");
  for (ral_int i = 0; i < n; i++) {
    ral_real err = fabs(x[i] - x_exp[i]);
    bool oki = err <= tol;
    printf("  %5.4f (%5.4f) %9.2e  %s\n", x[i], x_exp[i], err,
           (oki ? "Ok" : "FAIL!"));
    ok &= oki;
  }

  options.check_derivatives = 0; // Now set to 0
  x[0] = 1.0;
  x[1] = 2.0;
  nlls_solve(n, m, x, eval_F, NULL, eval_HF, &params, &options, &status, NULL,
             NULL, NULL, NULL);

  ok = status.status == 0;

  if (!ok)
    printf(" ** 3rd call to solver (FD) did not exit successfully.\n");

  printf("\nSolution:\n");
  for (ral_int i = 0; i < n; i++) {
    ral_real err = fabs(x[i] - x_exp[i]);
    bool oki = err <= tol;
    printf("  %5.4f (%5.4f) %9.2e  %s\n", x[i], x_exp[i], err,
           (oki ? "Ok" : "FAIL!"));
    ok &= oki;
  }

  exit_status += !ok;

  return exit_status;
}

/* Do a function evaluation */
ral_int eval_F(ral_int n, ral_int m, void *params, const ral_real *X,
               ral_real *f) {

  struct usertype *myparams = (struct usertype *)params;

  ral_int i;

  for (i = 0; i < m; i++) {
    f[i] = myparams->y_data[i] - exp(X[0] * myparams->x_data[i] + X[1]);
  }

  return 0;
}

/* Evaluate the Jacobian */
ral_int eval_J(const ral_int n, const ral_int m, void *params,
               const ral_real *X, ral_real *J) {

  struct usertype *myparams = (struct usertype *)params;

  ral_int i;

  for (i = 0; i < m; i++) {
    J[i] = -myparams->x_data[i] * exp(X[0] * myparams->x_data[i] + X[1]);
    J[m + i] = -exp(X[0] * myparams->x_data[i] + X[1]);
  }

  return 0;
}

/* Evaluate the Hessian */
ral_int eval_HF(const ral_int n, const ral_int m, void *params,
                const ral_real *X, const ral_real *f, ral_real *hf) {

  ral_int i;

  for (i = 0; i < n * n; i++) {
    hf[i] = 0.0;
  }

  return 0;
}

/* Generate some example data... */
void generate_data_example(ral_real *x_data, ral_real *y_data,
                           const ral_int m) {

  ral_int i;
  /* Note the 67 here needs to be hard-coded, and you can't initialize
     an array with a variable length in C */
  ral_real tempx[67] = {0.0,
                        0.075000000000000,
                        0.150000000000000,
                        0.225000000000000,
                        0.300000000000000,
                        0.375000000000000,
                        0.450000000000000,
                        0.525000000000000,
                        0.600000000000000,
                        0.675000000000000,
                        0.750000000000000,
                        0.825000000000000,
                        0.900000000000000,
                        0.975000000000000,
                        1.050000000000000,
                        1.125000000000000,
                        1.200000000000000,
                        1.275000000000000,
                        1.350000000000000,
                        1.425000000000000,
                        1.500000000000000,
                        1.575000000000000,
                        1.650000000000000,
                        1.725000000000000,
                        1.800000000000000,
                        1.875000000000000,
                        1.950000000000000,
                        2.025000000000000,
                        2.100000000000000,
                        2.175000000000000,
                        2.250000000000000,
                        2.325000000000000,
                        2.400000000000000,
                        2.475000000000000,
                        2.550000000000000,
                        2.625000000000000,
                        2.700000000000000,
                        2.775000000000000,
                        2.850000000000000,
                        2.925000000000000,
                        3.000000000000000,
                        3.075000000000000,
                        3.150000000000000,
                        3.225000000000001,
                        3.300000000000000,
                        3.375000000000000,
                        3.450000000000000,
                        3.525000000000000,
                        3.600000000000001,
                        3.675000000000000,
                        3.750000000000000,
                        3.825000000000000,
                        3.900000000000000,
                        3.975000000000000,
                        4.050000000000001,
                        4.125000000000000,
                        4.200000000000000,
                        4.275000000000000,
                        4.350000000000001,
                        4.425000000000000,
                        4.500000000000000,
                        4.575000000000000,
                        4.650000000000000,
                        4.725000000000001,
                        4.800000000000000,
                        4.875000000000000,
                        4.950000000000000};

  ral_real tempy[67] = {0.907946872110432, 1.199579396036134, 1.060092431384317,
                        1.298370500472354, 0.952768858414788, 1.209665290655204,
                        1.256912538155493, 1.163922146095987, 1.004877938808100,
                        1.205944250961060, 0.952693297695969, 1.449662692280761,
                        1.402015259144406, 1.378094012325746, 1.560882147577552,
                        1.437185539058121, 1.559853079888265, 1.877814947316832,
                        1.818781749024682, 1.375546045112591, 1.233967904388409,
                        1.887793124397751, 1.610237096463521, 1.787032484792262,
                        1.850015127982676, 2.120553361509177, 1.942913663511919,
                        2.106517132599766, 2.271787117356578, 1.727554346001754,
                        2.002909500898113, 1.975837413903495, 2.337446525801909,
                        1.960190841677278, 2.447097025572309, 2.161663720225506,
                        2.748798529374621, 2.507814238594416, 2.423769408403069,
                        2.578119353028746, 2.460310096221557, 2.638362783992324,
                        2.765540456237868, 2.837165966564409, 3.179711963042789,
                        3.245315453091675, 3.289631922410174, 3.360995198615834,
                        3.470489725998371, 3.169513520153466, 3.363740517933189,
                        3.665288099084969, 3.620334359722351, 4.018911445550667,
                        3.512715166706162, 3.874661411575566, 4.197746303653517,
                        3.703511523106007, 4.076351488309604, 4.056340365649961,
                        4.297751562451419, 4.373076571153739, 4.577093065941748,
                        4.856619059058190, 4.927350280596274, 4.703122139742729,
                        4.870205182453842};

  for (i = 0; i < m; i++) {
    x_data[i] = tempx[i];
    y_data[i] = tempy[i];
  }
}
