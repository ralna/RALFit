/*
 * Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

/* Fit the convolution model y_i = [Al * lognorma(a, b)]_i + [Ag * normal(mu,
 * sigma)]_i given the density observations at the measured diameter sizes.
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "ral_nlls.h"

#define CHK(OK) ((OK) ? "PASS" : "FAIL")
#ifdef SINGLE_PRECION
#define LOG logf
#define FABS fabsf
#define SQRT sqrtf
#define POW powf
#define EXP expf
#else
#define LOG log
#define FABS fabs
#define SQRT sqrt
#define POW pow
#define EXP exp
#endif

struct params_t {
  const int *diameter;
  const ral_real *density;
};

const ral_real pi = 3.14159265358979323846;

// Empirical data (histogram)
// Track diameter
int diameter[64] = {
    1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16,
    17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
    33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
    49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64};
// Diameter density count
ral_real density[64] = {
    0.0722713864, 0.0575221239, 0.0604719764, 0.0405604720, 0.0317109145,
    0.0309734513, 0.0258112094, 0.0228613569, 0.0213864307, 0.0213864307,
    0.0147492625, 0.0213864307, 0.0243362832, 0.0169616519, 0.0095870206,
    0.0147492625, 0.0140117994, 0.0132743363, 0.0147492625, 0.0140117994,
    0.0140117994, 0.0132743363, 0.0117994100, 0.0132743363, 0.0110619469,
    0.0103244838, 0.0117994100, 0.0117994100, 0.0147492625, 0.0110619469,
    0.0132743363, 0.0206489676, 0.0169616519, 0.0169616519, 0.0280235988,
    0.0221238938, 0.0235988201, 0.0221238938, 0.0206489676, 0.0228613569,
    0.0184365782, 0.0176991150, 0.0132743363, 0.0132743363, 0.0088495575,
    0.0095870206, 0.0073746313, 0.0110619469, 0.0036873156, 0.0051622419,
    0.0058997050, 0.0014749263, 0.0022123894, 0.0029498525, 0.0014749263,
    0.0007374631, 0.0014749263, 0.0014749263, 0.0007374631, 0.0000000000,
    0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000};

const struct params_t params = {diameter, density};

// scaled Log-Normal density distribution Al amplitude * Log-Normal(a, b)
ral_real lognormal(ral_real d, ral_real a, ral_real b, ral_real Al) {
  return Al / (d * b * SQRT(2.0 * pi)) *
         EXP(-(POW(LOG(d) - a, 2.0)) / (2.0 * POW(b, 2.0)));
}

// scaled normal density distribution Ag amplitude * Normal(mu, sigma)
ral_real gaussian(ral_real d, ral_real mu, ral_real sigma, ral_real Ag) {
  return Ag * EXP(-0.5 * POW((d - mu) / sigma, 2)) / (sigma * SQRT(2.0 * pi));
}

// residuals for the convolution model
ral_int eval_r(ral_int n, ral_int m, void *params, ral_real const *x, ral_real *r) {
  ral_real const a = x[0];
  ral_real const b = x[1];
  ral_real const Al = x[2];
  ral_real const mu = x[3];
  ral_real const sigma = x[4];
  ral_real const Ag = x[5];
  int const *d = ((struct params_t *)params)->diameter;
  ral_real const *y = ((struct params_t *)params)->density;
  for (int i = 0; i < m; ++i)
    r[i] = lognormal(d[i], a, b, Al) + gaussian(d[i], mu, sigma, Ag) - y[i];
  return 0;
}

// residual Jacobian for the convolution model
ral_int eval_j(ral_int n, ral_int m, void *params, ral_real const *x, ral_real *J) {
  ral_real const a = x[0];
  ral_real const b = x[1];
  ral_real const Al = x[2];
  ral_real const mu = x[3];
  ral_real const sigma = x[4];
  ral_real const Ag = x[5];
  int const *d = ((struct params_t *)params)->diameter;
  for (int i = 0; i < m; ++i) {
    ral_real l = lognormal(d[i], a, b, Al);
    J[i * n + 0] = (LOG(d[i]) - a) / POW(b, 2.0) * l;
    J[i * n + 1] = (POW(LOG(d[i]) - a, 2.0) - POW(b, 2)) / POW(b, 3) * l;
    J[i * n + 2] = lognormal(d[i], a, b, 1.0);
    ral_real g = gaussian(d[i], mu, sigma, Ag);
    J[i * n + 3] = (d[i] - mu) / POW(sigma, 2.0) * g;
    J[i * n + 4] =
        (POW(d[i] - mu, 2.0) - POW(sigma, 2.0)) / POW(sigma, 3.0) * g;
    J[i * n + 5] = gaussian(d[i], mu, sigma, 1.0);
  }
  return 0; // Success
}

int main(void) {
  printf("\nNuclear tracks: a nonlinear least-squares example (using finite "
         "differences)\n");

  const ral_int n = 6; /* vector (a, b, Al, mu, sigma, Ag) */
  const ral_int m = 64;
  ral_real SP1[6] = {1.65, 0.9, 1.0, 30.0, 1.5, 0.25};
  ral_real x[6];

  ral_real lower_bounds[m];
  ral_real weights[m];
  for (int j = 0; j < 55; ++j)
    weights[j] = 1.0;
  for (int j = 55; j <= 63; ++j)
    weights[j] = 5.0;
  ral_real wsum = 1.0 * (m - (63 - 55)) + 5.0 * (63 - 55);
  // normalize weights and add lower bounds
  for (int j = 0; j < 64; ++j) {
    lower_bounds[j] = 0.0;
    weights[j] /= wsum;
  }

  // Initialize options values
  struct ral_nlls_options options;
  ral_nlls_default_options(&options);
  options.print_level = 2;
  options.check_derivatives = 1;
  options.Fortran_Jacobian = false;
#ifdef SINGLE_PRECISION
  options.fd_step = 1.0e-3;
  options.derivative_test_tol = 5.0e-2;
  const ral_real x_exp[6] = {1.97, 1.36, 0.68, 36.98, 6.96, 0.34};
#else
  options.fd_step = 1.0e-8;
  const ral_real x_exp[6] = {1.99, 1.37, 0.68, 36.64, 7.08, 0.34};
#endif

  // initialize the workspace
  void *workspace;
  void *inner_workspace;

  ral_real tol = 0.0;
  bool ok = true, oki, okj;
  // init workspace allocates and links together workspace with inner_workspace
  ral_nlls_init_workspace(&workspace, &inner_workspace);

  struct ral_nlls_inform inform;

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < n; ++j)
      x[j] = SP1[j];

    if (i) {
      options.print_level = 1;
      printf("\n\n Solving using exact jacobian + check of derivatives\n\n");
      nlls_solve(n, m, x, eval_r, eval_j, NULL, &params, &options, &inform,
                 weights, NULL, lower_bounds, NULL);
      tol = 1.0e-2;

    } else {
      printf("\n\n Solving using finite-differences\n\n");
      nlls_solve(n, m, x, eval_r, NULL, NULL, &params, &options, &inform,
                 weights, NULL, lower_bounds, NULL);
      // relax testing tolerance for lower precision case
      tol = 5.0e-2;
    }

    printf("Coefficients:\n");

    oki = FABS(x[0] - x_exp[0]) <= tol;
    ok &= oki;
    okj = FABS(x[1] - x_exp[1]) <= tol;
    ok &= okj;
    printf("x[0] = %0.2f (%0.2f)  %s   x[1] = %0.2f (%2.2f) %s\n", x[0],
           x_exp[0], CHK(oki), x[1], x_exp[1], CHK(okj));
    oki = FABS(x[2] - x_exp[2]) <= tol;
    ok &= oki;
    okj = FABS(x[3] - x_exp[3]) <= tol;
    ok &= okj;
    printf("x[2] = %0.2f (%0.2f)  %s   x[3] = %0.2f (%2.2f) %s\n", x[2],
           x_exp[2], CHK(oki), x[3], x_exp[3], CHK(okj));
    oki = FABS(x[4] - x_exp[4]) <= tol;
    ok &= oki;
    okj = FABS(x[5] - x_exp[5]) <= tol;
    ok &= okj;
    printf("x[4] = %0.2f (%0.2f)  %s   x[5] = %0.2f (%2.2f) %s\n", x[4],
           x_exp[4], CHK(oki), x[5], x_exp[5], CHK(okj));

    printf("\n");
    ok &= inform.status == 0;
    if (ok) {
      printf("Regression computed successfully!\n");
      printf("Fit error                      : %e\n", inform.obj);
      printf("Norm of residual gradient      : %e (%e)\n", inform.norm_g,
             inform.scaled_g);
      printf("Objective fun calls (fin diff) : %i (%i)\n", inform.f_eval,
             inform.fd_f_eval);
      printf("Objective fun calls            : %i\n", inform.g_eval);
    } else {
      printf("Something wrong happened during the fit.\n");
      printf("%s\n", inform.error_message);
    }
  }

  ral_nlls_free_workspace(&workspace);
  ral_nlls_free_workspace(&inner_workspace);
  return ok ? 0 : 6;
}
