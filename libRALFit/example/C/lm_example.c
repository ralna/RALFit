/* Copyright (c) 2015, The Science and Technology Facilities Council (STFC)
 * All rights reserved.
 */
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include "ral_nlls.h"

#include "expfit.c"

#define N 40

void print_state (unsigned int iter, gsl_multifit_fdfsolver * s);

int
main (void)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int i, iter = 0;
  const int n = N;
  const int p = 3;

  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  double y[N], sigma[N];
  struct usertype d = { n, y, sigma};
  gsl_multifit_function_fdf f;
  double x_init[3] = { 1.0, 0.0, 0.0 };
  gsl_vector_view x = gsl_vector_view_array (x_init, p);
  const gsl_rng_type * type;
  gsl_rng * r;

  gsl_rng_env_setup();

  type = gsl_rng_default;
  r = gsl_rng_alloc (type);

  f.f = &expb_f;
  f.df = &expb_df;
  f.fdf = &expb_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

  /* This is the data to be fitted */

  for (i = 0; i < n; i++)
    {
      double t = i;
      y[i] = 1.0 + 5 * exp (-0.1 * t)
                 + gsl_ran_gaussian (r, 0.1);
      sigma[i] = 0.1;
      printf ("data: %u %g %g\n", i, y[i], sigma[i]);
    };

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);

  print_state (iter, s);

  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);

      printf ("status = %s\n", gsl_strerror (status));

      print_state (iter, s);

      if (status)
        break;

      status = gsl_multifit_test_delta (s->dx, s->x,
                                        1e-4, 1e-4);
    }
  while (status == GSL_CONTINUE && iter < 500);

  gsl_multifit_covar (s->J, 0.0, covar);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  {
    double chi = gsl_blas_dnrm2(s->f);
    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof));

    printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);

    printf ("A      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
    printf ("lambda = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
    printf ("b      = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
  }

  printf ("status = %s\n", gsl_strerror (status));
  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);
  gsl_rng_free (r);

  /*
    and herein lies the ral_nlls interface...
   */
  double x_ral[3] = {1.0, 0.0, 0.0 };
  struct nlls_control_type ral_options;
  struct nlls_inform_type ral_status;
  //  struct usertype ral_params;

  nlls_default_control(&ral_options);

  ral_options.print_level = 3;
  ral_options.maxit = 100;
  ral_options.model = 1;
  ral_options.nlls_method = 3;

  printf("p = %i n = %i", p, n);

  ral_nlls(p,n,
	   x_ral, eval_F, eval_J, eval_HF, &d,
	   &ral_status, &ral_options);

  printf(" A = %.5f \n lambda = %.5f \n b = %.5f\n ", x_ral[0], x_ral[1], x_ral[2]);

  return 0;
}

void
print_state (unsigned int iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
          "|f(x)| = %g\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
	  gsl_blas_dnrm2 (s->f));
}
