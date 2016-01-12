/* expfit.c -- model functions for exponential + background */

#include <math.h>
#include <stddef.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

struct usertype {
  size_t n;
  double * y;
  double * sigma;
};

int expb_f (const gsl_vector * x, void *data, 
        gsl_vector * f) {
  size_t n = ((struct usertype *)data)->n;
  double *y = ((struct usertype *)data)->y;
  double *sigma = ((struct usertype *) data)->sigma;

  double A = gsl_vector_get (x, 0);
  double lambda = gsl_vector_get (x, 1);
  double b = gsl_vector_get (x, 2);

  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Model Yi = A * exp(-lambda * i) + b */
      double t = i;
      double Yi = A * exp (-lambda * t) + b;
      gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
    }

  return GSL_SUCCESS;
}

int expb_df (const gsl_vector * x, void *params, 
         gsl_matrix * J) {

  struct usertype *myparams = (struct usertype *) params;
  
  int n = myparams->n;
  double *sigma = myparams->sigma;


  double A = gsl_vector_get (x, 0);
  double lambda = gsl_vector_get (x, 1);

  size_t i;
  
  /*for (i = 0; i<10; i++){
    printf("sigma[%z] = %5.3f\n",i,sigma[i]);
    }*/
  
  for (i = 0; i < n; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/sigma[i],      */
      /*       Yi = A * exp(-lambda * i) + b  */
      /* and the xj are the parameters (A,lambda,b) */
      double t = i;
      double s = sigma[i];
      double e = exp(-lambda * t);
      gsl_matrix_set (J, i, 0, e/s); 
      gsl_matrix_set (J, i, 1, -t * A * e/s);
      gsl_matrix_set (J, i, 2, 1/s);
    }
  return GSL_SUCCESS;
}

int expb_fdf (const gsl_vector * x, void *data,
          gsl_vector * f, gsl_matrix * J) {

  expb_f (x, data, f);
  expb_df (x, data, J);

  return GSL_SUCCESS;
}

void eval_F ( int fstatus, const int n, const int m, 
	      const double *x, double *f, const void *params){
  
  // first, convert x from an array into a gsl_vector...  
  gsl_vector * x_gsl = gsl_vector_alloc(n);
  int i;
  for (i = 0; i < n; i++) {
    gsl_vector_set(x_gsl,i,x[i]);
  }

  // then call the expb_f function
  gsl_vector * f_gsl =  gsl_vector_alloc(m);

  expb_f ( x_gsl,  params,  f_gsl);
  
  // then convert f back into an array...
  //  f = f_gsl->data;
  for (i = 0; i < m; i++){
    f[i] = gsl_vector_get( f_gsl,i);
    //    f[i] = f_gsl->data[i];
  }
  
  gsl_vector_free(x_gsl);
  gsl_vector_free(f_gsl);

}


void eval_J  ( int fstatus, const int n, const int m, 
	       const double *x, double *J, const void *params){

  // first, convert x from an array into a gsl_vector...
  gsl_vector * x_gsl = gsl_vector_alloc(n);
  int i;
  int jj;
  for (i = 0; i < n; i++) {
    gsl_vector_set(x_gsl,i,x[i]);
  }
  
  // then call the expb_f function
  gsl_matrix * J_gsl = gsl_matrix_alloc(40,3);
  expb_df (x_gsl, params, J_gsl);
  
  // then convert J into an array...
  // (find better way...)
  gsl_matrix * J_gsl_t = gsl_matrix_alloc(3,40);
  gsl_matrix_transpose_memcpy(J_gsl_t, J_gsl);
  for ( i = 0; i < m*n; i++){
    J[i] = J_gsl_t->data[i];
  }
  // the following code loses the pointer to J (maybe...)
  // J = J_gsl_t->data;
  gsl_vector_free(x_gsl);
  gsl_matrix_free(J_gsl);
  gsl_matrix_free(J_gsl_t);
  
  
  
}

void eval_HF ( int fstatus, const int n, const int m, 
	       const double *x, const double *f, double *hf, const void *params){
  
}
