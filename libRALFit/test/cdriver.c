#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/*#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>*/
#include "ral_nlls.h"

/* define the usertype */
struct usertype {
  double y_data[67];
  double x_data[67];
};

void generate_data_example ( double *x_data, double *y_data, const int m); // prototype
int eval_F  ( const int n, const int m,  void *params, 
      const double *X, double *f);
int eval_J  ( const int n, const int m,  void *params, 
      const double *X, double *f);
int eval_HF ( const int n, const int m,  void *params, 
	       const double *X, const double *f, double *hf);

/* A c driver for the ral_nlls program */
int main(void) {
  
  /* Problem data */
  const int n = 2;
  const int m = 67;
  
  /* Derived types */
  struct ral_nlls_options options;
  struct ral_nlls_inform status;
  struct usertype params;
  
  printf("===============\n");
  printf("RAL NLLS driver\n");
  printf("~  C version  ~\n");
  printf("===============\n");

  /* Generate the data... */
  generate_data_example(params.x_data,params.y_data,m);

  double x[n];
  x[0] = 1.0;
  x[1] = 2.0;

  ral_nlls_default_options(&options);
  
  options.print_level = 3;
  
  nlls_solve(n, m, x,
	   eval_F, eval_J, eval_HF, &params,
	     &options, &status, NULL, NULL, NULL, NULL );

  int i;
  printf("\nX = \n");
  for(i=0; i < n; i++) {
    printf("  %5.4f \n",x[i]);
  }

  return 0; /* success */
}

/* Do a function evaluation */
int eval_F(int n, int m, void *params, 
	      const double *X, double *f){
  
  struct usertype *myparams = (struct usertype *) params;

  int i;
  
  for(i=0; i<m; i++) {
    f[i] = myparams->y_data[i] - exp( X[0] * myparams->x_data[i] + X[1] );
  }

  return 0;
}

/* Evaluate the Jacobian */
int eval_J( const int n, const int m, void *params, 
	     const double *X, double *J){
  
  struct usertype *myparams = (struct usertype *) params;

  int i;

  for(i=0; i<m; i++) {
    J[i] = -myparams->x_data[i] * exp( X[0] * myparams->x_data[i] + X[1] );
    J[m + i] = - exp( X[0] * myparams->x_data[i] + X[1] );
  }
  
  return 0;
}

/* Evaluate the Hessian */
int eval_HF( const int n, const int m, void *params, 
	     const double *X, const double *f, 
	     double *hf){
  
  struct usertype *myparams = (struct usertype *) params;

  int i;

  for(i=0; i<n*n; i++) {
    hf[i] = 0.0;
  }
  
  return 0;
}

/* Generate some example data... */
void generate_data_example( double *x_data, double *y_data, const int m ) {
  
  int i;
  /* Note the 67 here needs to be hard-coded, and you can't initialize 
     an array with a variable length in C */
  double tempx[67] =  { 0.0, 
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
  
  double tempy[67] =  {0.907946872110432, 
		       1.199579396036134, 
		       1.060092431384317, 
		       1.298370500472354, 
		       0.952768858414788, 
		       1.209665290655204, 
		       1.256912538155493, 
		       1.163922146095987, 
		       1.004877938808100, 
		       1.205944250961060, 
		       0.952693297695969, 
		       1.449662692280761, 
		       1.402015259144406, 
		       1.378094012325746, 
		       1.560882147577552, 
		       1.437185539058121, 
		       1.559853079888265, 
		       1.877814947316832, 
		       1.818781749024682, 
		       1.375546045112591, 
		       1.233967904388409, 
		       1.887793124397751, 
		       1.610237096463521, 
		       1.787032484792262, 
		       1.850015127982676, 
		       2.120553361509177, 
		       1.942913663511919, 
		       2.106517132599766, 
		       2.271787117356578, 
		       1.727554346001754, 
		       2.002909500898113, 
		       1.975837413903495, 
		       2.337446525801909, 
		       1.960190841677278, 
		       2.447097025572309, 
		       2.161663720225506, 
		       2.748798529374621, 
		       2.507814238594416, 
		       2.423769408403069, 
		       2.578119353028746, 
		       2.460310096221557, 
		       2.638362783992324, 
		       2.765540456237868, 
		       2.837165966564409, 
		       3.179711963042789, 
		       3.245315453091675, 
		       3.289631922410174, 
		       3.360995198615834, 
		       3.470489725998371, 
		       3.169513520153466, 
		       3.363740517933189, 
		       3.665288099084969, 
		       3.620334359722351, 
		       4.018911445550667, 
		       3.512715166706162, 
		       3.874661411575566, 
		       4.197746303653517, 
		       3.703511523106007, 
		       4.076351488309604, 
		       4.056340365649961, 
		       4.297751562451419, 
		       4.373076571153739, 
		       4.577093065941748, 
		       4.856619059058190, 
		       4.927350280596274, 
		       4.703122139742729, 
		       4.870205182453842};
  
  for(i=0;i<m;i++){
    x_data[i] = tempx[i];
    y_data[i] = tempy[i];
  }

}
