/* COPYRIGHT (c) 2015 Science and Technology Facilities Council (STFC)
 * All rights reserved
 */
#ifdef __cplusplus
extern "C" {
#endif

#ifndef ral_nlls_h
#define ral_nlls_h

#ifndef nlls_control_type
#define nlls_default_control nlls_default_control_d
#define nlls_control_type nlls_control_type_d
#define nlls_inform_type nlls_inform_type_d
#define ral_nlls ral_nlls_d
#define ral_nlls_int_func ral_nlls_int_func_d
#define c_test_pass_f c_test_pass_f_d
#endif

typedef double nllspkgtype_d_;

/* Derived type to hold control parameters for ral_nlls */
struct nlls_control_type_d {
  int f_arrays; /* Use 1-based indexing if true(!=0) else 0-based */
  
  int error; /* Fortran output stream for error messages */
  int out;   /* Fortran output stream for general messages */
  int print_level; /* levels of print output */
  int maxit; /* maximum number of iterations */
  int model; /* what model to use? */
  int nlls_method; /* what nlls method to use? */
  int lls_solver; /* which lls solver to use? */
  nllspkgtype_d_ stop_g_absolute; /* absolute stopping tolerance */
  nllspkgtype_d_ stop_g_relative; /* relative stopping tolerance */
  nllspkgtype_d_ initial_radius; /* initial trust region radius */
  nllspkgtype_d_ maximum_radius; /* maximux trust region radius */
  nllspkgtype_d_ eta_successful; /* trust region step successful level */
  nllspkgtype_d_ eta_very_successful; /* trust region step very successful */
  nllspkgtype_d_ eta_too_successful; /* trust region step too successful */
  nllspkgtype_d_ radius_increase; /* how much to increase the radius by? */
  nllspkgtype_d_ radius_reduce; /* how much to reduce the radius by? */
  nllspkgtype_d_ radius_reduce_max; /* max amount to reduce the radius by */
};

struct nlls_inform_type_d {
  int status; /* flag */
};

/* Set default values of control */
void nlls_default_control_d( struct nlls_control_type *options );

/* define the eval_f_type */
typedef void (*eval_f_type) (int fstatus, 
			     int n, 
			     int m,
			     const nllspkgtype_d_ *x, 
			     nllspkgtype_d_ *f,
			     const void *params);
  
typedef void (*eval_j_type) (int fstatus, 
			     int n, 
			     int m,
			     const nllspkgtype_d_ *x, 
			     nllspkgtype_d_ *j,
			     const void *params);

/* Perform the nlls solve */
void ral_nlls_d( const int n, const int m, 
		 nllspkgtype_d_ X[],
		 eval_f_type eval_f,
		 eval_j_type eval_j,
		 void *params, 
		 struct nlls_inform_type *status,
		 struct nlls_control_type *options);

#endif

#ifdef __cplusplus
} /* extern "C" */
#endif
