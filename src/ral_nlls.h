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
/*#define ral_nlls ral_nlls_d*/
#define ral_nlls_int_func ral_nlls_int_func_d
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
  double stop_g_absolute; /* absolute stopping tolerance */
  double stop_g_relative; /* relative stopping tolerance */
  double initial_radius; /* initial trust region radius */
  double maximum_radius; /* maximux trust region radius */
  double eta_successful; /* trust region step successful level */
  double eta_very_successful; /* trust region step very successful */
  double eta_too_successful; /* trust region step too successful */
  double radius_increase; /* how much to increase the radius by? */
  double radius_reduce; /* how much to reduce the radius by? */
  double radius_reduce_max; /* max amount to reduce the radius by */
};

struct nlls_inform_type_d {
  int status; /* flag */
};

/* Set default values of control */
void nlls_default_control_d( struct nlls_control_type *options );

/* Perform the nlls solve */
/* 
   void ral_nlls_d( const int n, const int m, 
		 nllspkgtype_d_ X[],
		 int Work_int[], const int len_work_int,
		 nllspkgtype_d_ Work_real[], const int len_work_real,
		 void (*eval_F)(int,double,double),
		 void (*eval_J)(int,double,double),
		 struct nlls_inform_type *status,
		 struct nlls_control_type *options);
*/
void ral_nlls_int_func_d( const int n, const int m, 
		 nllspkgtype_d_ X[],
		 struct nlls_inform_type *status,
		 struct nlls_control_type *options);

#endif

#ifdef __cplusplus
} /* extern "C" */
#endif
