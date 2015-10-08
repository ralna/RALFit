/* COPYRIGHT (c) 2015 Science and Technology Facilities Council (STFC)
 * All rights reserved
 */

#ifndef ral_nlls_h
#define ral_nlls_h

#ifndef NLLS_control_type
#define NLLS_default_control NLLS_default_control_d
#define NLLS_control_type NLLS_control_type_d
#define NLLS_inform_type NLLS_inform_type_d
#define RAL_NLLS RAL_NLLS_d
#endif

typedef double NLLSpkgtype_d_;

/* Derived type to hold control parameters for ral_nlls */
struct NLLS_control_type_d {
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

struct NLLS_inform_type_d {
  int status; /* flag */
};

/* Set default values of control */
void NLLS_default_control_d( struct NLLS_control_type *options );

/* Perform the NLLS solve */
void RAL_NLLS_d( const int n, const int m, 
		 NLLSpkgtype_d_ X[],
		 int Work_int[], const int len_work_int,
		 NLLSpkgtype_d_ Work_real[], const int len_work_real,
		 void (*eval_F)(int,double,double),
		 void (*eval_J)(int,double,double),
		 struct NLLS_inform_type *status,
		 struct NLLS_control_type *options);

#endif
