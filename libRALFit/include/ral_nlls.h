/* COPYRIGHT (c) 2015 Science and Technology Facilities Council (STFC)
 * All rights reserved
 * Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
 */
#ifdef __cplusplus
extern "C" {
#else
#include <stdint.h>
#include <stdbool.h>
#endif

#ifndef ral_nlls_h
#define ral_nlls_h

#ifndef ral_real
#ifdef SINGLE_PRECISION
   #define ral_real float
#else
   #define ral_real double
#endif
#endif

#ifdef SINGLE_PRECISION
   #define PREC(fn) fn ## _s
#else
   #define PREC(fn) fn ## _d
#endif

#ifndef ral_int
#ifdef BUILD_ILP64
   // To build ilp64 also the default integer kind must be set to 8 bytes
   #define ral_int int64_t
#else
   #define ral_int int32_t
#endif
#endif


/* Derived type to hold option parameters for ral_nlls */
struct PREC(ral_nlls_options) {
  ral_int f_arrays; /* Use 1-based indexing if true(!=0) else 0-based */
  ral_int out;   /* Fortran output stream for general messages */
  ral_int print_level; /* levels of print output */
  bool print_options; /* print all options? */
  ral_int print_header; /*=k; print one-liner header every k-th iteration */
  ral_int maxit; /* maximum number of iterations */
  ral_int model; /* what model to use? */
  ral_int type_of_method; /* what method to use? */
  ral_int nlls_method; /* what nlls method to use? */
  bool allow_fallback_method; /* switch nlls method if chosen one fails? */
  ral_int lls_solver; /* which lls solver to use? */
  ral_real stop_g_absolute; /* absolute stopping tolerance for gradient*/
  ral_real stop_g_relative; /* relative stopping tolerance for gradient*/
  ral_real stop_f_absolute; /* absolute stopping tolerance for function*/
  ral_real stop_f_relative; /* relative stopping tolerance for function*/
  ral_real stop_s; /* stopping test on change in solution */
  ral_int relative_tr_radius; /* non-zero if we scale the initial tr radius */
  ral_real initial_radius_scale; /* what should we scale tr rad by? */
  ral_real initial_radius; /* initial trust region radius */
  ral_real base_regularization; /* if needed, the base reg parameter */
  ral_int regularization; /* what method used to solve the regularized nlls problem */
  ral_real regularization_term; /* reg weight used if regularization != 0 */
  ral_real regularization_power; /* reg order used if regularization != 0 */
  ral_real maximum_radius; /* maximux trust region radius */
  ral_real eta_successful; /* trust region step successful level */
  ral_real eta_success_but_reduce;
  ral_real eta_very_successful; /* trust region step very successful */
  ral_real eta_too_successful; /* trust region step too successful */
  ral_real radius_increase; /* how much to increase the radius by? */
  ral_real radius_reduce; /* how much to reduce the radius by? */
  ral_real radius_reduce_max; /* max amount to reduce the radius by */
  ral_int tr_update_strategy; /* 1: step function, 2: continuous */
  ral_real hybrid_switch;
  bool exact_second_derivatives;
  bool subproblem_eig_fact;
  bool use_ews_subproblem;
  bool force_min_eig_symm;
  ral_int scale; /* 0: don't scale, 1: norm of J, 2: norm of Hessian, 3: eigs */
  ral_real scale_max; /* max before we trim scaling */
  ral_real scale_min; /* min before we trim scaling */
  bool scale_trim_min; /* if min attained, trim? (or set to 1) */
  bool scale_trim_max; /* if max attained, trim? (or set to 1) */
  bool scale_require_increase; /* scaling matrix must increase to update */
  bool setup_workspaces;
  bool remove_workspaces;
  ral_int more_sorensen_maxits;
  ral_real more_sorensen_shift;
  ral_real more_sorensen_tiny;
  ral_real more_sorensen_tol;
  ral_real hybrid_tol;
  ral_int hybrid_switch_its;
  ral_real reg_order;
  ral_int inner_method;
  bool output_progress_vectors;
  bool update_lower_order;
  bool Fortran_Jacobian;
  ral_int box_nFref_max;
  ral_real box_gamma;
  ral_real box_decmin;
  ral_real box_bigbnd;
  ral_real box_wolfe_descent;
  ral_real box_wolfe_curvature;
  ral_real box_kanzow_power;
  ral_real box_kanzow_descent;
  ral_real box_quad_model_descent;
  bool box_tr_test_step;
  bool box_wolfe_test_step;
  ral_real box_tau_descent;
  ral_int box_max_ntrfail;
  ral_int box_quad_match;
  ral_real box_alpha_scale;
  ral_real box_Delta_scale;
  ral_real box_tau_min;
  ral_int box_ls_step_maxit;
  ral_int box_linesearch_type;
  ral_real fd_step;
  ral_int check_derivatives;
  ral_real derivative_test_tol;
};

struct PREC(ral_nlls_inform) {
  ral_int status; /* flag */
  char error_message[81];
  ral_int alloc_status;
  char bad_alloc[81];
  ral_int iter;
  ral_int inner_iter;
  bool inner_iter_success;
  ral_int f_eval;
  ral_int g_eval;
  ral_int h_eval;
  ral_int hp_eval;
  ral_int convergence_normf;
  ral_int convergence_normg;
  ral_int convergence_norms;
  ral_real resinf;
  ral_real gradinf;
  ral_real obj;
  ral_real norm_g;
  ral_real scaled_g;
  ral_int external_return;
  char external_name[81];
  ral_real step;
  ral_int ls_step_iter;
  ral_int f_eval_ls;
  ral_int g_eval_ls;
  ral_int pg_step_iter;
  ral_int f_eval_pg;
  ral_int g_eval_pg;
  ral_int fd_f_eval;
};

/* Set default values of options */
void PREC(ral_nlls_default_options)( struct PREC(ral_nlls_options) *options );

/* define the eval_f_type */
typedef ral_int (PREC(*ral_nlls_eval_r_type)) (
              ral_int n,
              ral_int m,
              void *params,
              const ral_real *x,
              ral_real *f
              );

/* define the eval_j_type */
typedef ral_int (PREC(*ral_nlls_eval_j_type)) (
              ral_int n,
              ral_int m,
              void *params,
              const ral_real *x,
              ral_real *j
              );

/* define the eval_hf_type */
typedef ral_int (PREC(*ral_nlls_eval_hf_type)) (
               ral_int n,
               ral_int m,
               void *params,
               const ral_real *x,
               const ral_real *f,
               ral_real *hf
               );

/* define the eval_hp_type */
typedef ral_int (PREC(*ral_nlls_eval_hp_type)) (
               ral_int n,
               ral_int m,
               const ral_real *x,
               const ral_real *y,
               ral_real *hp,
	       void *params
               );


/* Perform the nlls solve */
  void PREC(nlls_solve)( ral_int n, ral_int m,
		     ral_real X[],
		     PREC(ral_nlls_eval_r_type),
		     PREC(ral_nlls_eval_j_type),
		     PREC(ral_nlls_eval_hf_type),
		     void const* params,
		     struct PREC(ral_nlls_options) const* options,
		     struct PREC(ral_nlls_inform) *status,
		     ral_real weights[],
		     PREC(ral_nlls_eval_hp_type),
		     ral_real lower_bounds[],
		     ral_real upper_bounds[]
		     );
/* Initialise a workspace for use with ral_nlls_iterate_d() */
  void PREC(ral_nlls_init_workspace)(void **w, void **iw);
  /* Perform a single iteration */
  void PREC(ral_nlls_iterate)(
			ral_int n,
			ral_int m,
			ral_real X[],
			void *w,
			PREC(ral_nlls_eval_r_type),
			PREC(ral_nlls_eval_j_type),
			PREC(ral_nlls_eval_hf_type),
			void const* params,
			struct PREC(ral_nlls_options) const* options,
			struct PREC(ral_nlls_inform) *status,
			ral_real weights[],
			PREC(ral_nlls_eval_hp_type),
			ral_real lower_bounds[],
			ral_real upper_bounds[]
			);
  /* get the error string */
  void PREC(nlls_strerror)(
		       struct PREC(ral_nlls_inform) *status,
		       char error_string[81]
		       );
/* Free memory from a workspace */
  void PREC(ral_nlls_free_workspace)(void **w);

// Alias all public objects to selected precision
#define ral_nlls_default_options PREC(ral_nlls_default_options)
#define ral_nlls_options         PREC(ral_nlls_options)
#define ral_nlls_inform          PREC(ral_nlls_inform)
#define nlls_solve               PREC(nlls_solve)
#define ral_nlls_init_workspace  PREC(ral_nlls_init_workspace)
#define ral_nlls_iterate         PREC(ral_nlls_iterate)
#define nlls_strerror            PREC(nlls_strerror)
#define ral_nlls_free_workspace  PREC(ral_nlls_free_workspace)

#endif

#ifdef __cplusplus
} /* extern "C" */
#endif
