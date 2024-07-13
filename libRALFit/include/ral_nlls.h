/* COPYRIGHT (c) 2015 Science and Technology Facilities Council (STFC)
 * All rights reserved
 * Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
 */
#ifdef __cplusplus
extern "C" {
#else
#include <stdbool.h>
#endif

#ifndef ral_nlls_h
#define ral_nlls_h

#ifndef ral_int
#define ral_int int
#endif

typedef double ral_nllspkgtype_d_;

#ifndef ral_nlls_options_type
#define ral_nlls_default_options ral_nlls_default_options_d
#define ral_nlls_options ral_nlls_options_d
#define ral_nlls_inform ral_nlls_inform_d
#define nlls_solve nlls_solve_d
#define ral_nlls_init_workspace ral_nlls_init_workspace_d
#define ral_nlls_iterate ral_nlls_iterate_d
#define nlls_strerror nlls_strerror_d
#define ral_nlls_free_workspace ral_nlls_free_workspace_d
#endif


/* Derived type to hold option parameters for ral_nlls */
struct ral_nlls_options_d {
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
  ral_nllspkgtype_d_ stop_g_absolute; /* absolute stopping tolerance for gradient*/
  ral_nllspkgtype_d_ stop_g_relative; /* relative stopping tolerance for gradient*/
  ral_nllspkgtype_d_ stop_f_absolute; /* absolute stopping tolerance for function*/
  ral_nllspkgtype_d_ stop_f_relative; /* relative stopping tolerance for function*/
  ral_nllspkgtype_d_ stop_s; /* stopping test on change in solution */
  ral_int relative_tr_radius; /* non-zero if we scale the initial tr radius */
  ral_nllspkgtype_d_ initial_radius_scale; /* what should we scale tr rad by? */
  ral_nllspkgtype_d_ initial_radius; /* initial trust region radius */
  ral_nllspkgtype_d_ base_regularization; /* if needed, the base reg parameter */
  ral_int regularization; /* what method used to solve the regularized nlls problem */
  ral_nllspkgtype_d_ regularization_term; /* reg weight used if regularization != 0 */
  ral_nllspkgtype_d_ regularization_power; /* reg order used if regularization != 0 */
  ral_nllspkgtype_d_ maximum_radius; /* maximux trust region radius */
  ral_nllspkgtype_d_ eta_successful; /* trust region step successful level */
  ral_nllspkgtype_d_ eta_success_but_reduce;
  ral_nllspkgtype_d_ eta_very_successful; /* trust region step very successful */
  ral_nllspkgtype_d_ eta_too_successful; /* trust region step too successful */
  ral_nllspkgtype_d_ radius_increase; /* how much to increase the radius by? */
  ral_nllspkgtype_d_ radius_reduce; /* how much to reduce the radius by? */
  ral_nllspkgtype_d_ radius_reduce_max; /* max amount to reduce the radius by */
  ral_int tr_update_strategy; /* 1: step function, 2: continuous */
  ral_nllspkgtype_d_ hybrid_switch;
  bool exact_second_derivatives;
  bool subproblem_eig_fact;
  bool use_ews_subproblem;
  bool force_min_eig_symm;
  ral_int scale; /* 0: don't scale, 1: norm of J, 2: norm of Hessian, 3: eigs */
  ral_nllspkgtype_d_ scale_max; /* max before we trim scaling */
  ral_nllspkgtype_d_ scale_min; /* min before we trim scaling */
  bool scale_trim_min; /* if min attained, trim? (or set to 1) */
  bool scale_trim_max; /* if max attained, trim? (or set to 1) */
  bool scale_require_increase; /* scaling matrix must increase to update */
  bool setup_workspaces;
  bool remove_workspaces;
  ral_int more_sorensen_maxits;
  ral_nllspkgtype_d_ more_sorensen_shift;
  ral_nllspkgtype_d_ more_sorensen_tiny;
  ral_nllspkgtype_d_ more_sorensen_tol;
  ral_nllspkgtype_d_ hybrid_tol;
  ral_int hybrid_switch_its;
  ral_nllspkgtype_d_ reg_order;
  ral_int inner_method;
  bool output_progress_vectors;
  bool update_lower_order;
  bool Fortran_Jacobian;
  ral_int box_nFref_max;
  ral_nllspkgtype_d_ box_gamma;
  ral_nllspkgtype_d_ box_decmin;
  ral_nllspkgtype_d_ box_bigbnd;
  ral_nllspkgtype_d_ box_wolfe_descent;
  ral_nllspkgtype_d_ box_wolfe_curvature;
  ral_nllspkgtype_d_ box_kanzow_power;
  ral_nllspkgtype_d_ box_kanzow_descent;
  ral_nllspkgtype_d_ box_quad_model_descent;
  bool box_tr_test_step;
  bool box_wolfe_test_step;
  ral_nllspkgtype_d_ box_tau_descent;
  ral_int box_max_ntrfail;
  ral_int box_quad_match;
  ral_nllspkgtype_d_ box_alpha_scale;
  ral_nllspkgtype_d_ box_Delta_scale;
  ral_nllspkgtype_d_ box_tau_min;
  ral_int box_ls_step_maxit;
  ral_int box_linesearch_type;
  ral_nllspkgtype_d_ fd_step;
  ral_int check_derivatives;
  ral_nllspkgtype_d_ derivative_test_tol;
};

struct ral_nlls_inform_d {
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
  ral_nllspkgtype_d_ resinf;
  ral_nllspkgtype_d_ gradinf;
  ral_nllspkgtype_d_ obj;
  ral_nllspkgtype_d_ norm_g;
  ral_nllspkgtype_d_ scaled_g;
  ral_int external_return;
  char external_name[81];
  ral_nllspkgtype_d_ step;
  ral_int ls_step_iter;
  ral_int f_eval_ls;
  ral_int g_eval_ls;
  ral_int pg_step_iter;
  ral_int f_eval_pg;
  ral_int g_eval_pg;
  ral_int fd_f_eval;
};

/* Set default values of options */
void ral_nlls_default_options_d( struct ral_nlls_options_d *options );

/* define the eval_f_type */
typedef ral_int (*ral_nlls_eval_r_type) (
              ral_int n,
              ral_int m,
              void *params,
              const ral_nllspkgtype_d_ *x,
              ral_nllspkgtype_d_ *f
              );

/* define the eval_j_type */
typedef ral_int (*ral_nlls_eval_j_type) (
              ral_int n,
              ral_int m,
              void *params,
              const ral_nllspkgtype_d_ *x,
              ral_nllspkgtype_d_ *j
              );

/* define the eval_hf_type */
typedef ral_int (*ral_nlls_eval_hf_type) (
               ral_int n,
               ral_int m,
               void *params,
               const ral_nllspkgtype_d_ *x,
               const ral_nllspkgtype_d_ *f,
               ral_nllspkgtype_d_ *hf
               );

/* define the eval_hp_type */
typedef ral_int (*ral_nlls_eval_hp_type) (
               ral_int n,
               ral_int m,
               const ral_nllspkgtype_d_ *x,
               const ral_nllspkgtype_d_ *y,
               ral_nllspkgtype_d_ *hp,
	       void *params
               );


/* Perform the nlls solve */
  void nlls_solve_d( ral_int n, ral_int m,
		     ral_nllspkgtype_d_ X[],
		     ral_nlls_eval_r_type eval_r,
		     ral_nlls_eval_j_type eval_j,
		     ral_nlls_eval_hf_type eval_hf,
		     void const* params,
		     struct ral_nlls_options const* options,
		     struct ral_nlls_inform *status,
		     ral_nllspkgtype_d_ weights[],
		     ral_nlls_eval_hp_type eval_hp,
		     ral_nllspkgtype_d_ lower_bounds[],
		     ral_nllspkgtype_d_ upper_bounds[]
		     );
/* Initialise a workspace for use with ral_nlls_iterate_d() */
  void ral_nlls_init_workspace_d(void **w, void **iw);
  /* Perform a single iteration */
  void ral_nlls_iterate_d(
			ral_int n,
			ral_int m,
			ral_nllspkgtype_d_ X[],
			void *w,
			ral_nlls_eval_r_type eval_r,
			ral_nlls_eval_j_type eval_j,
			ral_nlls_eval_hf_type eval_hf,
			void const* params,
			struct ral_nlls_options const* options,
			struct ral_nlls_inform *status,
			ral_nllspkgtype_d_ weights[],
			ral_nlls_eval_hp_type eval_hp,
			ral_nllspkgtype_d_ lower_bounds[],
			ral_nllspkgtype_d_ upper_bounds[]
			);
  /* get the error string */
  void nlls_strerror_d(
		       struct ral_nlls_inform *status,
		       char error_string[81]
		       );
/* Free memory from a workspace */
  void ral_nlls_free_workspace_d(void **w);

#endif

#ifdef __cplusplus
} /* extern "C" */
#endif
