.. default-domain:: C

.. include:: macros.rst

How to use the package
======================

Overview
--------

Calling sequences
^^^^^^^^^^^^^^^^^

Function signatures are definied in a header file

.. code-block:: C 
 
   #include "ral_nlls.h"
 

The user can then call one of the procedures:

|nlls_default_options| initializes solver options to default values

|nlls_solve| solves the non-linear least squares problem.

|nlls_init_workspace| initializes a workspace for use with |nlls_iterate|.

|nlls_iterate| performs one iteration of the non-linear least squares
solver.

|nlls_free_workspace| frees memory allocated by a call to |nlls_init_workspace|.

The calling sequences of these subroutines are outlined in :ref:`arg-lists`.


The derived data types
^^^^^^^^^^^^^^^^^^^^^^

For each problem, the user must employ the derived types
defined by the module to declare scalars of the types |nlls_options| and 
|nlls_inform|.

The following pseudocode illustrates this.

.. code-block:: C

   #include "ral_nlls.h"
   ...
   struct ral_nlls_options options;
   struct ral_nlls_inform inform;
   ...


The members of |nlls_options| and |nlls_inform| are explained below in :ref:`data_types`.

.. _arg-lists:

Argument lists and calling sequences
------------------------------------

The term **package type** is used to mean ``double``. 

To initialize members of |nlls_options| to default values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To initialize the value of |nlls_options|, the user **must** make a call to the 
following subroutine.  Failure to do so will result in undefined behaviour.

.. c:function:: void ral_nlls_default_options(struct spral_lsmr_options *options)

   :param options: data structure to be initialised.

To solve the non-linear least squares problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../common/subroutines.rst

.. c:function:: void nlls_solve(int n, int n, int m, double X[], ral_nlls_eval_r_type eval_r,ral_nlls_eval_j_type eval_J, ral_nlls_eval_hf_type eval_Hf, void* params, struct nlls_options const* options, struct nlls_inform* inform, double weights[], nlls_eval_HP_type eval_HP, double lower_bounds[], double upper_bounds[])

   Solves the non-linear least squares problem.
   
   :param n: |n|

   :param m: |m|
		      
   :param X: |X|

   :param eval_r: |eval_r_desc| 
			
   :param eval_J: |eval_J_desc|

   :param eval_Hf: |eval_Hf_desc|

   :param params: |params|

   :param options: |options|

   :param inform:  |inform|

   :param weights: |weights|

   :param eval_HP: |eval_HP_desc|

   :param lower_bounds: |lower_bounds|

   :param upper_bounds: |upper_bounds|


To initialize a workspace for use with |nlls_iterate|
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Prior to the first call to |nlls_iterate|, the workspace must be
initialised by a call to the following subroutine:

.. c:function:: void ral_nlls_init_workspace(void **workspace)

   :param workspace:  will, on return, be allocated and initialised using Fortran intrinsics.  To avoid a memory leak, it must be freed through a call to |nlls_free_workspace|.


To iterate once
^^^^^^^^^^^^^^^

.. c:function:: void ral_nlls_iterate(int n, int m, double X[], void* workspace, ral_nlls_eval_r_type eval_r, ral_nlls_eval_j_type eval_J, ral_nlls_eval_hf_type eval_Hf, void* params, struct nlls_options const* options, struct nlls_inform* inform, double weights[], nlls_eval_HP_type eval_HP, double lower_bounds[], double upper_bounds[])
		  
   A call of this form allows the user to step through the solution process one
   iteration at a time.

   **n**, **m**, **eval_F**, **eval_J**, **eval_HF**, **params**, **info**,
   **options**, **weights**, **eval_HP**, **lower_bounds** and **upper_bounds**
   are as in the desciption of |nlls_solve|.

   :param  X: |iterate_X|

   :param workspace: is workspace allocated and initialised through a previous call to |nlls_init_workspace|.

The user may use the components ``convergence_normf`` and
``convergence_normg`` and ``converge_norms`` in |nlls_inform| to determine whether the iteration has
converged.

To free a workspace when it is no longer required
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Memory allocated by |nlls_init_workspace| may be freed by a call to the following subroutine:

.. c:function:: void ral_nlls_free_workspace(void **workspace)

   :param workspace:  is the workspace to be freed. On exit it will be set to ``NULL``.

.. _user-routines:

User-supplied function evaluation routines
------------------------------------------

In order to evaluate the function, Jacobian and Hessian at a point, the user
must supply callback functions that perform this operation that the code
**RALFit** will call internally.

In order to pass user-defined data into the evaluation calls, the parameter
``params`` is passed unaltered to the callback functions. Typically this
will be a pointer to a user defined structure that stores the data to be fitted.

For evaluating the function :math:`{\bm r} ( {\bm x} )`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A subroutine must be supplied to calculate :math:`{\bm r} ( {\bm x} )`
for a given vector :math:`{\bm x}`. It must have the following signature:

.. c:function:: int eval_r (int n, int m, void const* params, double const* x, double* r)

   :param n: |eval_r_n|

   :param m: |eval_r_m|
		      
   :param params: |eval_r_params|

   :param x: |eval_r_X|
	    
   :param r: |eval_r_r|

   :param status: |eval_r_status|

For evaluating the function :math:`{\bm J} = \nabla  {\bm r} ( {\bm x} )`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A subroutine must be supplied to calculate
:math:`{\bm J} = \nabla  {\bm r} ( {\bm x} )` for a given vector
:math:`{\bm x}`. 
It must have the following signature:

.. c:function:: int eval_J (int n, int m, void const* params, double const* x, double* J)
   
   :param n: |eval_J_n|
		      
   :param m: |eval_J_m|
		      
   :param params: |eval_J_params|

   :param x: |eval_J_X|

   :param J: |eval_J_r|

   :param status: |eval_J_status|


For evaluating the function :math:`Hf = \sum_{i=1}^m r_i( {\bm x} )  {\bm W} \nabla^2 r_i( {\bm x} )`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A subroutine must be supplied to calculate
:math:`Hf = \sum_{i=1}^m ( {\bm r} )_i \nabla^2 r_i( {\bm x} )` for
given vectors :math:`{\bm x} \in \mathbb{R}^n` and
:math:`{\bm r} \in \mathbb{R}^m`; here :math:`( {\bm r} )_i` denotes
the :math:`i`\ th component of the vector :math:`{\bm r}`. 
It must have the following signature:

.. c:function:: int eval_Hf (int n, int m, void const* params, double const* x, double const* r, double* Hf)
   
   :param n: |eval_Hf_n|

   :param m: |eval_Hf_m|
		      
   :param params: |eval_Hf_params|

   :param X: |eval_Hf_X|

   :param r: |eval_Hf_r|

   :param Hf: |eval_Hf_Hf|

   :param status: |eval_Hf_status|

For evaluating the function :math:`P({\bm x},{\bm y}) := ( H_1({\bm x}){\bm y} \dots  H_m({\bm x}){\bm y})`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A subroutine may be supplied to calculate
:math:`P({\bm x},{\bm y}) := ( H_1({\bm x}){\bm y} \dots  H_m({\bm x}){\bm y})` for
given vectors :math:`{\bm x}, {\bm y} \in \mathbb{R}^n`. The
subroutine must implement the following interface:

.. c:function:: int eval_HP (int n, int m, double const* x, double const* y, double* HP, void const*params)

   :param n: |eval_HP_n|

   :param m: |eval_HP_m|
		      
   :param x: |eval_HP_x|

   :param y: |eval_HP_y|

   :param HP: |eval_HP_HP|
			   
   :param params: |eval_HP_params|

.. _data_types:

Data types
----------

The derived data type for holding options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../common/options.rst


.. c:type:: nlls_options
	    
   This is used to hold controlling data. The components  are automatically given default values in the definition of the type.
   
   **Printing Controls**
   
   .. c:member:: int out 

		 |out|
		 Default is 6.

   .. c:member:: int print_level

		 |print_level|
		 
		 .. include:: ../common/options_print_level.txt
			      
		 Default is 0.

   .. c:member:: bool print_options

     |print_options|

     Default is false.

   .. c:member:: int print_header

     |print_header|

     Default is 30.
      
   **Choice of Algorithm**

   .. c:member:: int model

		 |model| 
		 
		 .. include:: ../common/options_model.txt
			      
		 Default is 3.

   .. c:member:: int type_of_method
		 
		 |type_of_method|
		 
		 .. include:: ../common/options_type_of_method.txt
			      
		 Default is 1.

   .. c:member:: int nlls_method

		 |nlls_method|
				       
		 .. include:: ../common/options_nlls_method.txt
		    
		 Default is 4.
						    
   .. c:member:: bool exact_second_derivatives 

		 |exact_second_derivatives|
		 Default is false.


   **Stopping rules**

   .. c:member:: int maxit 
		 
		 |maxit|
		 Default is 100

   .. c:member:: double stop_g_absolute
		 
		 |stop_g_absolute|
		 Defualt is 1e-5.
					   
   .. c:member:: double stop_g_relative

		 |stop_g_relative|
		 Default is 1e-8.
					   
   .. c:member:: double stop_f_absolute
		 
		 |stop_f_absolute|
		 Default is 1e-5.
					   
   .. c:member:: double stop_f_relative 
		 
		 |stop_f_relative|
		 Default is 1e-8.

   .. c:member:: double stop_s 
		 
		 |stop_s|
		 Default is ``eps``.
   
   **Trust region radius/regularization behaviour**
      
   .. c:member:: int relative_tr_radius 
		 
		 |relative_tr_radius|
		 Default is 0.
					      
   .. c:member:: double initial_radius_scale 

		 |initial_radius_scale|
		 Default is 1.0.

   .. c:member:: double initial_radius

		 |initial_radius|
		 Default is 100.0.

   .. c:member:: double maximum_radius

		 |maximum_radius|
		 Default is 1e8.
		 
   .. c:member:: double eta_successful 
      
		 |eta_successful|
		 Default is 1e-8.
   .. success_but_reduce is also available, but not documented

   .. c:member:: double eta_very_successful
		 
		 |eta_very_successful|
		 Default is 0.9.

   .. c:member:: double eta_too_successful

		 |eta_too_successful|
		 Default is 2.0.

   .. c:member:: double radius_increase 
		 
		 |radius_increase|
		 Default is 2.0.

   .. c:member:: double radius_reduce

		 |radius_reduce|
		 Default is 0.5.

   .. c:member:: int tr_update_strategy

		 |tr_update_strategy|

		 .. include:: ../common/options_tr_update_strategy.txt
		    
		 Default is 1.

   .. c:member:: double reg_order
		 
		 |reg_order|
		 Default is 0.0.
		
   **Scaling options**
   
   .. c:member:: int scale

		 |scale|
		 Default is 1.
			

   .. c:member:: bool scale_trim_max

		 |scale_trim_max|
		 Default is true.

   .. c:member:: double scale_max
      
		 |scale_max|
		 Default is 1e11.

   .. c:member:: bool scale_trim_min

		 |scale_trim_min|
		 Default is true.

   .. c:member:: double scale_min
		 
		 |scale_min|
		 Default is 1e-11.

   .. c:member:: bool scale_require_increase

		 |scale_require_increase|
		 Default is false.
      

   **Hybrid method options**  These options are used if ``model=3``

   .. c:member:: double hybrid_switch

		 |hybrid_switch|
		 Default is 0.1.

   .. c:member:: double hybrid_tol

		 |hybrid_tol|
		 Default is 2.0.	 

   .. c:member:: int hybrid_switch_its
		 
		 |hybrid_switch_its|
		 Default is 1.

   **Newton-Tensor options** These options are used if ``model=4``
					     
   .. c:member:: int inner_method
		 
		 |inner_method| 

		 .. include:: ../common/options_inner_method.txt		      
		 Default is 2.

   **More-Sorensen options**  These options are used if ``nlls_method=3``
   
   .. c:member:: int more_sorensen_maxits
		 
		 |more_sorensen_maxits|
		 Default is 500.

   .. c:member:: double more_sorensen_shift
		 		 
		 |more_sorensen_shift|
		 Default is 1e-13.

   .. c:member:: double more_sorensen_tiny
		 
		 |more_sorensen_tiny|
		 Default is 10.0 * ``eps``.

   .. c:member:: double more_sorensen_tol

		 |more_sorensen_tol|
		 Default is 1e-3.

   **Box bound options**  These options are used if box constraints are included.

   .. c:member:: integer box_nFref_max

		 |box_nFref_max|
		 Default is 4.

   .. c:member:: real box_gamma

		 |box_gamma|
		 Default is 0.9995.
		 
   .. c:member:: real box_decmin

		 |box_decmin|
		 Default is 2.0e-16.

   .. c:member:: real box_bigbnd

		 |box_bigbnd|
		 Default is 1.0e20.

   .. c:member:: real box_wolfe_descent

		 |box_wolfe_descent|
		 Default is 1.0e-4.

   .. c:member:: real box_wolfe_curvature

		 |box_wolfe_curvature|
		 Default is 0.9.

   .. c:member:: real box_kanzow_power

		 |box_kanzow_power|
		 Default is 2.1.

   .. c:member:: real box_kanzow_descent

		 |box_kanzow_descent|
		 Default is 1.0e-8.

   .. c:member:: real box_quad_model_descent

		 |box_quad_model_descent|
		 Default is 1.0e-8.

   .. c:member:: logical box_tr_test_step

		 |box_tr_test_step|
		 Default is true.

   .. c:member:: logical box_wolfe_test_step

		 |box_wolfe_test_step|
		 Default is true.

   .. c:member:: real box_tau_min

		 |box_tau_min|
		 Default is 0.25.

   .. c:member:: real box_tau_descent

		 |box_tau_descent|
		 Default is 1.0e-4.

   .. c:member:: integer box_max_ntrfail

		 |box_max_ntrfail|
		 Default is 2.

   .. c:member:: integer box_quad_match

		 |box_quad_match|
		 Default is 1.

   .. c:member:: real box_alpha_scale

		 |box_alpha_scale|
		 Default is 1.0.

   .. c:member:: real box_Delta_scale

		 |box_Delta_scale|
		 Default is 2.0.

   .. c:member:: real box_tau_wolfe

		 |box_tau_wolfe|
		 Default is 0.3.

   .. c:member:: real box_tau_tr_step

		 |box_tau_tr_step|
		 Default is 0.3.

   .. c:member:: integer box_ls_step_maxit

		 |box_ls_step_maxit|
		 Default is 20.

   .. c:member:: integer box_lineseach_type

		 |box_linesearch_type|
		 Default is 1.

					       .. include::  ../common/options_linesearch_type.txt
						  
   **Other options**
					     
   .. c:member:: bool output_progress_vectors
		 
		 |output_progress_vectors|
		 Default is false.

   **Internal options to help solving a regularized problem implicitly**

   .. c:member:: int regularization 
      
		 |regularization|
					  
		 .. include:: ../common/options_regularization.txt
		 Default is 0.
					  
   .. c:member:: double regularization_term 

		 |regularization_term|
		 Default is 0.0.
						
   .. c:member:: double regularization_power
      
		 |regularization_power|
		 Default is 0.0.


The derived data type for holding information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../common/inform.rst

.. c:type:: nlls_inform
	    
   This is used to hold information about the progress of the algorithm.

   .. c:member:: int status
		 
		 |status|

   .. c:member:: character error_message(81)

		 |error_message|

   .. c:member:: int alloc_status 
      
		 |alloc_status|

   .. c:member:: character bad_alloc(81)
		 
		 |bad_alloc|

   .. c:member:: int iter 

		 |iter|

   .. c:member:: int f_eval 
		 
		 |f_eval|
		      
   .. c:member:: int g_eval 
		 
		 |g_eval|

   .. c:member:: int h_eval
      
		 |h_eval|

   .. c:member:: int hp_eval
      
		 |hp_eval|
		 
   .. c:member:: int convergence_normf
		 
		 |convergence_normf|

   .. c:member:: int convergence_normg
		 
		 |convergence_normg|
				
   .. c:member:: int convergence_norms 

		 |convergence_norms|
				   
   .. c:member:: double resvec(iter+1)
		 
		 |resvec|

   .. c:member:: double gradvec(iter+1)
		 
		 |gradvec|

   .. c:member:: double obj
		 
		 |obj|

   .. c:member:: double norm_g
			       
		 |norm_g|

   .. c:member:: double scaled_g

		 |scaled_g|

   .. c:member:: int external_return
			       
		 |external_return|
   
   .. c:member:: character external_name(81)

		 |external_name|

   .. c:member:: double step

		 |step|
   

The workspace derived data type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:type:: nlls_workspace
	    
   This is used to save the state of the algorithm in between calls to |nlls_iterate|, 
   and must be used if that subroutine is required. It's components are not 
   designed to be accessed by the user.

.. _errors:

Warning and error messages
--------------------------

A successful return from a subroutine in the package is indicated by ``status``
in |nlls_inform| having the value zero.  
A non-zero value is asscociated with an error message,
which will be output on ``error`` in |nlls_inform|.

.. include:: ../common/errors.rst

   
