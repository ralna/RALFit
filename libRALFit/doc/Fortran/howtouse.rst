.. default-domain:: Fortran

.. include:: macros.rst

How to use the package
======================

Overview
--------

Calling sequences
^^^^^^^^^^^^^^^^^

Access to the package requires a ``USE`` statement

.. code-block:: fortran

   use ral_nlls_double


The user can then call one of the procedures:

|nlls_solve| solves the non-linear least squares problem.

|nlls_iterate| performs one iteration of the non-linear least squares
solver.

The calling sequences of these subroutines are outlined in :ref:`arg-lists`.


The derived data types
^^^^^^^^^^^^^^^^^^^^^^

For each problem, the user must employ the derived types
defined by the module to declare scalars of the types |nlls_inform| and |nlls_options|. 
If |nlls_iterate| is to be used, then a scalar of the type :f:type:`nlls_workspace` must also be defined. 
The following pseudocode illustrates this.

.. code-block:: fortran

   use nlls_module
   !...
   type (NLLS_inform) :: inform
   type (NLLS_options) :: options
   type (NLLS_workspace) :: work ! needed if nlls_iterate to be called
   !...


The components of |nlls_options| and |nlls_inform| are explained below in :ref:`data_types`.

.. _arg-lists:

Argument lists and calling sequences
------------------------------------

We use square brackets to indicate optional arguments, which
follow the argument :f:type:`inform`. Since we reserve the right to add
additional optional arguments in future releases of the code, **we
strongly recommend that all optional arguments be called by keyword, not
by position**.

The term **package type** is used to mean double precision.


To solve the non-linear least squares problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../common/subroutines.rst

.. f:subroutine:: nlls_solve(n,m,X,eval_r,eval_J,eval_Hf,params,options,inform[,weights,eval_HP,lower_bounds,upper_bounds])

   Solves the non-linear least squares problem.
   
   :p integer n [in]: |n|

   :p integer m [in]: |m|
		      
   :p real X(n) [inout]: |X|

   :p procedure eval_r: |eval_r_desc| 
			
   :p procedure eval_J: |eval_J_desc|

   :p procedure eval_Hf: |eval_Hf_desc|

   :p params_base_type params [in]: |params|

   :p nlls_options options [in]: |options|

   :p nlls_inform inform [out]:  |inform|

   :o real weights(n): |weights|

   :o procedure eval_HP: |eval_HP_desc|

   :o real lower_bounds(n): |lower_bounds|

   :o real upper_bounds(n): |upper_bounds|

To iterate once
^^^^^^^^^^^^^^^


.. f:subroutine:: nlls_iterate(n,m,X,eval_r,eval_J,eval_Hf,params,options,inform[,weights,eval_HP,lower_bounds,upper_bounds])
		  
   A call of this form allows the user to step through the solution process one
   iteration at a time.

   **n**, **m**, **eval_F**, **eval_J**, **eval_HF**, **params**, **inform**,
   **options**, **weights**, **eval_HP**, **lower_bounds** and **upper_bounds**
   are as in the desciption of |nlls_solve|.

   :p real X(n) [inout]: |iterate_X|

   :p nlls_workspace w [inout]: is used to store the current state of the iteration and should not be altered by the user.

The user may use the components ``convergence_normf`` and
``convergence_normg`` and ``converge_norms`` in |nlls_inform| to determine whether the iteration has
converged.
   

.. _user-routines:

User-supplied function evaluation routines
------------------------------------------

The user must supply routines to evaluate the residual, Jacobian and Hessian 
at a point.  **RALFit** will call these routines internally.

In order to pass user-defined data into the evaluation calls, :f:type:`params_base_type` is extended to a :f:type:`user_type`, as follows:

.. code-block:: fortran

       type, extends( params_base_type ) :: user_type
          ! code declaring components of user_type
       end type user_type

We recommend this type is wrapped in a module with the user-defined
routines for evaluating the function, Jacobian, and Hessian.

The components of the extended type are accessed through a
``select type`` construct:

.. code-block:: fortran

       select type(params)
       type is(user_type)
         ! code that accesses components of params that were defined within user_type
       end select

For evaluating the function :math:`{\bm r} ( {\bm x} )`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A subroutine must be supplied to calculate :math:`{\bm r} ( {\bm x} )`
for a given vector :math:`{\bm x}`. It must implement the following
interface:

.. code-block:: fortran

    abstract interface
       subroutine eval_r(status, n, m, x, r, params)
          integer, intent(inout) :: status
          integer, intent(in) :: n
          integer, intent(in) :: m
          double precision, dimension(n), intent(in) :: x
	  double precision, dimension(m), intent(out) :: r
          class(params_base_type), intent(in) :: params
       end subroutine eval_r
    end interface

.. f:subroutine:: eval_r(status, n, m, x, r, params)
   
   :p integer status [inout]: |eval_r_status|
			   
   :p integer n [in]: |eval_r_n|

   :p integer m [in]: |eval_r_m|
		      
   :p real X(n) [in]: |eval_r_X|
	    
   :p real r(m) [out]: |eval_r_r|
		    
   :p params_base_type params [in]: |eval_r_params|


For evaluating the function :math:`{\bm J} = \nabla  {\bm r} ( {\bm x} )`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A subroutine must be supplied to calculate
:math:`{\bm J} = \nabla  {\bm r} ( {\bm x} )` for a given vector
:math:`{\bm x}`. It must implement the following interface:

.. code-block:: fortran

    abstract interface
       subroutine eval_J(status, n, m, x, J, params)
          integer, intent(inout) :: status
	  integer, intent(in) :: n
          integer, intent(in) :: m
          double precision, dimension(n), intent(in)  :: x
	  double precision, dimension(n*m), intent(out) :: J
	  class(params_base_type), intent(in) :: params
      end subroutine eval_J
    end interface

.. f:subroutine:: eval_J(status,n,m,x,J,params)
   
   :p integer status [inout]: |eval_J_status|

   :p integer n [in]: |eval_J_n|
		      
   :p integer m [in]: |eval_J_m|

   :p real X(n) [in]: |eval_J_X|

   :p real J(m*n) [out]: |eval_J_r|

   :p params_base_type params [in]: |eval_J_params|


For evaluating the function :math:`Hf = \sum_{i=1}^m r_i( {\bm x} )  {\bm W} \nabla^2 r_i( {\bm x} )`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A subroutine must be supplied to calculate
:math:`Hf = \sum_{i=1}^m ( {\bm r} )_i \nabla^2 r_i( {\bm x} )` for
given vectors :math:`{\bm x} \in \mathbb{R}^n` and
:math:`{\bm r} \in \mathbb{R}^m`; here :math:`( {\bm r} )_i` denotes
the :math:`i`\ th component of the vector :math:`{\bm r}`. The
subroutine must implement the following interface:

.. code-block:: fortran

    abstract interface
       subroutine eval_Hf_type(status, n, m, x, r, Hf, params)
           integer, intent(inout) :: status           
	   integer, intent(in) :: n
           integer, intent(in) :: m
           double precision, dimension(n), intent(in)  :: x
           double precision, dimension(m), intent(in)  :: r
	   double precision, dimension(n*n), intent(out) :: Hf
	   class(params_base_type), intent(in) :: params
         end subroutine eval_Hf_type
    end interface
    :language: fortran

.. f:subroutine:: eval_Hf(status,n,m,x,r,Hf,params)
   
   :p integer status [inout]: |eval_Hf_status|

   :p integer n [in]: |eval_Hf_n|

   :p integer m [in]: |eval_Hf_m|

   :p real X(n) [in]: |eval_Hf_X|

   :p real r(m) [in]: |eval_Hf_r|

   :p real Hf(n*n) [out]: |eval_Hf_Hf|

   :p params_base_type params [in]: |eval_Hf_params|

For evaluating the function :math:`P({\bm x},{\bm y}) := ( H_1({\bm x}){\bm y} \dots  H_m({\bm x}){\bm y})`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A subroutine may be supplied to calculate
:math:`P({\bm x},{\bm y}) := ( H_1({\bm x}){\bm y} \dots  H_m({\bm x}){\bm y})` for
given vectors :math:`{\bm x}, {\bm y} \in \mathbb{R}^n`. The
subroutine must implement the following interface:

.. code-block:: fortran

    abstract interface
       subroutine eval_HP_type(status, n, m, x, y, HP, params)
           integer, intent(inout) :: status
           integer, intent(in) :: n
           integer, intent(in) :: m
           double precision, dimension(n), intent(in)  :: x
           double precision, dimension(n), intent(in)  :: y
           double precision, dimension(n*m), intent(out) :: HP
           class(params_base_type), intent(in) :: params
         end subroutine eval_HP_type
    end interface
    :language: fortran

.. f:subroutine:: eval_HP(status,n,m,x,y,HP,params)

   :p integer status [inout]: |eval_HP_status|
			   
   :p integer n [in]: |eval_HP_n|

   :p integer m [in]: |eval_HP_m|
		      
   :p real x(n) [in]: |eval_HP_x|

   :p real y(n) [in]: |eval_HP_y|

   :p real HP(n*m) [out]: |eval_HP_HP|
			   
   :p params_base_type params [in]: |eval_HP_params|
				 
.. _data_types:

Data types
----------

The derived data type for holding options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../common/options.rst


.. f:type:: nlls_options
	    
   This is used to hold controlling data. The components  are automatically given default values in the definition of the type.
   
   **Printing Controls**

   :f integer out [default=6]: |out|

   :f integer print_level [default=0]: |print_level|

				       .. include:: ../common/options_print_level.txt

   :f logical print_options [default=false]: |print_options|

   :f integer print_header [default=30]: |print_header|
      
   **Choice of Algorithm**

   :f integer model [default=3]: |model| 

	     			 .. include:: ../common/options_model.txt

   :f integer type_of_method [default=1]: |type_of_method|

					  .. include:: ../common/options_type_of_method.txt

   :f integer nlls_method [default=4]: |nlls_method|
				       
				       .. include:: ../common/options_nlls_method.txt
						    
   :f logical exact_second_derivatives [default=false]: |exact_second_derivatives|


   **Stopping rules**

   :f integer maxit [default=100]: |maxit|

   :f real stop_g_absolute [default=1e-5]: |stop_g_absolute|
					   
   :f real stop_g_relative [default=1e-8]: |stop_g_relative|
					   
   :f real stop_f_absolute [default=1e-5]: |stop_f_absolute|
					   
   :f real stop_f_relative [default=1e-8]: |stop_f_relative|

   :f real stop_s [default=eps]: |stop_s|
   
   **Trust region radius/regularization behaviour**
      
   :f integer relative_tr_radius [default=0]: |relative_tr_radius|
					      
   :f real initial_radius_scale [default=1.0]: |initial_radius_scale|

   :f real initial_radius [default=100.0]: |initial_radius|

   :f real maximum_radius [default=1e8]: |maximum_radius|

   :f real eta_successful [default=1e-8]: |eta_successful|
					  .. success_but_reduce is also available, but not documented

   :f real eta_very_successful [default=0.9]: |eta_very_successful|

   :f real eta_too_successful [default=2.0]: |eta_too_successful|

   :f real radius_increase [default=2.0]: |radius_increase|

   :f real radius_reduce [default=0.5]: |radius_reduce|

   :f integer tr_update_strategy [default=1]: |tr_update_strategy|

					      .. include:: ../common/options_tr_update_strategy.txt

   :f real reg_order [default=0.0]: |reg_order|
							   
   **Scaling options**
   
   :f integer scale [default=1]: |scale|
			

   :f logical scale_trim_max [default=true]: |scale_trim_max|

   :f real scale_max [default=1e11]: |scale_max|

   :f logical scale_trim_min [default=true]: |scale_trim_min|

   :f real scale_min [default=1e-11]: |scale_min|

   :f logical scale_require_increase [default=false]: |scale_require_increase|

   **Hybrid method options**  These options are used if ``model=3``

   :f real hybrid_switch [default=0.1]: |hybrid_switch|

   :f real hybrid_tol [default=2.0]: |hybrid_tol|

   :f integer hybrid_switch_its [default=1]: |hybrid_switch_its|

   **Newton-Tensor options** These options are used if ``model=4``
					     
   :f integer inner_method [default=2]: |inner_method| 

					.. include:: ../common/options_inner_method.txt

   **More-Sorensen options**  These options are used if ``nlls_method=3``
   
   :f integer more_sorensen_maxits [default=500]: |more_sorensen_maxits|

   :f real more_sorensen_shift [default=1e-13]: |more_sorensen_shift|

   :f real more_sorensen_tiny [default=10.0*eps]: |more_sorensen_tiny|

   :f real more_sorensen_tol [default=1e-3]: |more_sorensen_tol|

   **Box Bound Options** These options are used if box constraints are included.
   
   :f integer box_nFref_max [default=4]: |box_nFref_max|

   :f real box_gamma [default=0.9995]: |box_gamma|

   :f real box_decmin [default=2.0e-16]: |box_decmin|

   :f real box_bigbnd [default=1.0e20]: |box_bigbnd|

   :f real box_wolfe_descent [default=1.0e-4]: |box_wolfe_descent|

   :f real box_wolfe_curvature [default=0.9]: |box_wolfe_curvature|

   :f real box_kanzow_power [default=2.1]: |box_kanzow_power|

   :f real box_kanzow_descent [default=1.0e-8]: |box_kanzow_descent|

   :f real box_quad_model_descent [default=1.0e-8]: |box_quad_model_descent|

   :f logical box_tr_test_step [default=true]: |box_tr_test_step|

   :f logical box_wolfe_test_step [default=true]: |box_wolfe_test_step|

   :f real box_tau_min [default=0.25]: |box_tau_min|

   :f real box_tau_descent [default=1.0e-4]: |box_tau_descent|

   :f integer box_max_ntrfail [default=2]: |box_max_ntrfail|

   :f integer box_quad_match [default=1]: |box_quad_match|

   :f real box_alpha_scale [default=1.0]: |box_alpha_scale|

   :f real box_Delta_scale [default=2.0]: |box_Delta_scale|

   :f real box_tau_wolfe [default=0.3]: |box_tau_wolfe|

   :f real box_tau_tr_step [default=0.3]: |box_tau_tr_step|

   :f integer box_ls_step_maxit [default=20]: |box_ls_step_maxit|

   :f integer box_lineseach_type [default=1]: |box_linesearch_type|

					       .. include:: ../common/options_linesearch_type.txt

   **Other options**
					     
   :f logical output_progress_vectors [default=false]: |output_progress_vectors|

   :f integer save_covm [default=0]: |save_covm|
				     
				     .. include:: ../common/options_save_covm.txt
   
   **Internal options to help solving a regularized problem implicitly**

   :f integer regularization [default=0]: |regularization|
					  
					  .. include:: ../common/options_regularization.txt
					  
   :f real regularization_term [default=0.0]: |regularization_term|
						
   :f real regularization_power [default=0.0]: |regularization_power|



The derived data type for holding information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../common/inform.rst

.. f:type:: nlls_inform
	    
   This is used to hold information about the progress of the algorithm.

   :f integer status: |status|

   :f character error_message(80): |error_message|

   :f integer alloc_status: |alloc_status|

   :f character bad_alloc(80): |bad_alloc|

   :f integer iter: |iter|

   :f integer f_eval: |f_eval|
		      
   :f integer g_eval: |g_eval|

   :f integer h_eval: |h_eval|

   :f intgeer hp_eval: |hp_eval|

   :f integer convergence_normf: |convergence_normf|

   :f integer convergence_normf: |convergence_normg|
				
   :f integer convergence_normf: |convergence_norms|
				   
   :f real resvec(iter+1): |resvec|

   :f real resvec(iter+1): |gradvec|

   :f real obj: |obj|

   :f real norm_g: |norm_g|

   :f real scaled_g: |scaled_g|

   :f integer external_return: |external_return|
   
   :f character external_name(80): |external_name|

   :f real step: |step|

   :f real cov: |cov|

   :f real var: |var|
   

The workspace derived data type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. f:type:: nlls_workspace
	    
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

   
