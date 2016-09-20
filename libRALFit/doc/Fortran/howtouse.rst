.. include:: ../common/options.rst

.. default-domain:: Fortran

How to use the package
======================

Overview
--------

Calling sequences
^^^^^^^^^^^^^^^^^

Access to the package requires a ``USE`` statement

.. code-block:: Fortran

   use ral_nlls_double


The user can then call one of the procedures:

:f:subr:`nlls_solve()` solves the non-linear least squares problem.
(**TODO**: link)

:f:subr:`nlls_iterate()` performs one iteration of the non-linear least squares
solver.


The derived data types
^^^^^^^^^^^^^^^^^^^^^^

For each problem, the user must employ the derived types
defined by the module to declare scalars of the types :f:type:`nlls_inform` and :f:type:`nlls_options`. 
If :f:subr:`nlls_iterate()` is to be used, then a scalar of the type :f:type:`nlls_workspace` must also be defined. 
The following pseudocode illustrates this.

.. code-block:: Fortran

   use nlls_module
   ...
   type (NLLS_inform) :: inform
   type (NLLS_options) :: options
   type (NLLS_workspace) :: work ! needed if nlls_iterate to be called
   ...


The components of :f:type:`nlls_options` and 
:f:type:`nlls_inform` are explained in
Sections **TODO**

Argument lists and calling sequences
------------------------------------

We use square brackets to indicate  arguments. In each call, optional
arguments follow the argument inform. Since we reserve the right to add
additional optional arguments in future releases of the code, **we
strongly recommend that all optional arguments be called by keyword, not
by position**.

The term **package type** is used to mean default real if the single
precision version is being used and double precision real for the double
precision version.

To solve the non-linear least squares problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. f:subroutine:: nlls_solve(n,m,X,eval_r,eval_J,eval_Hf,params,options,inform[,weights])

   Solves the non-linear least squares problem.
   
   :p integer n [in]: holds the number :math:`n` of
		      variables to be fitted; i.e., :math:`n` is the length of the unknown
		      vector :math:`\bm x`. 
		      **Restriction:** **n** :math:`>` **0**.
   :p integer m [in]: holds the number :math:`m` of data points available; i.e., 
		      :math:`m` is the number of residuals :math:`r_i`. 
		      **Restriction:** **m** :math:`\geq` **0**.
   :p real X(n) [inout]: on entry, it must hold the initial guess for 
			   :math:`\bm x`, and on successful exit it holds the 
			   solution to the non-linear least squares problem.
   :p procedure eval_r: given a point :math:`{\bm x} _{k}^{}`, 
		 returns the vector :math:`{\bm r} ({\bm x} _{k}^{})`. 
		 Further details of the format required are given in 
		 Section **TODO**.
   :p procedure eval_J: given a point :math:`{\bm x} _{k}^{}`, 
			returns the :math:`m \times n` Jacobian matrix, 
			:math:`{\bm J} _{k}^{}`, of :math:`{\bm r}` evaluated at 
			:math:`{\bm x} _{k}^{}`. 
			Further details of the format required are given in
			Section **TODO**.
   :p procedure eval_Hf: given vectors :math:`{\bm x} \in \mathbb{R}^n` and
			 :math:`{\bm r} \in \mathbb{R}^m`, 
			 returns the quantity
			 :math:`\sum_{i=1}^m ( {\bm r} )_i \nabla^2  {\bm r} _i ( {\bm x} )`.
			 Further details of the format required are given in
			 Section **TODO**. 
			 If ``nlls_options%exact_second_derivative = .false.``, 
			 then this is not referenced.

   :p params_base_type params [in]: holds parameters to be
				    passed to the user-defined routines 
				    :f:subr:`eval_r()`, 
				    :f:subr:`eval_J()`, 
				    and :f:subr:`eval_Hf()`.
				    Further details of its use are given in
				    Section **todo**

   :p nlls_options options [in]: controls execution of algorithm.

   :p nlls_inform inform [out]:  components provide information
				 about the execution of the subroutine, as explained in
				 Section **todo**.

   :o real weights(n):  If present,
			this holds the square-roots of the diagonal entries of the weighting
			matrix, :math:`{\bm W}`, in **todo:ref**. 
			If absent, then the norm in **todo::ref** is taken to be the 2-norm,
			that is, :math:`{\bm W} = I`.

To iterate once
^^^^^^^^^^^^^^^


.. f:subroutine:: nlls_iterate(n,m,X,eval_r,eval_J,eval_Hf,params,options,inform[,weights])
		  
   A call of this form allows the user to step through the solution process one
   iteration at a time.

   **n**, **m**, **eval_F**, **eval_J**, **eval_HF**, **params**, **info**,
   **options** and **weights** are as in the desciption of :f:subr:`nlls_solve`.

   :p real X(n) [inout]: on the first call it must hold the initial guess for 
			 :math:`\bm x`. On return it holds the value of 
			 :math:`\bm x` at the current iterate, 
			 and must be passed unaltered to any subsequent 
			 call to :f:subr:`nlls_iterate()`.

   :p nlls_workspace w [inout]: is used to store the
				current state of the iteration and should not be altered by the
				user.

The user may use the components ``info%convergence_normf`` and
``info%convergence_normg`` to determine whether the iteration has
converged.

User-supplied function evaluation routines
------------------------------------------

To evaluate the residual, Jacobian and Hessian at a point, the user must
supply subroutines that perform this operation that **RALFit** will call internally.

In order to pass user-defined data into the evaluation calls, :f:type::`params_base_type` is extended to a :f:type::`user_type`, as follows:

.. code-block:: Fortran

       type, extends( params_base_type ) :: user_type
          ! code declaring components of user_type
       end type user_type

We recommend this type is wrapped in a module with the user-defined
routines for evaluating the function, Jacobian, and Hessian.

The components of the extended type are accessed through a
``select type`` construct:

.. code-block:: Fortran

       select type(params)
       type is(user_type)
         ! code that accesses components of params that were defined within user_type
       end select

For evaluating the function :math:`{\bm r} ( {\bm x} )`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A subroutine must be supplied to calculate :math:`{\bm r} ( {\bm x} )`
for a given vector :math:`{\bm x}`. It must implement the following
interface:

.. code-block:: Fortran

    abstract interface
       subroutine eval_r(n, m, params, x, r, status)
          integer, intent(in) :: n
          integer, intent(in) :: m
          class(params_base_type), intent(in) :: params
          double precision, dimension(n), intent(in) :: x
          double precision, dimension(m), intent(out) :: r
          integer, intent(inout) :: status
       end subroutine eval_r
    end interface

.. f:subroutine:: eval_r(n,m,params,x,r,status)
   
   :p integer n [in]: is passed unchanged as provided in the call to
		  :f:subr:`nlls_solve()`/:f:subr:`nlls_iterate()`.
   :p integer m [in]: is passed unchanged as provided in the call to
		  :f:subr:`nlls_solve()`/:f:subr:`nlls_iterate()`.
   :p params_base_type params [in]: is passed unchanged as provided in the call to
      :f:subr:`nlls_solve()`/:f:subr:`nlls_iterate()`.
   :p real X(n) [in]: holds the current point :math:`{\bm x}_{k}^{}` at which to evaluate
      :math:`{\bm r} ( {\bm x} _{k}^{})`.
   :p real r(m) [out]: must be set by the routine to hold the residual function evaluated at
      the current point :math:`{\bm x} _{k}^{}`, :math:`{\bm r} ({\bm x} _{k}^{})`.
   :p integer status [inout]: is initialised to ``0`` before the routine is called. 
			   If it is set to a non-zero value by the routine, 
			   then :f:subr:`nlls_solve()` / :f:subr:`nlls_iterate()`
			   will exit with error.

For evaluating the function :math:`{\bm J} = \nabla  {\bm r} ( {\bm x} )`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A subroutine must be supplied to calculate
:math:`{\bm J} = \nabla  {\bm r} ( {\bm x} )` for a given vector
:math:`{\bm x}`. It must implement the following interface:

.. code-block:: Fortran

    abstract interface
       subroutine eval_J(n, m, params, x, J, status)
          integer, intent(in) :: n
          integer, intent(in) :: m
          class(params_base_type), intent(in) :: params
          double precision, dimension(n), intent(in)  :: x
          double precision, dimension(n*m), intent(out) :: J
          integer, intent(inout) :: status
      end subroutine eval_J
    end interface

.. f:subroutine:: eval_J(n,m,params,x,J,status)
   
   :p integer n [in]: is passed unchanged as provided in the call to
		  :f:subr:`nlls_solve()`/:f:subr:`nlls_iterate()`.
   :p integer m [in]: is passed unchanged as provided in the call to
		  :f:subr:`nlls_solve()`/:f:subr:`nlls_iterate()`.
   :p params_base_type params [in]: is passed unchanged as provided in the call to
      :f:subr:`nlls_solve()`/:f:subr:`nlls_iterate()`.
   :p real X(n) [in]: holds the current point :math:`{\bm x}_{k}^{}` at which to evaluate
      :math:`{\bm J} (  {\bm x} _{k}^{})`.
   :p real J(m*n) [out]: must be set by the routine to 
			 hold the Jacobian of the residual function
			 evaluated at the current point :math:`{\bm x}_{k}^{}`,
			 :math:`{\bm r} (  {\bm x} _{k}^{})`. 
			 ``J(i*m+j)`` must be set to hold 
			 :math:`\nabla_{x_j} r_i(  {\bm x} _{k}^{})`.
   :p integer status [inout]: is initialised to ``0`` before the routine is called. 
			   If it is set to a non-zero value by the routine, 
			   then :f:subr:`nlls_solve()`/:f:subr:`nlls_iterate()`
			   will exit with error.


For evaluating the function :math:`Hf = \sum_{i=1}^m r_i( {\bm x} )  {\bm W} \nabla^2 r_i( {\bm x} )`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A subroutine must be supplied to calculate
:math:`Hf = \sum_{i=1}^m ( {\bm r} )_i \nabla^2 r_i( {\bm x} )` for
given vectors :math:`{\bm x} \in \mathbb{R}^n` and
:math:`{\bm r} \in \mathbb{R}^m`; here :math:`( {\bm r} )_i` denotes
the :math:`i`\ th component of the vector :math:`{\bm r}`. The
subroutine must implement the following interface:

.. code-block:: Fortran

    abstract interface
       subroutine eval_Hf_type(n, m, params, x, r, Hf, status)
           integer, intent(in) :: n
           integer, intent(in) :: m
           class(params_base_type), intent(in) :: params
           double precision, dimension(n), intent(in)  :: x
           double precision, dimension(m), intent(in)  :: r
           double precision, dimension(n*n), intent(out) :: Hf
           integer, intent(inout) :: status
         end subroutine eval_Hf_type
    end interface
    :language: Fortran

.. f:subroutine:: eval_Hf(n,m,params,x,r,Hf,status)
   
   :p integer n [in]: is passed unchanged as provided in the call to
		  :f:subr:`nlls_solve()`/:f:subr:`nlls_iterate()`.
   :p integer m [in]: is passed unchanged as provided in the call to
		  :f:subr:`nlls_solve()`/:f:subr:`nlls_iterate()`.
   :p params_base_type params [in]: is passed unchanged as provided in the call to
      :f:subr:`nlls_solve()`/:f:subr:`nlls_iterate()`.
   :p real X(n) [in]: holds the current point :math:`{\bm x}_{k}^{}` at which to evaluate
      :math:`\sum_{i=1}^m ( {\bm r} )_i \nabla^2 r_i( {\bm x} )`.
   :p real r(m) [in]: holds :math:`{\bm W}  {\bm r} ( {\bm x} )`, the (weighted) 
		      residual, as computed from vector returned by the last call 
		      to :f:subr::`eval_r()`.
   :p real Hf(n*n) [out]: must be set by the routine to hold the matrix
      :math:`\sum_{i = 1}^m ( {\bm r} )_{i}\nabla^2 r_{i}^{}(  {\bm x} _{k}^{})`,
      held by columns as a vector, where :math:`( {\bm r} )_i`
      denotes the :math:`i`\ th component of :math:`\texttt{r}`, 
      the vector passed to the routine.
   :p integer status [inout]: is initialised to ``0`` before the routine is called. 
			   If it is set to a non-zero value by the routine, 
			   then :f:subr:`nlls_solve()`/:f:subr:`nlls_iterate()`
			   will exit with error.

Data types
----------

The derived data type for holding options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. f:type:: nlls_options
	    
   This is used to hold controlling data. The components  are automatically given default values in the definition of the type.
   
   **Printing Controls**
   
   :f integer error [default=6]: |error|
   :f integer out [default=6]: |out|
   :f integer print_level [default=0]: |print_level|

				       .. include:: ../common/options_print_level.txt
      
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
	 
   :f real regularization_weight [default=0.0]: |regularization_weight|
						
   :f real regularization_power [default=0.0]: |regularization_power|

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
					     
   :f real reg_order [default=0.0]: |reg_order|
		
   :f integer inner_method [default=2]: |inner_method| 

					.. include:: ../common/options_inner_method.txt

   **TODO** check this!

   **More-Sorensen options**  These options are used if ``nlls_method=3``
   
   :f integer more_sorensen_maxits [default=500]: |more_sorensen_maxits|

   :f integer more_sorensen_maxits [default=3]: |more_sorensen_maxits|

   :f real more_sorensen_shift [default=1e-13]: |more_sorensen_shift|

   :f real more_sorensen_tiny [default=10.0*eps]: |more_sorensen_tiny|

   :f real more_sorensen_tol [default=1e-3]: |more_sorensen_tol|
						  
   **Other options**
					     
   :f logical calculate_svd_J [default=false]: |calculate_svd_J|

   :f logical output_progress_vectors [default=false]: |output_progress_vectors|
   


The derived data type for holding information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. f:type:: nlls_inform
	    
   This is used to hold information about the progress of the algorithm.


Warning and error messages
--------------------------


   
