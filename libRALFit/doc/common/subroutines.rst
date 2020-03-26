.. |n| replace:: holds the number :math:`n` of variables to be fitted; i.e., :math:`n` is the length of the unknown vector :math:`\bm x`. **Restriction:** **n** :math:`>` **0**.

.. |m| replace:: holds the number :math:`m` of data points available; i.e., :math:`m` is the number of residuals :math:`r_i`. **Restriction:** **m** :math:`\geq` **0**.

.. |X| replace:: on entry, it must hold the initial guess for :math:`\bm x`, and on successful exit it holds the solution to the non-linear least squares problem.

.. |eval_r_desc| replace:: given a point :math:`{\bm x} _{k}^{}`, returns the vector :math:`{\bm r} ({\bm x} _{k}^{})`. Further details of the format required are given in |eval_r| in :ref:`user-routines`.

.. |eval_J_desc| replace:: given a point :math:`{\bm x} _{k}^{}`, returns the :math:`m \times n` Jacobian matrix, :math:`{\bm J} _{k}^{}`, of :math:`{\bm r}` evaluated at :math:`{\bm x} _{k}^{}`. Further details of the format required are given in |eval_J| in :ref:`user-routines`.

.. |eval_Hf_desc| replace:: given vectors :math:`{\bm x} \in \mathbb{R}^n` and :math:`{\bm r} \in \mathbb{R}^m`, returns the quantity :math:`\sum_{i=1}^m ( {\bm r} )_i \nabla^2  {\bm r} _i ( {\bm x} )`. Further details of the format required are given in |eval_Hf| in :ref:`user-routines`.  If ``exact_second_derivative = .false.`` in |nlls_options|, then this is not referenced.

.. |params| replace:: holds parameters to be passed to the user-defined routines |eval_r|, |eval_J|, and |eval_Hf|. Further details of its use are given in :ref:`user-routines`.

.. |options| replace:: controls execution of algorithm.

.. |inform| replace:: components provide information about the execution of the subroutine.

.. |weights| replace:: If present, this holds the square-roots of the diagonal entries of the weighting matrix, :math:`{\bm W}`. If absent, then the norm in the least squares problem is taken to be the 2-norm, that is, :math:`{\bm W} = I`.

.. |eval_HP_desc| replace:: If present, this is a routine that, given vectors :math:`{\bm x}, {\bm y} \in \mathbb{R}^m`, returns the matrix :math:`P({\bm x},{\bm y}) := ( H_1({\bm x}){\bm y} \dots  H_m({\bm x}){\bm y})`. Further details of the format required are given in |eval_HP| in :ref:`user-routines`.  This is only referenced if ``model = 4`` in |nlls_options|.

.. |lower_bounds| replace:: If present, this contains the lower bounds on the solution.

.. |upper_bounds| replace:: If present, this contains the upper bounds on the solution.
			    
.. |iterate_X| replace:: on the first call it must hold the initial guess for :math:`\bm x`. On return it holds the value of :math:`\bm x` at the current iterate,  and must be passed unaltered to any subsequent call to |nlls_iterate|.

.. |w| replace:: is used to store the current state of the iteration and should not be altered by the user.

.. |eval_r_n| replace:: is passed unchanged as provided in the call to |nlls_solve|/|nlls_iterate|.

.. |eval_r_m| replace:: is passed unchanged as provided in the call to |nlls_solve|/|nlls_iterate|.

.. |eval_r_params| replace:: is passed unchanged as provided in the call to |nlls_solve|/|nlls_iterate|.

.. |eval_r_X| replace::  holds the current point :math:`{\bm x}_{k}^{}` at which to evaluate :math:`{\bm r} ( {\bm x} _{k}^{})`.

.. |eval_r_r| replace:: must be set by the routine to hold the residual function evaluated at the current point :math:`{\bm x} _{k}^{}`, :math:`{\bm r} ({\bm x} _{k}^{})`.

.. |eval_r_status| replace:: is initialised to ``0`` before the routine is called. If it is set to a non-zero value by the routine, then |nlls_solve| / |nlls_iterate| will exit with error.

.. |eval_J_n| replace:: is passed unchanged as provided in the call to |nlls_solve|/|nlls_iterate|.

.. |eval_J_m| replace:: is passed unchanged as provided in the call to |nlls_solve|/|nlls_iterate|.

.. |eval_J_params| replace:: is passed unchanged as provided in the call to |nlls_solve|/|nlls_iterate|.

.. |eval_J_X| replace::  holds the current point :math:`{\bm x}_{k}^{}` at which to evaluate :math:`{\bm J} (  {\bm x} _{k}^{})`.

.. |eval_J_r| replace:: must be set by the routine to hold the Jacobian of the residual function evaluated at the current point :math:`{\bm x}_{k}^{}`, :math:`{\bm r} (  {\bm x} _{k}^{})`. ``J(i*m+j)`` must be set to hold :math:`\nabla_{x_j} r_i(  {\bm x} _{k}^{})`.

.. |eval_J_status| replace:: is initialised to ``0`` before the routine is called. If it is set to a non-zero value by the routine, then |nlls_solve| / |nlls_iterate| will exit with error.

.. |eval_Hf_n| replace:: is passed unchanged as provided in the call to |nlls_solve|/|nlls_iterate|.

.. |eval_Hf_m| replace:: is passed unchanged as provided in the call to |nlls_solve|/|nlls_iterate|.

.. |eval_Hf_params| replace:: is passed unchanged as provided in the call to |nlls_solve|/|nlls_iterate|.

.. |eval_Hf_X| replace::  holds the current point :math:`{\bm x}_{k}^{}` at which to evaluate :math:`\sum_{i=1}^m ( {\bm r} )_i \nabla^2 r_i( {\bm x} )`.

.. |eval_Hf_r| replace:: holds :math:`{\bm W}  {\bm r} ( {\bm x} )`, the (weighted) residual, as computed from vector returned by the last call to |eval_r|.

.. |eval_Hf_Hf| replace:: must be set by the routine to hold the matrix :math:`\sum_{i = 1}^m ( {\bm r} )_{i}\nabla^2 r_{i}^{}(  {\bm x} _{k}^{})`, held by columns as a vector, where :math:`( {\bm r} )_i` denotes the :math:`i`\ th component of :math:`\texttt{r}`, the vector passed to the routine.

.. |eval_Hf_status| replace:: is initialised to ``0`` before the routine is called. If it is set to a non-zero value by the routine, then |nlls_solve| / |nlls_iterate| will exit with error.

.. |eval_HP_n| replace:: is passed unchanged as provided in the call to |nlls_solve|/|nlls_iterate|.

.. |eval_HP_m| replace:: is passed unchanged as provided in the call to |nlls_solve|/|nlls_iterate|.

.. |eval_HP_params| replace:: is passed unchanged as provided in the call to |nlls_solve|/|nlls_iterate|.

.. |eval_HP_x| replace::  holds the current point :math:`{\bm x}_{k}^{}` at which to evaluate the Hessians :math:`\nabla^2 r_i( {\bm x_k} )`.

.. |eval_HP_y| replace:: holds :math:`{\bm y}`, the vector which multiplies each Hessian.

.. |eval_HP_HP| replace:: must be set by the routine to hold the matrix :math:`P({\bm x},{\bm y}) := ( H_1({\bm x}){\bm y} \dots  H_m({\bm x}){\bm y})`, held by columns as a vector.

.. |eval_HP_status| replace:: is initialised to ``0`` before the routine is called. If it is set to a non-zero value by the routine, then |nlls_solve| / |nlls_iterate| will exit with error.
