.. |out| replace::  the Fortran unit number for general messages. If it is negative, these messages will be suppressed.

.. |print_level| replace:: controls the level of output required.  Options are:

.. !+----+------------------------------------------------------------------------+
.. !|  0 | No informational output will occur.                                    |
.. !+----+------------------------------------------------------------------------+
.. !|  1 | Prints a brief summary when finished.                                  |
.. !+----+------------------------------------------------------------------------+
.. !|  2 | Gives a one-line summary for each iteration.                           |
.. !+----+------------------------------------------------------------------------+
.. !|  3 | As 2, but with more details.                                           |
.. !+----+------------------------------------------------------------------------+
.. !|  4 | As 3, plus gives a summary of the inner iteration for each iteration.  |
.. !+----+------------------------------------------------------------------------+
.. !|  5 | As 4, plus gives more verbose (debugging) output.                      |
.. !+----+------------------------------------------------------------------------+

.. |print_options| replace:: determines whether to print a list of all options and their values at the beggining of the solve.

.. |print_header| replace:: it will print the column header every 30 iterations when ``print_level > 1``.

.. |maxit| replace:: gives the number of iterations the algorithm is allowed to take before being stopped.  This is not accessed if |nlls_iterate| is used. 

.. |model| replace:: specifies the model, :math:`m_k(\cdot)`, used.  Possible values are:

.. |type_of_method| replace:: specifies the type of globalization method used.  Possible values are:

.. |nlls_method| replace:: specifies the method used to solve (or approximate the solution to) the trust-region sub problem.  Possible values are:

.. |stop_g_absolute| replace:: specifies the absolute tolerance used in the convergence test on :math:`\|{\iter{\vJ}}^T\vr(\iter{\vx}))\|/\|\vr(\iter{\vx})\|`.

.. |stop_g_relative| replace:: specifies the relative tolerance used in the convergence test on :math:`\|{\iter{\vJ}}^T\vr(\iter{\vx})\|/\|\vr(\iter{\vx})\|`.

.. |stop_f_absolute| replace:: specifies the absolute tolerance used in the convergence test on :math:`\|\vr(\iter{\vx})\|`.

.. |stop_f_relative| replace:: specifies the relative tolerance used in the convergence test on :math:`\|\vr(\iter{\vx})\|`.

.. |stop_s| replace:: specifies the tolerance used in the convergence test on :math:`\|\iter{\vs}\|`.

.. |relative_tr_radius| replace:: specifies whether the initial trust region radius should be scaled.

.. |initial_radius_scale| replace:: specifies the scaling parameter for the initial trust region radius, which is only used if ``relative_tr_radius = 1``.

.. |initial_radius| replace:: specifies the initial trust-region radius, :math:`\Delta`.

.. |regularization| replace:: specifies if a regularized non-linear least squares 
			      problem needs to be solved, and if so, what method is 
			      used to solve it.  Possible values are:

.. |regularization_term| replace:: specifies the regularization weight, :math:`\sigma`, used in the least-squares problem.

.. |regularization_power| replace:: specifies the regularization index, :math:`p`, used in the least-squares problem.

.. |maximum_radius| replace:: specifies the maximum size permitted for the trust-region radius.

.. |eta_successful| replace:: specifies the smallest value of :math:`\rho` such that the step is accepted.

.. success_but_reduce is also available, but not documented

.. |eta_very_successful| replace:: specifies the value of :math:`\rho` after which the trust-region radius is increased.

.. |eta_too_successful| replace:: specifies that value of :math:`\rho` after which the step is accepted, but keep the trust-region radius unchanged.

.. |radius_increase| replace:: specifies the factor to increase the trust-region radius by.

.. |radius_reduce| replace:: specifies the factor to decrease the trust-region radius by.

.. |tr_update_strategy| replace:: specifies the strategy used to update :math:`\Delta_k`.  Possible values are:

.. |hybrid_switch| replace:: specifies the value, if ``model=3``, at which second derivatives are used.

.. |exact_second_derivatives| replace:: if ``true``, signifies that the exact second derivatives are available (and, if ``false``, approximates them using a secant method).

.. |scale| replace:: specifies how, if at all, we scale the Jacobian.  We calculate a diagonal scaling matrix, :math:`{\tt D}`, as follows: 

.. |scale_trim_max| replace:: specifies whether or not to trim large values of the scaling matrix, :math:`D`. If ``true``, :math:`{\tt D}_{i,i} \leftarrow min({\tt D}_{i,i}, {\tt scale\_max})`.

.. |scale_max| replace:: specifies the maximum value allowed if ``scale_trim_max = true``.

.. |scale_trim_min| replace:: specifies whether or not to trim small values of the scaling matrix, :math:`{\tt D}`. If ``true``, :math:`{\tt D}_{i,i} \leftarrow max({\tt D}_{i,i}, {\tt scale_max})`.

.. |scale_min| replace:: specifies the minimum value allowed if ``scale_trim_max = true``.

.. |scale_require_increase| replace:: specifies whether or not to require :math:`{\tt D}_{i,i}` to increase before updating it.

.. |calculate_svd_J| replace:: specifies whether or not to calculate the singular value decomposition of :math:`{\tt J}` at each iteration.  

.. |more_sorensen_maxits| replace:: if ``nlls_method = 3``, specifies the maximum number of iterations allowed in the More-Sorensen method.

.. |more_sorensen_shift| replace:: if ``nlls_method = 3``, specifies the shift to be used in the More-Sorensen method. 

.. |more_sorensen_tiny| replace:: if ``nlls_method = 3``, specifies the value below which numbers are considered to be essentially zero.

.. |more_sorensen_tol| replace:: if ``nlls_method = 3``, specifies the tolerance to be used in the More-Sorensen method.

.. |hybrid_tol| replace:: if ``model=3``, specifies the value such that if  :math:`\|{\iter{\vJ}}^T \vW \vr(\vx_k) \|_2 < \mathtt{hybrid\_tol} * 0.5 \|\vr(\vx_k)\|_\vW^2` the method switches to a quasi-Newton method.

.. |hybrid_switch_its| replace:: if ``model=3``, sets how many iterates in a row must the condition in the definition of ``hybrid_tol`` hold before a switch.

.. |reg_order| replace:: if ``nlls_method = 4``, the order of the regularization used (:math:`p` in **TODO** (eq:: reg_subproblem)).   If ``reg_order = 0.0``, then the algorithm chooses an appropriate value of :math:`p`. 

.. |inner_method| replace::  if ``nlls_method = 4``, specifies the method used to pass in the regularization parameter to the inner non-linear least squares solver.   Possible values are:

.. |output_progress_vectors| replace:: if true, outputs the progress vectors ``nlls_inform%resvec`` and ``nlls_inform%gradvec`` at the end of the routine.
