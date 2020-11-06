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

.. |print_header| replace:: prints the column header every ``print_header`` iterations when ``print_level > 1``.

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

.. |regularization| replace:: specifies the method by which a regularized non-linear least squares 
			      problem is solved implicitly.  Is designed to be used when solving the
			      nonlinear least-squares problem recursively.  Possible values are:

.. |regularization_term| replace:: specifies the regularization weight, :math:`\sigma`, used when implicitly solving the least-squares problem.

.. |regularization_power| replace:: specifies the regularization index, :math:`p`, used when implicitly solving  the least-squares problem.

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

.. |more_sorensen_maxits| replace:: if ``nlls_method = 3``, specifies the maximum number of iterations allowed in the More-Sorensen method.

.. |more_sorensen_shift| replace:: if ``nlls_method = 3``, specifies the shift to be used in the More-Sorensen method. 

.. |more_sorensen_tiny| replace:: if ``nlls_method = 3``, specifies the value below which numbers are considered to be essentially zero.

.. |more_sorensen_tol| replace:: if ``nlls_method = 3``, specifies the tolerance to be used in the More-Sorensen method.

.. |box_nFref_max| replace:: Memory size for the non-monotone projected gradient linesearch.

.. |box_gamma| replace:: Sufficient decrease parameter (:math:`0 < \gamma < 1`). A step is deemed successful if :math:`F(x_{k+1}) \leq \gamma F(x_k)`.

.. |box_decmin| replace:: Defines a safe :math:`\epsilon_{\rm machine}`.
			  Reserved for compatibility use.

.. |box_bigbnd| replace:: Value used as a proxy for :math:`\pm \infty` in the box bound.

.. |box_wolfe_descent| replace:: Wolfe descent condition parameter :math:`0<\sigma_1<1/2`.

.. |box_wolfe_curvature| replace:: Wolfe curvature condition parameter :math:`0<\sigma_2<1`.

.. |box_kanzow_power| replace:: Parameter setting :math:`\nu > 1` in equation :eq:`ls-eqn`. 

.. |box_kanzow_descent| replace:: Parameter setting :math:`\kappa > 0` in equation :eq:`ls-eqn` See (LS STEP, Section 4, page 392, Kanzow 2014).  The descent test in equation :eq:`ls-eqn` is only accepted if the projected trust region step was not too severe, specifically, when the projection ratio :math:`\tau` is greater than :math:`\tau_d`.

.. |box_quad_model_descent| replace:: Error tolerance, :math:`\epsilon_q > 0`, required for
   the error :math:`e_k := |f(x_{k+1} - q(x_{k+1})|` in a quadratic model,
   :math:`q(x_{k+1})` of :math:`F(x_{k+1})` to be deemed sufficient.

.. |box_tr_test_step| replace::  If true, then a trust region step is taken if the trust region loop is successful and projection ratio :math:`\tau` is about :math:`\tau_{\rm TR}`.

.. |box_wolfe_test_step| replace:: If true, and the current point satisfies the weak Wolfe descent conditions, and the projection ratio :math:`\tau` is above :math:`\tau_W`, then the step is taken.

.. |box_tau_min| replace:: Threshold, :math:`0 < \tau_{min} < \min(\tau_W,\tau_{TR})`, to determine if the projection of the trust region direction is too severe.  If this is the case, the trust region loop is terminated and a linesearch step is forced.

.. |box_tau_descent| replace:: Tolerance :math:`0 < \tau_d < 1` to test if the projected trust region step descends. 

.. |box_max_ntrfail| replace:: Number of unsuccessful trust region iterations to allow without passing the various descent tests. Ignored when :math:`\mathrm{proj}(x)=x`.

.. |box_quad_match| replace:: Number of consecutive times :math:`e_k < \epsilon_q` required before setting ``\alpha_0 = \mu_\alpha * alpha_k-1`` in the projected gradient step.

.. |box_alpha_scale| replace:: Scale factor :math:`\mu_\alpha` for initial step length in the projected gradient linesearch.
			       
.. |box_Delta_scale| replace:: Scaling factor :math:`\mu_\Delta`, used when updating :math:`\Delta_{k+1} = \mu_\Delta\|s_k\|` from linesearch/project gradient step.

.. |box_tau_wolfe| replace:: Tolerance that defines the value of :math:`0 < \tau_W < 1`.

.. |box_tau_tr_step| replace:: Tolerance :math:`0 < \tau_{TR} < 1` to allow a successful trust region step.

.. |box_ls_step_maxit| replace:: Maximum number of iterations to perform in the linesearch step.

.. |box_linesearch_type| replace:: Linesearch type -- available options are:

.. |hybrid_tol| replace:: if ``model=3``, specifies the value such that if  :math:`\|{\iter{\vJ}}^T \vW \vr(\vx_k) \|_2 < \mathtt{hybrid\_tol} * 0.5 \|\vr(\vx_k)\|_\vW^2` the method switches to a quasi-Newton method.

.. |hybrid_switch_its| replace:: if ``model=3``, sets how many iterates in a row must the condition in the definition of ``hybrid_tol`` hold before a switch.

.. |reg_order| replace:: if ``type_of_method = 2``, the order of the regularization used (:math:`q` in  (:eq:`regsub`)).   If ``reg_order = 0.0``, then the algorithm chooses an appropriate value of :math:`q`. 

.. |inner_method| replace::  if ``nlls_method = 4``, specifies the method used to pass in the regularization parameter to the inner non-linear least squares solver.   Possible values are:

.. |output_progress_vectors| replace:: if true, outputs the progress vectors ``nlls_inform%resvec`` and ``nlls_inform%gradvec`` at the end of the routine.

.. |save_covm| replace:: Determines whether to return information about the covariance matrix in ``nlls_informm``.  Options are:

