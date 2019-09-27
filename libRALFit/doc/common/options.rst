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

.. |box_nFref_max| replace:: Memory size for the non-monotone linesearch.

.. |box_gamma| replace:: Kanzow sufficient decrease ratio (see equtaion 25 in Kanzow [2004])

.. |box_decmin| replace:: FIXME! macheps

.. |box_bigbnd| replace:: Magic number to consider box bound (+/-) infinity

.. |box_wolfe_descent| replace:: Wolfe descent condition :math:`0<\sigma_1<1/2`

.. |box_wolfe_curvature| replace:: Wolfe curvature condition :math:`0<\sigma_2`

.. |box_kanzow_power| replace:: Tolerance to consider projected direction a descent direction.  See LS STEP, Section 4, p392, Kanzow [2014]

.. |box_kanzow_descent| replace:: FIXME! sqrt(mcheps)

.. |box_quad_model_descent| replace:: FIXME! sqrt(mcheps)

.. |box_tr_test_step| replace:: Take a projected trust region step when trust region test is OK?  If false, force a linesearch or projected gradient step.

.. |box_wolfe_test_step| replace:: Take a project trust region step when Wolfe test is OK?  If false, force a linesearch or projected gradient step.

.. |box_tau_min| replace:: Threshold to determine if the projection of the trust region direction is too severe (:math:0 < \tau_min 1)

.. |box_tau_descent| replace:: :math:`\tau >=` ``tau_descent`` in order to test for descent

.. |box_max_ntrfail| replace:: Max times trust region iterations can fail without passing the various descent tests. Ignored when :math:`\mathrm{proj}(x)==x`.

.. |box_quad_match| replace:: Number of consecutive times quadratic model matches :math:`f(x_{k+1})` required before setting initial alpha step for projected gradient step equal to ``scale_alpha*alpha_k-1`` (FIXME)

.. |box_alpha_scale| replace:: Initial step scale (if :math:`\mathrm{quad}_i >= \mathrm{box\_quad}_i`)

.. |box_Delta_scale| replace:: Scaling factor to use when updating Delta from linesearch/project gradient step

.. |box_tau_wolfe| replace:: FIXME

.. |box_tau_tr_wolfe| replace:: FIXME

.. |box_tau_tr_step| replace:: FIXME

.. |box_ls_step_maxit| replace:: FIXME

.. |box_linesearch_type| replace:: Linesearch type: 1 is Dennis-Schnabel, 2 is Hager-Zhang (FIXME - link to description)

.. |hybrid_tol| replace:: if ``model=3``, specifies the value such that if  :math:`\|{\iter{\vJ}}^T \vW \vr(\vx_k) \|_2 < \mathtt{hybrid\_tol} * 0.5 \|\vr(\vx_k)\|_\vW^2` the method switches to a quasi-Newton method.

.. |hybrid_switch_its| replace:: if ``model=3``, sets how many iterates in a row must the condition in the definition of ``hybrid_tol`` hold before a switch.

.. |reg_order| replace:: if ``type_of_method = 2``, the order of the regularization used (:math:`p` in  (:eq:`regsub`)).   If ``reg_order = 0.0``, then the algorithm chooses an appropriate value of :math:`p`. 

.. |inner_method| replace::  if ``nlls_method = 4``, specifies the method used to pass in the regularization parameter to the inner non-linear least squares solver.   Possible values are:

.. |output_progress_vectors| replace:: if true, outputs the progress vectors ``nlls_inform%resvec`` and ``nlls_inform%gradvec`` at the end of the routine.
