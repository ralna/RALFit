.. |status| replace:: gives the exit status of the subroutine.  See :ref:`errors` for details.

.. |error_message| replace:: holds the error message corresponding to the exit status.

.. |alloc_status| replace:: gives the status of the last attempted allocation/deallocation.

.. |bad_alloc| replace:: holds the name of the array that was being allocated when an error was flagged.

.. |iter| replace:: gives the total number of iterations performed.

.. |f_eval| replace:: gives the total number of evaluations of the objective function.

.. |g_eval| replace:: gives the total number of evaluations of the gradient of the objective function.

.. |h_eval| replace:: gives the total number of evaluations of the Hessian of the objective function using ``eval_hf``.

.. |hp_eval| replace:: gives the total number of evaluations of the Hessian of the objective function using ``eval_hp``.

.. |fd_f_eval| replace:: gives the total number of evaluations of the objective function due to approximating or checking the gradient of the objective function.

.. |convergence_normf| replace:: tells us if the test on the size of :math:`\vr` is satisfied.

.. |convergence_normg| replace:: that tells us if the test on the size of the gradient is satisfied.

.. |convergence_norms| replace:: that tells us if the test on the step length is satisfied.

.. |resvec| replace:: if ``output_progress_vectors=true`` in |nlls_options|, holds the vector of residuals.

.. |gradvec| replace:: if ``output_progress_vectors=true`` in |nlls_options|, holds the vector of gradients.

.. |obj| replace:: holds the value of the objective function at the best estimate of the solution determined by the algorithm.

.. |norm_g| replace:: holds the gradient of the objective function at the best estimate of the solution determined by the package.

.. |scaled_g| replace:: holds a scaled version of the gradient of the objective function at the best estimate of the solution determined by the package.

.. |external_return| replace:: gives the error code that was returned by a call to an external routine.

.. |external_name| replace:: holds the name of the external code that flagged an error.

.. |step| replace:: holds the size of the last step taken.

.. |cov| replace:: On exit, optionally contains information about the covariance matrix, as requested by ``nlls_options%save_covm``.

.. |var| replace::  On exit, optionally contains information about the diagonal of the covariance matrix, as requested by ``nlls_options%save_covm``.
