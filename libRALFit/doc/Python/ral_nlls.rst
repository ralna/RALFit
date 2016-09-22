==================
:py:mod:`ral_nlls`
==================

.. py:module:: ral_nlls

This module is a binary module that provides non-linear least squares solvers.

.. py:function:: (x, info) = solve(x0, r, J, Hr=None, params=None, options=None)

   Returns a solution to the least squares problem

   .. math::
      \min_x F(x):=\frac{1}{2}\|r(x)\|_2^2

   The user must supply an initial guess ``x0``, a callable function
   ``r(x[, params])`` that returns residual vector :math:`r(x)`, and a callable
   ``J(x[, params])`` that returns the :math:`m\times n` Jacobian matrix
   :math:`J(x)=\nabla r(x)`. The user may also supply a callable
   ``Hf(x, r[, params])`` that returns the second order information

   .. math::
      Hf(x) = \sum_{i=1}^m r_i(x) \nabla^2 r_i(x).

   Should the ``params`` parameter not be ``None``, then it is passed unaltered
   to the usable callable routines, and is otherwise not used by the routine.
   If it is ``None``, it is not passed.

   The ``options`` argument is a dictionary specifying options for running the
   routine. Possible options, with their default values, are as follows:

   * ``print_level=0`` integer specifying amount of printing by the routine.

     * ``<= 0`` No output;
     * ``= 1`` One line summary for each iteration;
     * ``= 2`` As above, plus summary of inner iterations; or
     * ``>= 3`` Increasingly verbose (debugging) output.

   * ``maxit=100`` integer specifying maximum number of iterations to perform
   * ``model=1`` integer specifying mathematical model to use

     * ``1`` Gauss-Newton (no Hessian);
     * ``2`` Newton (exact Hessian); or
     * ``9`` Hybrid method (mixture of Gauss-Newton/Newton as appropriate).

   * ``nlls_method=1`` integer specifying solution method for trust-region
     sub-problem:

     * ``1`` Powell's dogleg (approximates the solution);
     * ``2`` Adachi-Iwata-Nakatsukasa-Takeda (AINT);
     * ``3`` More-Sorensen; or
     * ``4`` Galahad DTRS.

   * ``stop_g_absolute=1e-5`` float specifying absolute tolerance for
     convergence.
   * ``stop_g_relative=1e-8`` float specifying relative tolerance for
     convergence.
   * ``relative_tr_radius=0`` integer specifying whether initial trust region
     should be scaled.

Example
=======

The below code fits the model

.. math::
   y_i = x_1 e^{x_2t_i}

to the data:

+---+-----+-----+-----+------+------+
| t | 1.0 | 2.0 | 4.0 |  5.0 |  8.0 |
+---+-----+-----+-----+------+------+
| y | 3.0 | 4.0 | 6.0 | 11.0 | 20.0 |
+---+-----+-----+-----+------+------+

.. literalinclude:: ../../example/Python/solve.py
