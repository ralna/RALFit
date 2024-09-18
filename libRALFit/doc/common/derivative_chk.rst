Writing the derivative call-back function can be error-prone, and incorrect derivatives are
one of the most frequent reasons for lack of convergence.
To address this, setting the option ``check_derivatives`` activates a checker for the
the call-back that evaluates :math:`J=\nabla r(x)`.
The checker produces a table similar to

.. code::

   Begin Derivative Checker

       Jacobian storage scheme (Fortran_Jacobian) = Fortran (column-major)

       Jac[     1,     1] =   2.013752707470E+00 ~   0.000000000000E+00  [ 2.014E+04], ( 0.224E-06)  X  Skip
       Jac[     1,     2] =   4.055199966845E+00 ~   0.000000000000E+00  [ 4.055E+04], ( 0.224E-06)  XT  Skip
       Jac[     1,     3] =   3.311545195869E+01 ~   0.000000000000E+00  [ 3.312E+05], ( 0.224E-06)  XT  Skip
       Jac[     2,     1] =   3.020629061206E+00 ~   3.020628836751E+00  [ 7.431E-08], ( 0.149E-06)
       ...
       Jac[    20,     3] =   9.866788062658E+01 ~   9.866785123348E+01  [ 2.979E-07], ( 0.149E-06)

       Warning: derivative checker skipped      5 entries that have too tight bounds on the variable(s).

       Note: derivative checker detected that     66 entries may correspond to the transpose.

   End Derivative Checker


The initial line indicates the storage option set for the Jacobian matrix.
The first column after the equal sign (``=``), is the derivative returned by the user-supplied call-back, the column
after the ``~`` sign is the approximated finite-difference derivative, the value inside the brackets is the relative
threshold
:math:`\frac{|\mathrm{approx} - \mathrm{exact}|}{\max(|\mathrm{approx}|,\; \mathrm{fd_{ttol}})}`,
(option :math:`\mathrm{fd_{ttol}}` is defined by the option ``derivative_test_tol``).
The value inside the parenthesis is
the relative tolerance to compare the relative threshold against. The last column provides some flags: ``X`` to
indicate that the threshold is larger than the tolerance and is deemed likely to be wrong. ``T`` indicates that
the value stored in :math:`J(i,j)` corresponds the to the value belonging to the transposed Jacobian matrix,
providing a hint that possibly the storage sequence is likely wrong. It also hints to check that the matrix is
being stored row-major and that the solver option ``Fortran_Jacobian`` is set to column-major or vice-versa.
Finally, ``Skip`` indicates that either the associated variable is fixed (constrained to a fixed value) or the
bounds on it are too tight to perform a finite-difference approximation and thus the check for this entry cannot be
checked and is skipped.

The derivative checker uses finite-differences to compare with the user provided derivatives and such the quality of
the approximation depends on the finite-difference step used (see option ``finite_difference_step``).

The option ``derivative_test_tol`` is involved in defining the relative tolerance to decide if the user-supplied
derivative is correct, a smaller value implies a more stringent test.

Under certain circumstances the checker may signal false-positives, tweaking the options ``finite_difference_step``
and ``derivative_test_tol`` guard against this happening.

.. Note::

   It is highly recommended that the option
   ``check_derivatives`` is activated during the writing or development of the derivative call-back.
   After validating the Jacobian matrix the option should be reset to avoid performance impact.

