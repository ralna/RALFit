**RALFit** computes a solution :math:`\vx` to the non-linear least-squares problem

.. math::

   \min_\vx \  F(\vx) := \frac{1}{2}\| \vr(\vx) \|_{\vW}^2 + \frac{\sigma}{p}\| \vx\|_2^p,


where :math:`\vW\in\mathbb{R}^{m\times m}` is a diagonal, non-negative, weighting matrix, 
and :math:`\vr(\vx) =(\comp[1]{r}(\vx), \comp[2]{r}(\vx),...,\comp[m]{r}(\vx))^T` 
is a non-linear function.

A typical use may be to fit a function :math:`f(\vx,t)` to the data :math:`y_i, \ t_i`, 
weighted by the uncertainty of the data, :math:`\sigma_i`:

.. math:: 

   \min_\vx \  \frac{1}{2} \sum_{i=1}^m \left(\frac{y_i - f(\vx;t_i)}{\sigma_i}\right)^2,

which corresponds to taking :math:`r_i(\vx) := y_i - f(\vx;t_i)` and
:math:`\vW` such that :math:`\vW_{ii} = (1/{\sigma_i^2}).`
For this reason we refer to the function :math:`\vr` as the *residual* function.

Various algorithms for solving this problem are implemented -- see :doc:`method`.


