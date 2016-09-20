Summary
=======

**RALFit** computes a solution :math:`\vx` to the non-linear least-squares problem

.. math::

   \min_\vx \  F(\vx) := \frac{1}{2}\| \vr(\vx) \|_{\vW}^2 + \frac{\sigma}{p}\| \vx\|_2^p,


where :math:`\vW\in\mathbb{R}^{m\times m}` is a diagonal, non-negative, weighting matrix, 
and :math:`\vr(\vx) =(\comp[1]{r}(\vx), \comp[2]{r}(\vx),...,\comp[m]{r}(\vx))^T` 
is a non-linear function.

A typical use may be to fit a function :math:`f(\vx)` to the data :math:`y_i, \ t_i`, 
weighted by the uncertainty of the data, :math:`\sigma_i`, so that

.. math::
   
   r_i(\vx) := y_i - f(\vx;t_i),

and :math:`\vW` is the diagonal matrix such that 
:math:`\vW_{ii} = (1/\sqrt{\sigma_i}).`
For this reason we refer to the function :math:`\vr` as the *residual* function.

The algorithm is iterative.
At each point, :math:`\iter{\vx}`, the algorithm builds a model of the function at the next step, :math:`F({\iter{\vx}+\iter{\vs}})`, which we refer to as :math:`m_k(\cdot)`.  We allow either a Gauss-Newton model, a (quasi-)Newton model, or a Newton-tensor model; 
see Section **TODO** \ref{sec:model_description} 
for more details about these models.  

Once the model has been formed we find a candidate for the next step by either solving a trust-region sub-problem of the form

.. math:: 

  \iter{\vs} = \arg \min_{\vs} \ \iter{m} (\vs) \quad \mathrm{s.t.} \quad  \|\vs\|_B \leq \Delta_k,

or by solving the regularized problem 

.. math::

   \iter{\vs} = \arg \min_{\vs} \ \iter{m} (\vs)  + \frac{1}{\Delta_k}\cdot \frac{1}{p} \|\vs\|_B^p,

where :math:`\Delta_k` is a parameter of the algorithm 
(the trust region radius or the inverse of the regularization parameter respectively), 
:math:`p` is a given integer, and 
:math:`B` is a symmetric positive definite weighting matrix that is calculated by the algorithm.
The quantity

.. math::
   
   \rho = \frac{F(\iter{\vx}) - F(\iter{\vx} + \iter{\vs})}{\iter{m}(\iter{\vx}) - \iter{m}(\iter{\vx} + \iter{\vs})}

is then calculated.
If this is sufficiently large we accept the step, and 
:math:`\iter[k+1]{\vx}` is set to 
:math:`\iter{\vx} + \iter{\vs}`; if not, the parameter 
:math:`\Delta_k` is reduced and  the resulting new trust-region sub-problem is solved.  
If the step is very successful -- in that 
:math:`\rho` is close to one --
:math:`\Delta_k` is increased.

This process continues until either the residual, 
:math:`\|\vr(\iter{\vx})\|_\vW`, or a measure of the gradient,
:math:`\|{\iter{\vJ}}^T\vW\vr(\iter{\vx})\|_2 / \|\vr(\iter{\vx})\|_\vW`, 
is sufficiently small.
