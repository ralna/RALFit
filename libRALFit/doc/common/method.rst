**RALFit** computes a solution :math:`\vx` to the non-linear least-squares problem

.. math::	   
   :label: lsq

   \min_\vx \  F(\vx) := \frac{1}{2}\| \vr(\vx) \|_{\vW}^2 + \frac{\sigma}{p}\| \vx\|_2^p,


Here we describe the method used to solve :eq:`lsq`. **RALFit** implements an iterative method that, at each iteration, calculates and returns a step :math:`\vs` that reduces the model by an acceptable amount by solving (or approximating a solution to) either the trust-region subproblem

.. math:: 
   :label: trsub
  
   \iter{\vs} = \arg \min_{\vs} \ \iter{m} (\vs) \quad \mathrm{s.t.} \quad  \|\vs\|_B \leq \Delta_k,
   
or a regularized problem 

.. math:: 
   :label: regsub
      
   \iter{\vs} = \arg \min_{\vs} \ \iter{m} (\vs)  + \frac{1}{\Delta_k}\cdot \frac{1}{p} \|\vs\|_B^p,

The algorithm is iterative.
At each point, :math:`\iter{\vx}`, the algorithm builds a model of the function at the next step, :math:`F({\iter{\vx}+\iter{\vs}})`, which we refer to as :math:`m_k(\cdot)`.  We allow either a Gauss-Newton model, a (quasi-)Newton model, or a Newton-tensor model; 
see :ref:`models` for more details.

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



Incorporating the regularization term
-------------------------------------

If a non-zero regularization term is required in
:eq:`lsq`, then this is handled by transforming the
problem internally into a least squares problem. The formulation used
will depend on the value of :math:`p`.

**If** :math:`\bf p = 2`,
we solve a least squares problem with
:math:`n` additional degrees of freedom. The new function,
:math:`\widehat{\vr} : \mathbb{R}^{n}\rightarrow\mathbb{R}^{m+n}`, takes
:math:`\widehat{\vr}_i(\vx) = \vr_i(\vx)`, for :math:`i = 1,\dots, m`,
and :math:`\widehat{\vr}_{m+j}(\vx) =
\sqrt{\sigma}[\vx]_j` for :math:`j = 1,\dots, n`, where :math:`[\vx]_j`
denotes the :math:`j`\ th component of :math:`\vx`. We therefore have
that :math:`\nabla \widehat{\vr}_{m+j}(\vx) = \sqrt{\sigma}\ve^j` (where
:math:`[\ve^j]_i = \delta_{ij}`), and the second derivatives vanish.

*Solved implicitly...* **todo**

**If** :math:`\bf p \ne 2`, 
then we solve a least squares problem with
one additional degree of freedom. In this case the new function,
:math:`\widehat{\vr} : \mathbb{R}^{n}\rightarrow\mathbb{R}^{m+1}`, again
takes :math:`\widehat{\vr}_i(\vx) = \vr_i(\vx)`, for
:math:`i = 1,\dots, m`, but now
:math:`\widehat{\vr}_{m+1}(\vx) = \left(\frac{2\sigma}{p}\right)^{\frac{1}{2}}\|\vx\|^{\frac{p}{2}}.`
We therefore have that
:math:`\nabla \widehat{\vr}_{m+1}(\vx) = \left(\frac{2\sigma}{p}\right)^{\frac{1}{2}}\|\vx\|^{\frac{p-4}{2}}\vx^T`.
The second derivative is given by :math:`\nabla^2\widehat{\vr}_{m+1} = \left(\frac{2\sigma}{p}\right)^{\frac{1}{2}}\|\vx\|^{\frac{p-4}{2}}\left(I + \frac{\vx\vx^T}{\|\vx\|^2}\right).`

*Solved implicitly...* **todo**

.. _models:

The models
----------

A vital component of the algorithm is the choice of model employed.
There are four choices available, controlled by the parameter
``model`` of |nlls_options|.

``model = 1``
  this implements the **Gauss-Newton** model. Here we
  replace :math:`\vr(\iter[k]{\vx} + \vs)` by its first-order Taylor
  approximation, :math:`\vr(\iter{\vx}) + \iter{\vJ}\vs`. The model is
  therefore given by

  .. math::
     :label: gauss-newton-model

     m_k^{GN}(\vs) = \frac{1}{2} \|\vr(\iter{\vx}) + \iter{\vJ}\vs\|_\vW^2.

``model = 2`` 
  this implements the **Newton** model. Here, instead of
  approximating the residual, :math:`\vr(\cdot)`, we take as our model the
  second-order Taylor approximation of the function,
  :math:`F(\iter[k+1]{\vx}).` Namely, we use

  .. math::
     :label: newton-model
       m_k^{N}(\vs) = F(\iter{\vx}) + {\iter{\vg}}^T\vs + \frac{1}{2}\vs^T\left( {\iter{\vJ}}^T \vW \iter{\vJ} + \iter{\vH}\right) \vs,

  where :math:`\iter{\vg} = {\iter{\vJ}}^T\vW \vr(\iter{\vx})` and
  :math:`\iter{\vH} = \sum_{i=1}^m\iter[i]{r}(\iter{\vx}) \vW \nabla^2 \iter[i]{r}(\iter{\vx}).`
  Note that
  :math:`m_k^{N}(\vs) = m_k^{GN}(\vs) + \frac{1}{2}\vs^T\iter{\vH} \vs`.

  If the second derivatives of :math:`\vr(\cdot)` are not available (i.e.,
  the option ``exact_second_derivatives`` is set to ``false``, then the method approximates
  the matrix :math:`\iter{\vH}`; see :ref:`approx_hess`.

``model = 3``
  This implements a **hybrid** model. In practice the
  Gauss-Newton model tends to work well far away from the solution,
  whereas Newton performs better once we are near to the minimum
  (particularly if the residual is large at the solution). This option
  will try to switch between these two models, picking the model that is
  most appropriate for the step. In particular, we start using
  :math:`m_k^{GN}(\cdot)`, and switch to :math:`m_k^{N}(\cdot)` if
  :math:`\|{\iter{\vg}}\|_2 \leq \mathtt{hybrid\_tol} \frac{1}{2}\|\vr(\iter{\vx})\|^2_\vW`
  for more than ``hybrid_switch_its`` iterations in a row. If, in
  subsequent iterations, we fail to get a decrease in the function value,
  then the algorithm interprets this as being not sufficiently close to
  the solution, and thus switches back to using the Gauss-Newton model.

  **todo** :: make the switching algorithm cleraer
  
``model = 4`` 
   this implements a **Newton-tensor** model. This uses a
   second order Taylor approximation to the residual, namely

   .. math:: 

      r_i(\iter{\vx} + \vs) \approx (\iter{\vt}(\vs))_i := r_i(\iter{\vx}) + (\iter{\vJ})_i\vs + \frac{1}{2}\vs^T B_{ik}\vs,

   where :math:`(\iter{\vJ})_i` is the ith row of :math:`\iter{\vJ}`, and
   :math:`B_{ik}` is :math:`\nabla^2 r_i(\iter{\vx})`. We use this to
	 define our model

	 .. math::
	    :label: newton-tensor-model

	    m_k^{NT}(\vs) = \frac{1}{2}\|\vt_k(\vs)\|_\vW^2.


.. _approx_hess:

Approximating the Hessian
-------------------------

If the exact Hessian is not available, we 
approximate it using the method of Dennis, Gay, and Welsch.

**todo**



Subproblem solves
-----------------

The main algorithm  calls a number
of subroutines. The most vital is the subroutine ``calculate_step``, which
finds a step that minimizes the model chosen, subject to a globalization
strategy. The algorithm supports the use of two such strategies: using a
trust-region, and regularization. If Gauss-Newton, (quasi-)Newton, or a
hybrid method is used (``model = 1,2,3`` in |nlls_options|), 
then the model function is
quadratic, and the methods available to solve the subproblem are
described in :ref:`sec_tr` and :ref:`sec_reg`. 
If the Newton-Tensor model is selected (``model = 4`` in |nlls_options|), then this model
is not quadratic, and the methods available are described in
:ref:`newton-tensor-subproblem`.

Note that, when calculating the step, if the initial regularization
parameter :math:`\sigma` in :eq:`lsq` is non-zero,
then we must modify :math:`{\iter[k]{\tJ}}^T\iter[k]{\tJ}` to take into
account the Jacobian of the modified least squares problem being solved.
Practically, this amounts to making the change

.. math::

   {\iter[k]{\tJ}}^T\iter[k]{\tJ} = {\iter[k]{\tJ}}^T\iter[k]{\tJ} + 
    \begin{cases}
      \sigma I & \text{if }p = 2\\
      \frac{\sigma p}{2} \|\iter[k]{\vx}\|^{p-4}\iter[k]{\vx}{\iter[k]{\vx}}^T & \text{otherwise}
    \end{cases}.

.. _sec_tr:

The trust region method
^^^^^^^^^^^^^^^^^^^^^^^

If ``model = 1, 2,`` or ``3``, and ``type_of_method=1``, then we solve the subproblem :eq:`trsub`, and we take
as our next step the minimum of the model within some radius of the
current point. The method used to solve this is dependent on the control
parameter optionsnlls\_method. The algorithms called for each of the
options are listed below:

``nlls_method = 1``
    approximates the solution to :eq:`trsub`
    by using Powell’s dogleg method. This takes
    as the step a linear combination of the Gauss-Newton step and the
    steepest descent step, and the method used is described in Algorithm
    [alg:dogleg]. **todo**

``nlls_method = 2``
    solves the trust region subproblem using
    the trust region solver of Adachi, Iwata, Nakatsukasa, and Takeda. This
    reformulates the problem :eq:`trsub` as a generalized
    eigenvalue problem, and solves that. See [1]_ for more details.

``nlls_method = 3``
    this solves the trust region subproblem using a
    variant of the More-Sorensen method. In particular, we implement
    Algorithm 7.3.6 in Trust Region Methods by Conn, Gould and Toint [2]_.

``nlls_method = 4``
    this solves the trust region subproblem by
    first converting the problem into the form

    .. math:: \min_\vp \vw^T \vp + \frac{1}{2} \vp^T \vD \vp \quad {\rm s.t.} \quad \|\vp\| \leq \Delta,

    where :math:`\vD` is a diagonal matrix. We do this by performing an
    eigen-decomposition of the Hessian in the model. Then, we call the
    Galahad routine DTRS; see the Galahad [3] documentation for further
    details.

.. _sec_reg:

Regularization
^^^^^^^^^^^^^^

If ``model = 1, 2,`` or ``3``, and ``type_of_method=2``, then the next step is taken to be the
minimum of the model with a regularization term added
:eq:`regsub`. At present, only one method of solving
this subproblem is supported:

``nlls_method = 4``: 
  this solves the regularized subproblem by first
  converting the problem into the form

  .. math:: \min_\vp \vw^T \vp + \frac{1}{2} \vp^T \vD \vp + \frac{1}{p}\|\vp\|_2^p,

  where :math:`\vD` is a diagonal matrix. We do this by performing an
  eigen-decomposition of the Hessian in the model. Then, we call the
  Galahad routine DRQS; see the Galahad [3]_ documentation for further
  details.

.. _newton-tensor-subproblem:

Newton-Tensor subproblem
^^^^^^^^^^^^^^^^^^^^^^^^

If ``model=4``, then the non-quadratic Newton-Tensor model is used.
As such, none of the established subproblem solvers described in
:ref:`sec_tr` or :ref:`sec_reg` can be used.

If we use regularization (with :math:`p=2`), then the subproblem we need
to solve is of the form

.. math::
   :label: reg_newton_tensor_subproblem

   \min_\vs \frac{1}{2}\sum_{i=1}^mW_{ii}{(\vt_k(\vs))_i}^2 + \frac{1}{2\Delta_k}\|\vs\|_2^2

Note that :eq:`reg_newton_tensor_subproblem` is a
sum-of-squares, and as such can be solved by calling ral\_nlls
recursively. We support two options:

``inner_method = 1``
  if this option is selected, then ``nlls_solve()``
  is called to solve :eq:`newton-tensor-model` directly. The current
  regularization parameter of the ‘outer’ method is used as a base
  regularization in the ‘inner’ method, so that the (quadratic) subproblem
  being solved in the ‘inner’ call is of the form

  .. math:: \min_\vs \, m_k(\vs) + \frac{1}{2}\left(\frac{1}{\Delta_k} + \frac{1}{\delta_k}\right)\|\vs\|_B^2,

  where :math:`m_k(\vs)` is a quadratic model of
  :eq:`newton-tensor-model`, :math:`\Delta_k` is the (fixed)
  regularization parameter of the outer iteration, and :math:`\delta_k`
  the regularization parameter of the inner iteration, which is free to be
  updated as required by the method.

``inner_method = 2``
  in this case we use ``ral_nlls()`` to solve the
  regularized model :eq:`reg_newton_tensor_subproblem`)
  directly. The number of parameters for this subproblem is :math:`n+m`.
  Specifically, we have a problem of the form

  .. math::
     
     \min_\vs \frac{1}{2} \|\widehat{\vr}(\vs)\|_\vW^2,
     \quad \text{where }   
     (\widehat{\vr}(\vs))_i =
     \begin{cases}
     (\vt_k(\vs))_i &  1 \leq i \leq m \\
     \frac{1}{\sqrt{\Delta_k}}s_i& m+1 \leq i \leq n+m
     \end{cases}.

  This subproblem can then be solved using any of the methods described in
  :ref:`sec_tr` or :ref:`sec_reg`.

Accepting the step and updating the parameter
---------------------------------------------

Once a step has been suggested, we must decide whether or not to accept
the step, and whether the trust region radius or regularization
parameter, as appropriate, should grow, shrink, or remain the same.

These decisions are made with reference to a parameter, :math:`\rho`,
which measures the ratio of the actual reduction in the model to the
predicted reduction in the model. If this is larger than
``eta_successful`` in |nlls_options|, then the step 
ise accepted  **TODO** (see Line 28 of Algorithm[ alg:nlls_solve]).

The value of :math:`\Delta_k` then needs to be updated, if appropriate.
The package supports two options:

``tr_update_strategy = 1`` 
  a step-function is used to
  decide whether or not to increase or decrease :math:`\Delta_k`.

``tr_update_strategy = 2`` 
  a continuous function is used to make the decision.

The method used is outlined in Algorithm **todo** [alg:update_tr].



.. [1] Adachi, Satoru and Iwata, Satoru and Nakatsukasa, Yuji and Takeda, Akiko (2015). Solving the trust region subproblem by a generalized eigenvalue problem. Technical report, Mathematical Engineering, The University of Tokyo.
.. [2] Conn, A. R., Gould, N. I., & Toint, P. L. (2000). Trust region methods. SIAM.
.. [3] Gould, N. I., Orban, D., & Toint, P. L. (2003). GALAHAD, a library of thread-safe Fortran 90 packages for large-scale nonlinear optimization. ACM Transactions on Mathematical Software (TOMS), 29(4), 353-372.
.. [4] Nocedal, J., & Wright, S. (2006). Numerical optimization. Springer Science & Business Media.
