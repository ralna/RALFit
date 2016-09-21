Consider fitting the function :math:`y(t) = x_1e^{x_2 t}` to data :math:`(\bm{t}, \bm{y})`
using a non-linear least squares fit.\


The residual function is given by

.. math:: 

   r_i(\vx)  = x_1 e^{x_2 t_i} - y_i.


We can calculate the Jacobian and Hessian of the residual as

.. math::

   \nabla r_i(\vx) = \left(\begin{array}{cc}
      e^{x_2 t_i} &
      t_i x_1 e^{x_2 t_i}
      \end{array}\right),

.. math:: 
   \nabla^2 r_i(\vx) = \left(\begin{array}{cc}
      0                 & t_i e^{x_2 t_i}    \\
      t_i e^{x_2 t_i}     & x_1 t_i^2 e^{x_2 t_i}
   \end{array}\right).

For some data points, :math:`y_i`, :math:`t_i`, :math:`(i = 1,\dots,m)` the user must return

.. math::
   
   \vr(\vx) = \begin{bmatrix}
      r_1(\vx)\\
      \vdots \\
      r_m(\vx)
    \end{bmatrix}, \quad   \vJ(\vx) = 
    \begin{bmatrix}
      \nabla r_1(\vx) \\
      \vdots \\
      \nabla r_m(\vx) \\
    \end{bmatrix}, \quad 
    Hf(\vx) = 
    \sum_{i=1}^m
    (\vr)_i \nabla^2 r_i(\vx),

where, in the case of the Hessian, :math:`(\vr)_i` is the :math:`i\text{th}` component of a residual vector passed to the user.

.. list-table::
   
    * - :math:`i`
      - 1 
      - 2 
      - 3  
      - 4  
      - 5
    * - :math:`t_i`
      - 1
      - 2
      - 4 
      - 5 
      - 8
    * - :math:`y_i`
      - 3
      - 4
      - 6
      - 11
      - 20

and initial guess :math:`\vx = (2.5, 0.25)`, the following code performs the fit (with no weightings, i.e., :math:`\vW = \vI`).
