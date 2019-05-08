.. default-domain:: Python

.. include:: macros.rst

How to use the package
======================

Overview
--------

Calling sequences
^^^^^^^^^^^^^^^^^

.. py:module:: ral_nlls

This module is a binary module that provides
non-linear least squares solvers.
It requires the :py:mod:`numpy` module.

The module is imported by issuing the command

.. code-block:: Python

   import numpy
   import ral_nlls
 

The user can then call the procedure:

|nlls_solve| solves the non-linear least squares problem.

At present the Python interface only supports the option :math:`W=I`.


Argument list and calling sequence
----------------------------------



To solve the non-linear least squares problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../common/subroutines.rst

.. py:function:: solve(x0, r, J, Hr=None, params=None, options=None) -> (x, info)

   Solves the non-linear least squares problem.
   
   :param x0:  holds the initial guess for the solution
   :type x0: array_like with shape (n,)

   :param r: given a point :math:`{\bm x} _{k}^{}`, returns the vector :math:`{\bm r} ({\bm x} _{k}^{})`. The function must have the signature ``r(x[,params])``, where ``x`` is an ndarray of shape (n,) and `params` is as defined below.  It must return a 1-d array_like of shape (m,).
   :type r: callable
			
   :param J: given a point :math:`{\bm x} _{k}^{}`, returns the :math:`m \times n` Jacobian matrix, :math:`{\bm J} _{k}^{}`, of :math:`{\bm r}` evaluated at :math:`{\bm x} _{k}^{}`. The function must have the signature ``J(x[,params])``, where ``x`` is an ndarray of shape (n,) and `params` is as defined below.  It must return an array_like of shape (m,n).
   :type J: callable

   :param Hf: given vectors :math:`{\bm x} \in \mathbb{R}^n` and :math:`{\bm r} \in \mathbb{R}^m`, returns the quantity :math:`\sum_{i=1}^m ( {\bm r} )_i \nabla^2  {\bm r} _i ( {\bm x} )`. The function must have the signature ``Hf(x[,params])``, where ``x`` is an ndarray of shape (n,).  It must return an array_like of shape (n,n) and `params` is as defined below.  If ``'exact_second_derivative' : F`` in |nlls_options|, then this is not referenced.
   :type Hf: None or callable

   :param params: holds parameters to be passed to the routines ``r``, ``J`` and ``Hf``.  
   :type params: None or tuple

   :param options: |options|
   :type options: None or dict


   :returns: The first component contains the approximate solution. The second component provides information about the execution of the subroutine.
   :rtype: tuple( nparray(n,), dict)

Data types
----------

.. _options:

The dictionary for holding options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../common/options.rst


Controlling data is sent to the subroutine using a Python dictionary.  A description of the possible options, together with their default values, follows below:
   
   **Printing Controls**
   
   * out (int)

		 |out|
		 Default is 6.

   * print_level (int)

		 |print_level|
		 
		 .. include:: ../common/options_print_level.txt
			      
		 Default is 0.

   * print_options (bool)

     |print_options|

     Default is false.

   * print_header (int)

     |print_header|

     Default is 30.
      
   **Choice of Algorithm**

   *  model (int)

		 |model| 
		 
		 .. include:: ../common/options_model.txt
			      
		 Default is 3.

   *   type_of_method (int)
		 
		 |type_of_method|
		 
		 .. include:: ../common/options_type_of_method.txt
			      
		 Default is 1.

   *  nlls_method (int)

		 |nlls_method|
				       
		 .. include:: ../common/options_nlls_method.txt
		    
		 Default is 4.
						    
   *  exact_second_derivatives (bool)

		 |exact_second_derivatives|
		 Default is false.


   **Stopping rules**

   *  maxit (int)
		 
		 |maxit|
		 Default is 100

   *  stop_g_absolute (float)
		 
		 |stop_g_absolute|
		 Defualt is 1e-5.
					   
   *  stop_g_relative (float)

		 |stop_g_relative|
		 Default is 1e-8.
					   
   *  stop_f_absolute (float)
		 
		 |stop_f_absolute|
		 Default is 1e-5.
					   
   *  stop_f_relative  (float)
		 
		 |stop_f_relative|
		 Default is 1e-8.

   *  stop_s (float)
		 
		 |stop_s|
		 Default is ``eps``.
   
   **Trust region radius/regularization behaviour**
      
   *  relative_tr_radius (int)
		 
		 |relative_tr_radius|
		 Default is 0.
					      
   *  initial_radius_scale (float)

		 |initial_radius_scale|
		 Default is 1.0.

   *  initial_radius (float)

		 |initial_radius|
		 Default is 100.0.

   *  maximum_radius (float)

		 |maximum_radius|
		 Default is 1e8.
		 
   *  eta_successful  (float)
      
		 |eta_successful|
		 Default is 1e-8.
		 .. success_but_reduce is also available, but not documented

   *  eta_very_successful (float)
		 
		 |eta_very_successful|
		 Default is 0.9.

   *  eta_too_successful (float)

		 |eta_too_successful|
		 Default is 2.0.

   *  radius_increase  (float)
		 
		 |radius_increase|
		 Default is 2.0.

   *  radius_reduce (float)

		 |radius_reduce|
		 Default is 0.5.

   *  tr_update_strategy (int)

		 |tr_update_strategy|

		 .. include:: ../common/options_tr_update_strategy.txt
		    
		 Default is 1.

    *  reg_order (float)
		 
		 |reg_order|
		 Default is 0.0.

	 
   **Scaling options**
   
   *  scale (int)

		 |scale|
		 Default is 1.
			

   *  scale_trim_max (bool)

		 |scale_trim_max|
		 Default is true.

   *  scale_max (float)
      
		 |scale_max|
		 Default is 1e11.

   *   scale_trim_min (bool)

		 |scale_trim_min|
		 Default is true.

   *  scale_min (float)
		 
		 |scale_min|
		 Default is 1e-11.

   *  scale_require_increase (bool)

		 |scale_require_increase|
		 Default is false.
      

   **Hybrid method options**  These options are used if ``model=3``

   *  hybrid_switch (float)

		 |hybrid_switch|
		 Default is 0.1.

   *  hybrid_tol (float)

		 |hybrid_tol|
		 Default is 2.0.	 
 
   *   hybrid_switch_its (int)
		 
		 |hybrid_switch_its|
		 Default is 1.

   **Newton-Tensor options** These options are used if ``model=4``
					     		 
   *  inner_method (int)
		 
		 |inner_method| 

		 .. include:: ../common/options_inner_method.txt		      
		 Default is 2.

   **More-Sorensen options**  These options are used if ``nlls_method=3``
   
   *  more_sorensen_maxits (int)
		 
		 |more_sorensen_maxits|
		 Default is 500.

   *   more_sorensen_maxits (int)
      
		 |more_sorensen_maxits|
		 Default is 3.

   *  more_sorensen_shift (int)
		 		 
		 |more_sorensen_shift|
		 Default is 1e-13.

   *  more_sorensen_tiny (float)
		 
		 |more_sorensen_tiny|
		 Default is 10.0 * ``eps``.

   *  more_sorensen_tol (float)

		 |more_sorensen_tol|
		 Default is 1e-3.
						  
   **Other options**
					     
   *  calculate_svd_J (bool)
		 
		 |calculate_svd_J|
		 Default is false.

   *  output_progress_vectors (bool)
		 
		 |output_progress_vectors|
		 Default is false.

   **Internal options to help solving a regularized problem implicitly**

   *  regularization (int)
      
		 |regularization|
					  
		 .. include:: ../common/options_regularization.txt
		 Default is 0.
					  
   *  regularization_term (float)

		 |regularization_term|
		 Default is 0.0.
						
   *  regularization_power (float)
      
		 |regularization_power|
		 Default is 0.0.


.. _info:

The dictionary for holding information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../common/inform.rst

This is used to hold information about the progress of the algorithm.

   *  iter (int)

		 |iter|

   *  f_eval  (int)
		 
		 |f_eval|
		      
   *  g_eval (int)
		 
		 |g_eval|

   *  h_eval (int)
      
		 |h_eval|

   *  convergence_normf (int)
		 
		 |convergence_normf|

   *  convergence_normf (int)
		 
		 |convergence_normg|
				
   *  convergence_normf (int)

		 |convergence_norms|
				   
   *  obj (float)
		 
		 |obj|

   *  norm_g (float)
			       
		 |norm_g|

   *  scaled_g (float)

		 |scaled_g|

   *  step (float)
     |step|

