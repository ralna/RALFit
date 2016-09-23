.. include:: ../common/install.rst

Building the Python interface
-----------------------------

From the directory ``RALFit/libRALFit/``, run the ``distutils`` build script:

.. code::

   python setup.py build_ext --inplace

This will create the python module in a file ``ral_nlls.so``.

If you want to install it on a system level, instead use the following build
and installation commands

.. code::

   python setup.py build_ext
   python setup.py install
