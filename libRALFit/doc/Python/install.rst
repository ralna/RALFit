.. include:: ../common/install.rst

Building the Python interface
-----------------------------

Run the `distutils` build script:

.. code::

   python setup.py build_ext --inplace

This will create a ``ral_nlls.so`` that is the python module.

If you want to install it on a system level, instead use the following build
and installation commands

.. code::

   python setup.py build_ext
   python setup.py install
