![STFC logo](http://www.stfc.ac.uk/stfc/includes/themes/MuraSTFC/assets/legacy/2473_web_2.png)

# RALFit

[![travis status](https://travis-ci.org/ralna/RALFit.svg?branch=master)](https://travis-ci.org/ralna/RALFit)
[![codecov](https://codecov.io/gh/ralna/RALFit/branch/master/graph/badge.svg)](https://codecov.io/gh/ralna/RALFit)

A non-linear least squares solver that is primarily developed by the Numerical Analysis group at STFC Rutherford Appleton Laboratory hsl@stfc.ac.uk.

## Installation

### Requirements

RALFit has been tested on Linux, although it should work on other platforms.  It requires Fortran and C compilers, and BLAS and LAPACK must be installed on the system

### Compilation instructions

To compile, move to the `libRALFit` directory and issue the commands:
```
mkdir build
cd build
cmake ..
make
```

The package is written in modern Fortran, and we provide a number of interfaces for other languages.  For documentation, follow the links below:

[![Fortran doc](https://readthedocs.org/projects/ralfit-fortran/badge/?version=latest)](http://ralfit.readthedocs.io/projects/Fortran/en/latest/?badge=latest) Documentation for Fortran version.

[![C doc](https://readthedocs.org/projects/ralfit-c/badge/?version=latest)](http://ralfit.readthedocs.io/projects/C/en/latest/?badge=latest) Documentation for C interface.

[![Python doc](https://readthedocs.org/projects/ralfit-python/badge/?version=latest)](http://ralfit.readthedocs.io/projects/Python/en/latest/?badge=latest) Documentation for Python interface.
