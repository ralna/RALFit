<table border="0"><tr style="background-color: white;">
<td><img alt="STFC logo" src=".github/assets/ukri-stfc-square-logo.png" width=100></td>
<td vertical-align='center'><img alt="AMD logo" src=".github/assets/amd-official-logo.jpg" width=200></td></tr></table>

# RALFit

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

[![Fortran doc](https://readthedocs.org/projects/ralfit-fortran/badge/?version=latest)](http://ralfit.readthedocs.io/projects/Fortran/en/latest/?badge=latest)
[http://ralfit.readthedocs.io/projects/Fortran/en/latest/](Fortran documentation)

[![C doc](https://readthedocs.org/projects/ralfit-c/badge/?version=latest)](http://ralfit.readthedocs.io/projects/C/en/latest/?badge=latest)
[http://ralfit.readthedocs.io/projects/C/en/latest/](C interface documentation)

[![Python doc](https://readthedocs.org/projects/ralfit-python/badge/?version=latest)](http://ralfit.readthedocs.io/projects/Python/en/latest/?badge=latest)
[http://ralfit.readthedocs.io/projects/Python/en/latest/] (Python interface documentation)
