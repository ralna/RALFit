<img alt="STFC logo" src="https://www.ukri.org/wp-content/uploads/2022/03/ukri-stfc-square-logo.png" width=100>

# RALFit

This readme provides instructions to compile the RALFit solver.

To read more about the solver, see
[the documentation](https://ralfit.readthedocs.io/projects/Fortran/en/latest/), the
[C interface](https://ralfit.readthedocs.io/projects/C/en/latest/), or the
[Python interface](https://ralfit.readthedocs.io/projects/Python/en/latest/).


## Requirements

* CMake (tested on 3.18 and 3.31) and GNU Make
* Modern Fortran compiler
* C compiler
* Python 3 (tested on 3.11.4) if building the Python package

## Configure

The package uses CMake to configure and Make build. The following
CMake options are briefly described

  * `-DCMAKE_BUILD_TYPE={Release|Debug}` defines the type of build, default `Release`
  * `-DBUILD_PYTHON={On|Off}` option to build the Python interface (package), default `On`
  * `-DSINGLE_PRECISION={Off|On}` selects which floating point type to
    use, default `Off` (double precision)

### Notes

 * Building the `Debug` version implies the _internal reference LAPACK
   and BLAS libraries_ are used, not the user-provided one.

 * Some of the example sources are only compiled when building the double
   precision floating point version of the library and are not available
   when building the single precision version.

## Build the library

### Download or clone
```{terminal}
git clone https://github.com/ralna/RALFit.git
```

### Python package
If building Python interface, then configure or load the virtual environment
```{terminal}
python3 -mvenv p3
. p3/bin/activate
```

### Configure
```{terminal}
cd RALFit/libRALFit`
cmake -S . -B build
cmake --build build --parallel 4
```

### Running the examples

The examples sources are located in the sub-folders of `libRALFit/examples/` and there are Fortran, C and Python examples.

The compiled example executables can be found in `build/examples/`


## Building with custom compiler and libraries

This section shows how to build RALFit using an accelerated compiler (AMD,
Intel, etc) and linking to performant LAPACK/BLAS libraries (e.g. AOCL, MKL).
We illustrate how to build using AMD's AOCC compiler and AOCL accelerated BLAS
and LAPACK libraries.

### Requirements

Prior to building an accelerated version of RALFit, install compiler and libraries.
Refer [here](https://www.amd.com/en/developer/aocc.html) to install and load [AMD AOCC compilers](https://www.amd.com/en/developer/aocc.html) and [here](https://www.amd.com/en/developer/aocl.html) to install and configure [AOCL BLAS and LAPACK](https://www.amd.com/en/developer/aocl.html).

Make sure the correct compiler version is available on the path:
```{terminal}
$ clang --version
AMD clang version 17.0.6 (CLANG: AOCC_5.1.0)
Target: x86_64-unknown-linux-gnu
Thread model: posix
InstalledDir: /home/u/amd/aocc-compiler-rel-5.1.0/bin
```
and a the same banner is printed for `flang --version` command.

### Configuring and building

To instruct CMake to compile to the native CPU architecture and use the AMD accelerated versions BLAS and LAPACK,
pass the flags `-DBLA_VENDOR=AOCL -DCMAKE_C_FLAGS_RELEASE="-O2 -march=native -DNDEBUG" -DCMAKE_Fortran_FLAGS_RELEASE="-O2 -march=native -DNDEBUG"`.

Make sure the correct compiler, BLAS, and LAPACK versions are found by CMake:
```{terminal}
$ CC=clang FC=flang cmake -S . -B build --fresh -DBLA_VENDOR=AOCL -DCMAKE_C_FLAGS_RELEASE="-O2 -march=native -DNDEBUG" -DCMAKE_Fortran_FLAGS_RELEASE="-O2 -march=native -DNDEBUG"
```
Should produce an output similar to
```{terminal}
...
-- Check for working Fortran compiler: /home/u/amd/aocc-compiler-rel-5.0.0/bin/flang - skipped
-- Check for working C compiler: /home/u/amd/aocc-compiler-rel-5.0.0/bin/clang - skipped
...
-- Found LAPACK: /home/u/amd/5.0.0/aocl/lib/libflame.so;/home/u/amd/5.0.0/aocl/lib/libblis.so;-fopenmp
-- The LAPACK library was found at/home/u/amd/5.0.0/aocl/lib/libflame.so/home/u/amd/5.0.0/gcc/lib/libblis.so-fopenmp
...
```

Now it is ready to build
```{terminal}
cmake --build build --parallel 4
```