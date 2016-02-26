from distutils.core import setup, Extension

ral_nlls_module = Extension(
        'ral_nlls',
        sources = ['src/ral_nlls_pyiface.c'],
        libraries = ['ral_nlls', 'lapack', 'blas'],
        library_dirs = ['build'],
        extra_compile_args = ['-std=c99'],
        include_dirs = ['include'],
        )

setup (name='RalNLLS', version='1.0', description='Nonlinear least squares',
        ext_modules=[ral_nlls_module])
