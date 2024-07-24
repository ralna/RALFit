# Building the documentation

The documentation for RALFit is built using Sphinx.

To install the requirements for building the documentation, run
```
python -m pip install libRALFit/doc/requirements.txt
```

RALFit's docs are designed so that each supported language has a separate set of documentation, with common elements shared. To build the docs for a given Language:
```
cd libRALFit/doc/<Language>/
make html
```
where `<Language>` is Fortran, C, or Python.