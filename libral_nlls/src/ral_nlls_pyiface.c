#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

#include "ral_nlls.h"

struct callback_data {
   PyObject* f;
   PyObject* J;
   PyObject* Hr;
   PyObject* params;
};

static
int eval_f(int n, int m, const void *params, const double *x, double *f) {
   // Recover our datatype
   const struct callback_data *data = (struct callback_data*) params;

   // Copy x into Python array
   npy_intp xdim[] = {n};
   PyArrayObject* xpy = (PyArrayObject*) PyArray_SimpleNew(1, xdim, NPY_DOUBLE);
   double* xval = (double*) PyArray_DATA(xpy);
   for(int i=0; i<n; ++i)
      xval[i] = x[i];

   // Call routine
   PyObject *arglist;
   if(data->params)  arglist = Py_BuildValue("(OO)", xpy, data->params);
   else              arglist = Py_BuildValue("(O)", xpy);
   PyObject *result = PyObject_CallObject(data->f, arglist);
   Py_DECREF(arglist);
   Py_DECREF(xpy);
   if(!result) return -1;

   // Extract result
   PyArrayObject* farray = (PyArrayObject*) PyArray_FROM_OTF(result, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
   if(farray == NULL) {
      PyErr_SetString(PyExc_RuntimeError, "Error extracting array from r call");
      Py_DECREF(result);
      return -1;
   }
   if(PyArray_NDIM(farray) != 1) {
      PyErr_SetString(PyExc_RuntimeError, "f() must return rank-1 array");
      Py_DECREF(farray);
      Py_DECREF(result);
      return -2;
   }
   const double *fval = (double*) PyArray_DATA(farray);
   for(int i=0; i<m; ++i) {
      f[i] = fval[i];
   }
   Py_DECREF(farray);
   Py_DECREF(result);

   return 0; // Success
}

static
int eval_J(int n, int m, const void *params, const double *x, double *J) {
   // Recover our datatype
   const struct callback_data *data = (struct callback_data*) params;

   // Copy x into Python array
   npy_intp xdim[] = {n};
   PyArrayObject* xpy = (PyArrayObject*) PyArray_SimpleNew(1, xdim, NPY_DOUBLE);
   double* xval = (double*) PyArray_DATA(xpy);
   for(int i=0; i<n; ++i)
      xval[i] = x[i];

   // Call routine
   PyObject *arglist;
   if(data->params)  arglist = Py_BuildValue("(OO)", xpy, data->params);
   else              arglist = Py_BuildValue("(O)", xpy);
   PyObject *result = PyObject_CallObject(data->J, arglist);
   Py_DECREF(arglist);
   Py_DECREF(xpy);
   if(!result) return -1;

   // Extract result
   PyArrayObject* Jarray = (PyArrayObject*) PyArray_FROM_OTF(result, NPY_FLOAT64, NPY_ARRAY_IN_FARRAY);
   if(Jarray == NULL) {
      PyErr_SetString(PyExc_RuntimeError, "Error extracting array from J call");
      Py_DECREF(result);
      return -1;
   }
   if(PyArray_NDIM(Jarray) != 2) {
      PyErr_SetString(PyExc_RuntimeError, "J() must return rank-2 array");
      Py_DECREF(Jarray);
      Py_DECREF(result);
      return -2;
   }
   const double *Jval = (double*) PyArray_DATA(Jarray);
   for(int i=0; i<m*n; ++i) {
      J[i] = Jval[i];
   }
   Py_DECREF(Jarray);
   Py_DECREF(result);

   return 0; // Success
}

static
int eval_Hr(int n, int m, const void *params, const double *x, const double *r, double *Hr) {
   // Recover our datatype
   const struct callback_data *data = (struct callback_data*) params;

   // Copy x into Python array
   npy_intp xdim[] = {n};
   PyArrayObject* xpy = (PyArrayObject*) PyArray_SimpleNew(1, xdim, NPY_DOUBLE);
   double* xval = (double*) PyArray_DATA(xpy);
   for(int i=0; i<n; ++i)
      xval[i] = x[i];

   // Copy r into Python array
   npy_intp rdim[] = {m};
   PyArrayObject* rpy = (PyArrayObject*) PyArray_SimpleNew(1, rdim, NPY_DOUBLE);
   double* rval = (double*) PyArray_DATA(rpy);
   for(int i=0; i<m; ++i)
      rval[i] = r[i];

   // Call routine
   PyObject *arglist;
   if(data->params)  arglist = Py_BuildValue("(OOO)", xpy, rpy, data->params);
   else              arglist = Py_BuildValue("(OO)", xpy, rpy);
   PyObject *result = PyObject_CallObject(data->Hr, arglist);
   Py_DECREF(arglist);
   Py_DECREF(xpy);
   Py_DECREF(rpy);
   if(!result) return -1;

   // Extract result
   PyArrayObject* Hrarray = (PyArrayObject*) PyArray_FROM_OTF(result, NPY_FLOAT64, NPY_ARRAY_IN_FARRAY);
   if(Hrarray == NULL) {
      PyErr_SetString(PyExc_RuntimeError, "Error extracting array from Hr call");
      Py_DECREF(result);
      return -1;
   }
   if(PyArray_NDIM(Hrarray) != 2) {
      PyErr_SetString(PyExc_RuntimeError, "Hr() must return rank-2 array");
      Py_DECREF(Hrarray);
      Py_DECREF(result);
      return -2;
   }
   const double *Hrval = (double*) PyArray_DATA(Hrarray);
   for(int i=0; i<n*n; ++i) {
      Hr[i] = Hrval[i];
   }
   Py_DECREF(Hrarray);
   Py_DECREF(result);

   return 0; // Success
}

static
bool set_opts(struct ral_nlls_options *options, PyObject *pyoptions) {
   ral_nlls_default_options_d(options);

   if(!pyoptions) return true; // Just use defaults

   PyObject *key, *value;
   Py_ssize_t pos = 0;
   while(PyDict_Next(pyoptions, &pos, &key, &value)) {
      const char* key_name = PyString_AsString(key);
      if(!key_name) {
         PyErr_SetString(PyExc_RuntimeError, "Non-string option, can't interpret!");
         return false;
      }
      if(strcmp(key_name, "print_level")==0) {
         long v = PyInt_AsLong(value);
         if(v==-1 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['print_level'] must be an integer.");
            return false;
         }
         options->print_level = (int) v;
         continue;
      }
      if(strcmp(key_name, "maxit")==0) {
         long v = PyInt_AsLong(value);
         if(v==-1 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['maxit'] must be an integer.");
            return false;
         }
         options->maxit = (int) v;
         continue;
      }
      if(strcmp(key_name, "model")==0) {
         long v = PyInt_AsLong(value);
         if(v==-1 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['model'] must be an integer.");
            return false;
         }
         options->model = (int) v;
         continue;
      }
      if(strcmp(key_name, "nlls_method")==0) {
         long v = PyInt_AsLong(value);
         if(v==-1 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['nlls_method'] must be an integer.");
            return false;
         }
         options->nlls_method = (int) v;
         continue;
      }
      if(strcmp(key_name, "stop_g_absolute")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['stop_g_absolute'] must be a float.");
            return false;
         }
         options->stop_g_absolute = v;
         continue;
      }
      if(strcmp(key_name, "stop_g_relative")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['stop_g_relative'] must be a float.");
            return false;
         }
         options->stop_g_relative = v;
         continue;
      }
      if(strcmp(key_name, "relative_tr_radius")==0) {
         long v = PyInt_AsLong(value);
         if(v==-1 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['relative_tr_radius'] must be an integer.");
            return false;
         }
         options->relative_tr_radius = (int) v;
         continue;
      }
      // If we reach this point, unreconised option
      char errmsg[200];
      snprintf(errmsg, 200, "Bad key options['%s']\n", key_name);
      PyErr_SetString(PyExc_RuntimeError, errmsg);
      return false;
   }

   return true; // success
}

static PyObject*
make_info_dict(const struct ral_nlls_inform *inform) {
   PyObject *pyinfo = PyDict_New();

   PyDict_SetItemString(pyinfo, "iter", PyInt_FromLong(inform->iter));
   PyDict_SetItemString(pyinfo, "f_eval", PyInt_FromLong(inform->f_eval));
   PyDict_SetItemString(pyinfo, "g_eval", PyInt_FromLong(inform->g_eval));
   PyDict_SetItemString(pyinfo, "h_eval", PyInt_FromLong(inform->h_eval));
   PyDict_SetItemString(pyinfo, "convergence_normf",
         PyInt_FromLong(inform->convergence_normf)
         );
   PyDict_SetItemString(pyinfo, "resinf", PyFloat_FromDouble(inform->resinf));
   PyDict_SetItemString(pyinfo, "gradinf", PyFloat_FromDouble(inform->gradinf));
   PyDict_SetItemString(pyinfo, "obj", PyFloat_FromDouble(inform->obj));
   PyDict_SetItemString(pyinfo, "norm_g", PyFloat_FromDouble(inform->norm_g));
   PyDict_SetItemString(pyinfo, "scaled_g",
         PyFloat_FromDouble(inform->scaled_g)
         );

   return pyinfo;
}

/*
 * x = ral_nlls.solve(x0, f, J=None, Hr=None, params=(), options={})
 */
static PyObject*
ral_nlls_solve(PyObject* self, PyObject* args, PyObject* keywds)
{
   PyObject *x0ptr=NULL, *options_ptr=NULL;
   PyObject *arglist=NULL, *result=NULL;
   PyArrayObject *x0=NULL, *f=NULL, *x=NULL;

   struct callback_data data;
   data.f = NULL; data.J = NULL; data.Hr = NULL; data.params = NULL;

   static char *kwlist[] = {"x0", "r", "J", "Hr", "params", "options", NULL};
   if(!PyArg_ParseTupleAndKeywords(args, keywds, "OOO|OOO!", kwlist,
            &x0ptr,
            &data.f,
            &data.J,
            &data.Hr,
            &data.params,
            &PyDict_Type, &options_ptr)
         )
      return NULL;

   /* x0 */
   x0 = (PyArrayObject*) PyArray_FROM_OTF(x0ptr, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
   if(x0 == NULL) return NULL;
   if(PyArray_NDIM(x0) != 1) {
      PyErr_SetString(PyExc_RuntimeError, "x0 must be a rank-1 array");
      goto fail;
   }
   npy_intp* xdim = PyArray_DIMS(x0);
   int n = xdim[0];

   /* Determine m by making call to f */
   if(data.params)   arglist = Py_BuildValue("(OO)", x0, data.params);
   else              arglist = Py_BuildValue("(O)", x0);
   result = PyObject_CallObject(data.f, arglist);
   if(!result) goto fail;
   Py_DECREF(arglist); arglist=NULL;
   f = (PyArrayObject*) PyArray_FROM_OTF(result, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
   if(f == NULL) {
      PyErr_SetString(PyExc_RuntimeError, "Error extracting array from f call");
      goto fail;
   }
   if(PyArray_NDIM(f) != 1) {
      PyErr_SetString(PyExc_RuntimeError, "r() must return rank-1 array");
      goto fail;
   }
   npy_intp* fdim = PyArray_DIMS(f);
   int m = fdim[0];
   Py_DECREF(f); f=NULL;
   Py_DECREF(result); result=NULL;

   /* Construct return array x[] and set equal to x0 */
   x = (PyArrayObject*) PyArray_SimpleNew(1, xdim, NPY_DOUBLE);
   double* xval = (double*) PyArray_DATA(x);
   const double* x0val = (double*) PyArray_DATA(x0);
   for(int i=0; i<n; ++i)
      xval[i] = x0val[i];

   /* Call RAL_NLLS */
   struct ral_nlls_options options;
   if(!set_opts(&options, options_ptr)) goto fail;
   struct ral_nlls_inform inform;
   if(data.Hr) {
      nlls_solve_d(n, m, xval, eval_f, eval_J, eval_Hr, &data, &options, &inform);
   } else {
      options.exact_second_derivatives = false;
      nlls_solve_d(n, m, xval, eval_f, eval_J, NULL, &data, &options, &inform);
   }
   switch(inform.status) {
      case 0: // Clean exit
         break;
      case -1: // Exceeded max itr
         PyErr_SetString(PyExc_RuntimeError,
               "Exceeded maximum number of iterations");
         goto fail;
      case -2: // Error return from evaluation of f/J/Hr
         // No error msg, allow existing one to propagate
         goto fail;
      case -3: // Unsupported choice of model
         PyErr_SetString(PyExc_RuntimeError,
               "Bad model");
         goto fail;
      default: ; // empty statement for language conformatity.
         char errmsg[100];
         sprintf(errmsg, "NLLS_SOLVE failed with unrecognised error code %d\n",
               inform.status);
         PyErr_SetString(PyExc_RuntimeError, errmsg);
         goto fail;
   }

   /* Free references and return solution */
   Py_DECREF(x0); x0=NULL;
   PyObject *pyinfo = make_info_dict(&inform);
   return Py_BuildValue("(OO)", x, pyinfo);

   fail:
   Py_XDECREF(arglist); Py_XDECREF(result);
   Py_XDECREF(x0); Py_XDECREF(f); Py_XDECREF(x);
   return NULL;
}

static PyMethodDef RalNllsMethods[] = {
   {"solve", (PyCFunction)ral_nlls_solve, METH_VARARGS | METH_KEYWORDS,
    "Solve a non-linear least squares problem.\n"
    "   x = ral_nlls.solve(x0, r, J, Hr=None, params=None, options={})"
   },
   {NULL, NULL, 0, NULL} /* Sentinel */
};

PyMODINIT_FUNC
initral_nlls(void)
{
   (void) Py_InitModule("ral_nlls", RalNllsMethods);
   import_array();
}
