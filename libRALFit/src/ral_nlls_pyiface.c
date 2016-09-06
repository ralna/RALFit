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
PyObject* build_arglist(Py_ssize_t sz, PyObject* extra) {
   Py_ssize_t als = sz;
   if(extra) als += PyTuple_Size(extra);
   PyObject *arglist = PyTuple_New(als);
   if(extra)
      for(int i=0; i<PyTuple_Size(extra); ++i) {
         PyObject *obj = PyTuple_GET_ITEM(extra, i);
         PyTuple_SET_ITEM(arglist, i+sz, obj); Py_INCREF(obj);
      }
   return arglist;
}

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
   PyObject *arglist = build_arglist(1, data->params);
   PyTuple_SET_ITEM(arglist, 0, (PyObject*) xpy); Py_INCREF(xpy);
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
   PyObject *arglist = build_arglist(1, data->params);
   PyTuple_SET_ITEM(arglist, 0, (PyObject*) xpy); Py_INCREF(xpy);
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
   PyObject *arglist = build_arglist(2, data->params);
   PyTuple_SET_ITEM(arglist, 0, (PyObject*) xpy); Py_INCREF(xpy);
   PyTuple_SET_ITEM(arglist, 1, (PyObject*) rpy); Py_INCREF(rpy);
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
      if(strcmp(key_name, "initial_radius_scale")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['initial_radius_scale'] must be a float.");
            return false;
         }
         options->initial_radius_scale = v;
         continue;
      }
      if(strcmp(key_name, "initial_radius")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['initial_radius'] must be a float.");
            return false;
         }
         options->initial_radius = v;
         continue;
      }
      if(strcmp(key_name, "maximum_radius")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['maximum_radius'] must be a float.");
            return false;
         }
         options->maximum_radius = v;
         continue;
      }
      if(strcmp(key_name, "eta_successful")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['eta_successful'] must be a float.");
            return false;
         }
         options->eta_successful = v;
         continue;
      }
      if(strcmp(key_name, "eta_success_but_reduce")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['eta_success_but_reduce'] must be a float.");
            return false;
         }
         options->eta_success_but_reduce = v;
         continue;
      }
      if(strcmp(key_name, "eta_very_successful")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['eta_very_successful'] must be a float.");
            return false;
         }
         options->eta_very_successful = v;
         continue;
      }
      if(strcmp(key_name, "eta_too_successful")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['eta_too_successful'] must be a float.");
            return false;
         }
         options->eta_too_successful = v;
         continue;
      }
      if(strcmp(key_name, "radius_increase")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['radius_increase'] must be a float.");
            return false;
         }
         options->radius_increase = v;
         continue;
      }
      if(strcmp(key_name, "radius_reduce")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['radius_reduce'] must be a float.");
            return false;
         }
         options->radius_reduce = v;
         continue;
      }
      if(strcmp(key_name, "radius_reduce_max")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['radius_reduce_max'] must be a float.");
            return false;
         }
         options->radius_reduce_max = v;
         continue;
      }
      if(strcmp(key_name, "tr_update_strategy")==0) {
	long v = PyInt_AsLong(value);
	if(v==-1 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['tr_update_strategy'] must be a float.");
            return false;
	}
	options->tr_update_strategy = v;
	continue;
      }
      if(strcmp(key_name, "hybrid_switch")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['hybrid_switch'] must be a float.");
            return false;
         }
         options->hybrid_switch = v;
         continue;
      }
      // bool: exact_second_derivatives
      // bool: subproblem_eig_fact
      //      if (strcmp(key_name, "subproblem_eig_fact")==0) {
      //	bool v = PyBool
      //      }
      if(strcmp(key_name, "scale")==0) {
	long v = PyInt_AsLong(value);
	if(v==-1 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['scale'] must be a float.");
            return false;
	}
	options->scale = v;
	continue;
      }
      if(strcmp(key_name, "scale_max")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['scale_max'] must be a float.");
            return false;
         }
         options->scale_max = v;
         continue;
      }
      if(strcmp(key_name, "scale_min")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['scale_min'] must be a float.");
            return false;
         }
         options->scale_min = v;
         continue;
      }
      // bool :: scale_trim_min
      // bool :: scale_trim_max
      // bool :: scale_require_increase
      // bool :: calculate_svd
      if(strcmp(key_name, "more_sorensen_maxits")==0) {
	long v = PyInt_AsLong(value);
	if(v==-1 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['more_sorensen_maxits'] must be a float.");
            return false;
	}
	options->more_sorensen_maxits = v;
	continue;
      }
      if(strcmp(key_name, "more_sorensen_shift")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['more_sorensen_shift'] must be a float.");
            return false;
         }
         options->more_sorensen_shift = v;
         continue;
      }
      if(strcmp(key_name, "more_sorensen_tiny")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['more_sorensen_tiny'] must be a float.");
            return false;
         }
         options->more_sorensen_tiny = v;
         continue;
      }
      if(strcmp(key_name, "more_sorensen_tol")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['more_sorensen_tol'] must be a float.");
            return false;
         }
         options->more_sorensen_tol = v;
         continue;
      }
      if(strcmp(key_name, "more_sorensen_tol")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['more_sorensen_tol'] must be a float.");
            return false;
         }
         options->more_sorensen_tol = v;
         continue;
      }
      if(strcmp(key_name, "hybrid_tol")==0) {
         double v = PyFloat_AsDouble(value);
         if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['hybrid_tol'] must be a float.");
            return false;
         }
         options->hybrid_tol = v;
         continue;
      }
      if(strcmp(key_name, "hybrid_switch_its")==0) {
	long v = PyInt_AsLong(value);
	if(v==-1 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['hybrid_switch_its'] must be a float.");
            return false;
	}
	options->hybrid_switch_its = v;
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
   //   PyDict_SetItemString(pyinfo, "resinf", PyFloat_FromDouble(inform->resinf));
   //   PyDict_SetItemString(pyinfo, "gradinf", PyFloat_FromDouble(inform->gradinf));
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
   arglist = build_arglist(1, data.params);
   PyTuple_SET_ITEM(arglist, 0, (PyObject*) x0); Py_INCREF(x0);
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
     nlls_solve_d(n, m, xval, eval_f, eval_J, eval_Hr, &data, &options, &inform,NULL);
   } else {
      options.exact_second_derivatives = false;
      nlls_solve_d(n, m, xval, eval_f, eval_J, NULL, &data, &options, &inform,NULL);
   }
   switch(inform.status) {
   case 0: // Clean exit
     break;
   case -1: // Exceeded max itr
     PyErr_SetString(PyExc_RuntimeError,
		     "Exceeded maximum number of iterations");
     goto output;
   case -2: // Error return from evaluation of f/J/Hr
     // No error msg, allow existing one to propagate
     goto fail;
   case -3: // Unsupported choice of model
     PyErr_SetString(PyExc_RuntimeError,
		     "Bad model");
     goto fail;
   case -4: // Error from external
     PyErr_SetString(PyExc_RuntimeError,
		     "External routine gave an error");
     goto fail;
   case -5: // Unsupported method
     PyErr_SetString(PyExc_RuntimeError,
		     "Bad method");
     goto fail;
   case -6: // Allocation error
     PyErr_SetString(PyExc_RuntimeError,
		     "Error allocating memory");
     goto fail;
   case -7: // max_tr_reductions
     PyErr_SetString(PyExc_RuntimeError,
		     "Max number of TR reductions");
     goto output;
   case -8: // no progress in x
     PyErr_SetString(PyExc_RuntimeError,
		     "No progress in x");
     goto output;
   case -9: //n_gt_m
     PyErr_SetString(PyExc_RuntimeError,
		     "n > m");
     goto fail;
   case -10: //bad_tr_strategy
     PyErr_SetString(PyExc_RuntimeError,
		     "Bad trust region strategy");
     goto fail;
   case -11: // find_beta
     PyErr_SetString(PyExc_RuntimeError,
		     "Error in find_beta");
     goto output;
   case -12: // bad_scaling
     PyErr_SetString(PyExc_RuntimeError,
		     "Bad trust region strategy");
     goto fail;
   case -101: // dogleg_model
     PyErr_SetString(PyExc_RuntimeError,
		     "Bad model for dogleg");
     goto fail;
   case -201: // aint_eig_imag
     PyErr_SetString(PyExc_RuntimeError,
		     "Error in aint_eig_imag");
     goto output;
   case -202: // aint_eig_odd
     PyErr_SetString(PyExc_RuntimeError,
		     "Error in aint_eig_odd");
     goto output;
   case -301: // ms_maxits
     PyErr_SetString(PyExc_RuntimeError,
		     "Max iters reached in More Sorensen");
     goto output;
   case -302: // ms_too_many_shifts
     PyErr_SetString(PyExc_RuntimeError,
		     "Too many shifts in More-sorensen");
     goto output;
   case -303: // ms_no_progress
     PyErr_SetString(PyExc_RuntimeError,
		     "No progress in More Sorensen");
     goto output;
   default: ; // empty statement for language conformatity.
     char errmsg[100];
     sprintf(errmsg, "NLLS_SOLVE failed with unrecognised error code %d\n",
	     inform.status);
     PyErr_SetString(PyExc_RuntimeError, errmsg);
     goto fail;
   }

   /* Free references and return solution */
 output:
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
