// Copyright (c) 2016, The Science and Technology Facilities Council (STFC)
// All rights reserved.
// Copyright (C) 2024 Advanced Micro Devices, Inc. All rights reserved.
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

#include "ral_nlls.h"

struct ral_nlls_state {
  PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct ral_nlls_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct ral_nlls_state _state;
#endif

struct callback_data {
  PyObject* f;
  PyObject* J;
  PyObject* Hr;
  PyObject* Hp;
  PyObject* params;
};

///
/// get the argument list of a
///
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

///
/// the eval_f subroutine
///
static
int eval_f(int n, int m, void *params, const double *x, double *f) {
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

///
/// the eval_J subroutine
///
static
int eval_J(int n, int m, void *params, const double *x, double *J) {
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

///
/// the eval_Hr subroutine
///
static
int eval_Hr(int n, int m, void *params, const double *x, const double *r,
	    double *Hr) {
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

///
/// the eval_Hp subroutine
///
static
int eval_Hp(int n, int m, const double *x, const double *y,
	    double *Hp, void *params) {
   // Recover our datatype
   const struct callback_data *data = (struct callback_data*) params;

   // Copy x into Python array
   npy_intp xdim[] = {n};
   PyArrayObject* xpy = (PyArrayObject*) PyArray_SimpleNew(1, xdim, NPY_DOUBLE);
   double* xval = (double*) PyArray_DATA(xpy);
   for(int i=0; i<n; ++i)
      xval[i] = x[i];

   // Copy y into Python array
   npy_intp ydim[] = {n};
   PyArrayObject* ypy = (PyArrayObject*) PyArray_SimpleNew(1, ydim, NPY_DOUBLE);
   double* yval = (double*) PyArray_DATA(ypy);
   for(int i=0; i<n; ++i)
      yval[i] = y[i];

   // Call routine
   PyObject *arglist = build_arglist(2, data->params);
   PyTuple_SET_ITEM(arglist, 0, (PyObject*) xpy); Py_INCREF(xpy);
   PyTuple_SET_ITEM(arglist, 1, (PyObject*) ypy); Py_INCREF(ypy);
   PyObject *result = PyObject_CallObject(data->Hp, arglist);
   Py_DECREF(arglist);
   Py_DECREF(xpy);
   Py_DECREF(ypy);
   if(!result) return -1;

   // Extract result
   PyArrayObject* Hparray = (PyArrayObject*) PyArray_FROM_OTF(result, NPY_FLOAT64, NPY_ARRAY_IN_FARRAY);
   if(Hparray == NULL) {
      PyErr_SetString(PyExc_RuntimeError, "Error extracting array from Hp call");
      Py_DECREF(result);
      return -1;
   }
   if(PyArray_NDIM(Hparray) != 2) {
      PyErr_SetString(PyExc_RuntimeError, "Hp() must return rank-2 array");
      Py_DECREF(Hparray);
      Py_DECREF(result);
      return -2;
   }
   const double *Hpval = (double*) PyArray_DATA(Hparray);
   for(int i=0; i<n*m; ++i) {
      Hp[i] = Hpval[i];
   }
   Py_DECREF(Hparray);
   Py_DECREF(result);

   return 0; // Success
}

///
/// set the options up...
/// pick up the defaults from the C code, and update any
/// passed via Python
///
static
bool set_opts(struct ral_nlls_options *options, PyObject *pyoptions) {
   ral_nlls_default_options_d(options);

   if(!pyoptions) return true; // Just use defaults

   #if PY_MAJOR_VERSION >= 3 || (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION > 6)
   #define PyInt_AsLong   PyLong_AsLong
   #define PyInt_FromLong   PyLong_FromLong
   #endif
   PyObject *key, *value;
   Py_ssize_t pos = 0;
   const char* key_name;

   while(PyDict_Next(pyoptions, &pos, &key, &value)) {
     if (PyUnicode_Check(key)){
       // change key from Unicode to ascii
       key = PyUnicode_AsASCIIString(key);
     }
     key_name = PyBytes_AsString(key);
     if(!key_name) {
       PyErr_SetString(PyExc_RuntimeError, "The keys in the options dictionary must be strings.");
       return false;
     }
     if(strcmp(key_name, "out")==0) {
       long v = PyInt_AsLong(value);
       if(v==-1 && PyErr_Occurred()) {
	 PyErr_SetString(PyExc_RuntimeError, "options['out'] must be an integer.");
	 return false;
       }
       options->out = (int) v;
       continue;
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
      if(strcmp(key_name, "print_options")==0) {
  	    int vint = PyObject_IsTrue(value); // 1 if true, 0 otherwise
		    if (vint == 1) {
	        options->print_options=true;
      	} else if (vint == 0) {
      	  options->print_options=false;
      	} else {
    	    PyErr_SetString(PyExc_RuntimeError, "options['print_options'] must be a bool.");
      	  return false;
	      }
      	continue;
      }
      if(strcmp(key_name, "print_header")==0) {
         long v = PyInt_AsLong(value);
         if(v==-1 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['print_header'] must be an integer.");
            return false;
         }
         options->print_header = (int) v;
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
      if(strcmp(key_name, "type_of_method")==0) {
	long v = PyInt_AsLong(value);
	if(v==-1 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['type_of_method'] must be an integer.");
	  return false;
	}
	options->type_of_method = (int) v;
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
      // bool: allow_fallback_method

      if(strcmp(key_name, "allow_fallback_method")==0) {
	int vint = PyObject_IsTrue(value); // 1 if true, 0 otherwise
	printf("%d\n",vint);
	if (vint == 1){
	  options->allow_fallback_method=true;
	}else if (vint == 0){
	  options->allow_fallback_method=false;
	}else{
	  PyErr_SetString(PyExc_RuntimeError, "options['allow_fallback_method'] must be a bool.");
	  return false;
	}
	continue;
      }

      if(strcmp(key_name, "lls_solver")==0) {
	long v = PyInt_AsLong(value);
	if(v==-1 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['lls_solver'] must be an integer.");
	  return false;
	}
	options->lls_solver = (int) v;
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


      if(strcmp(key_name, "stop_f_absolute")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['stop_f_absolute'] must be a float.");
	  return false;
	}
	options->stop_f_absolute = v;
	continue;
      }

      if(strcmp(key_name, "stop_f_relative")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['stop_f_relative'] must be a float.");
	  return false;
	}
	options->stop_f_relative = v;
	continue;
      }

      if(strcmp(key_name, "stop_s")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['stop_s'] must be a float.");
	  return false;
	}
	options->stop_s = v;
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

      if(strcmp(key_name, "base_regularization")==0) {
	long v = PyInt_AsLong(value);
	if(v==-1 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['base_regularization'] must be an integer.");
	  return false;
	}
	options->base_regularization = (int) v;
	continue;
      }


      if(strcmp(key_name, "regularization")==0) {
	long v = PyInt_AsLong(value);
	if(v==-1 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['regularization'] must be an integer.");
	  return false;
	}
	options->regularization = (int) v;
	continue;
      }


      if(strcmp(key_name, "regularization_term")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['regularization_term'] must be a float.");
	  return false;
	}
	options->regularization_term = v;
	continue;
      }


      if(strcmp(key_name, "regularization_power")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['regularization_power'] must be a float.");
	  return false;
	}
	options->regularization_power = v;
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
            PyErr_SetString(PyExc_RuntimeError, "options['tr_update_strategy'] must be an integer.");
            return false;
	}
	options->tr_update_strategy = (int) v;
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
      if(strcmp(key_name, "exact_second_derivatives")==0) {
	int vint = PyObject_IsTrue(value); // 1 if true, 0 otherwise
	if (vint == 1){
	  options->exact_second_derivatives=true;
	}else if (vint == 0){
	  options->exact_second_derivatives=false;
	}else{
	  PyErr_SetString(PyExc_RuntimeError, "options['exact_second_derivatives'] must be a bool.");
	  return false;
	}
	continue;
      }

      if(strcmp(key_name, "subproblem_eig_fact")==0) {
	int vint = PyObject_IsTrue(value); // 1 if true, 0 otherwise
	printf("%d\n",vint);
	if (vint == 1){
	  options->subproblem_eig_fact=true;
	}else if (vint == 0){
	  options->subproblem_eig_fact=false;
	}else{
	  PyErr_SetString(PyExc_RuntimeError, "options['subproblem_eig_fact'] must be a bool.");
	  return false;
	}
	continue;
      }

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

      if(strcmp(key_name, "scale_trim_max")==0) {
	int vint = PyObject_IsTrue(value); // 1 if true, 0 otherwise

	if (vint == 1){
	  options->scale_trim_max=true;
	}else if (vint == 0){
	  options->scale_trim_max=false;
	}else{
	  PyErr_SetString(PyExc_RuntimeError, "options['scale_trim_max'] must be a bool.");
	  return false;
	}
	continue;
      }

      if(strcmp(key_name, "scale_trim_min")==0) {
	int vint = PyObject_IsTrue(value); // 1 if true, 0 otherwise

	if (vint == 1){
	  options->scale_trim_min=true;
	}else if (vint == 0){
	  options->scale_trim_min=false;
	}else{
	  PyErr_SetString(PyExc_RuntimeError, "options['scale_trim_min'] must be a bool.");
	  return false;
	}
	continue;
      }

      if(strcmp(key_name, "scale_require_increase")==0) {
	int vint = PyObject_IsTrue(value); // 1 if true, 0 otherwise

	if (vint == 1){
	  options->scale_require_increase=true;
	}else if (vint == 0){
	  options->scale_require_increase=false;
	}else{
	  PyErr_SetString(PyExc_RuntimeError, "options['scale_require_increase'] must be a bool.");
	  return false;
	}
	continue;
      }


      if(strcmp(key_name, "setup_workspaces")==0) {
	int vint = PyObject_IsTrue(value); // 1 if true, 0 otherwise

	if (vint == 1){
	  options->setup_workspaces=true;
	}else if (vint == 0){
	  options->setup_workspaces=false;
	}else{
	  PyErr_SetString(PyExc_RuntimeError, "options['setup_workspaces'] must be a bool.");
	  return false;
	}
	continue;
      }


      if(strcmp(key_name, "remove_workspaces")==0) {
	int vint = PyObject_IsTrue(value); // 1 if true, 0 otherwise

	if (vint == 1){
	  options->remove_workspaces=true;
	}else if (vint == 0){
	  options->remove_workspaces=false;
	}else{
	  PyErr_SetString(PyExc_RuntimeError, "options['remove_workspaces'] must be a bool.");
	  return false;
	}
	continue;
      }

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
            PyErr_SetString(PyExc_RuntimeError, "options['hybrid_switch_its'] must be an integer.");
            return false;
	}
	options->hybrid_switch_its = (int) v;
	continue;
      }

      if(strcmp(key_name, "reg_order")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['reg_order'] must be a float.");
	  return false;
	}
	options->reg_order = v;
	continue;
      }

      if(strcmp(key_name, "inner_method")==0) {
	long v = PyInt_AsLong(value);
	if(v==-1 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['inner_method'] must be an integer.");
	  return false;
	}
	options->inner_method = (int) v;
	continue;
      }


      if(strcmp(key_name, "output_progress_vectors")==0) {
	int vint = PyObject_IsTrue(value); // 1 if true, 0 otherwise

	if (vint == 1){
	  options->output_progress_vectors=true;
	}else if (vint == 0){
	  options->output_progress_vectors=false;
	}else{
	  PyErr_SetString(PyExc_RuntimeError, "options['output_progress_vectors'] must be a bool.");
	  return false;
	}
	continue;
      }

      // bool: Fortran_Jacobian
      if(strcmp(key_name, "Fortran_Jacobian")==0) {
	int vint = PyObject_IsTrue(value); // 1 if true, 0 otherwise
	printf("%d\n",vint);
	if (vint == 1){
	  options->Fortran_Jacobian=true;
	}else if (vint == 0){
	  options->Fortran_Jacobian=false;
	}else{
	  PyErr_SetString(
			  PyExc_RuntimeError,
			  "options['Fortran_Jacobian'] must be a bool."
			  );
	  return false;
	}
	continue;
      }

      /* integer(c_int) :: box_nFref_max */
      if(strcmp(key_name, "box_nFref_max")==0) {
	long v = PyInt_AsLong(value);
	if(v==-1 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['box_nFref_max'] must be an integer.");
	  return false;
	}
	options->box_nFref_max = (int) v;
	continue;
      }


      /* 	real(wp) :: box_gamma */
      if(strcmp(key_name, "box_gamma")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['box_gamma'] must be a float.");
            return false;
	}
	options->box_gamma = v;
	continue;
      }
      /* 	real(wp) :: box_decmin */
      if(strcmp(key_name, "box_decmin")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['box_decmin'] must be a float.");
            return false;
	}
	options->box_decmin = v;
	continue;
      }

      /* 	real(wp) :: box_bigbnd */
      if(strcmp(key_name, "box_bigbnd")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['box_bigbnd'] must be a float.");
	  return false;
	}
	options->box_bigbnd = v;
	continue;
      }
      /* 	real(wp) :: box_wolfe_descent */
      if(strcmp(key_name, "box_wolfe_descent")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['box_wolfe_descent'] must be a float.");
	  return false;
	}
	options->box_wolfe_descent = v;
	continue;
      }
      /* 	real(wp) :: box_wolfe_curvature */
      if(strcmp(key_name, "box_wolfe_curvature")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['box_wolfe_curvature'] must be a float.");
	  return false;
	}
	options->box_wolfe_curvature = v;
	continue;
      }

      /* 	real(wp) :: box_kanzow_power */
      if(strcmp(key_name, "box_kanzow_power")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['box_kanzow_power'] must be a float.");
	  return false;
         }
	options->box_kanzow_power = v;
	continue;
      }
      /* 	real(wp) :: box_kanzow_descent */
      if(strcmp(key_name, "box_kanzow_descent")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['box_kanzow_descent'] must be a float.");
	  return false;
	}
	options->box_kanzow_descent = v;
	continue;
      }

      /* 	real(wp) :: box_quad_model_descent */
      if(strcmp(key_name, "box_quad_model_descent")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['box_quad_model_descent'] must be a float.");
	  return false;
         }
	options->box_quad_model_descent = v;
	continue;
      }

      /* 	Logical(c_bool):: box_tr_test_step */
      if(strcmp(key_name, "box_tr_test_step")==0) {
	int vint = PyObject_IsTrue(value); // 1 if true, 0 otherwise
	printf("%d\n",vint);
	if (vint == 1){
	  options->box_tr_test_step=true;
	}else if (vint == 0){
	  options->box_tr_test_step=false;
	}else{
	  PyErr_SetString(PyExc_RuntimeError, "options['box_tr_test_step'] must be a bool.");
	  return false;
	}
	continue;
      }

      /* 	Logical(c_bool):: box_wolfe_test_step */
      if(strcmp(key_name, "box_wolfe_test_step")==0) {
	int vint = PyObject_IsTrue(value); // 1 if true, 0 otherwise
	printf("%d\n",vint);
	if (vint == 1){
	  options->box_wolfe_test_step=true;
	}else if (vint == 0){
	  options->box_wolfe_test_step=false;
	}else{
	  PyErr_SetString(PyExc_RuntimeError, "options['box_wolfe_test_step'] must be a bool.");
	  return false;
	}
	continue;
      }

      /* 	real(wp) :: box_tau_descent */
      if(strcmp(key_name, "box_tau_descent")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['box_tau_descent'] must be a float.");
            return false;
	}
	options->box_tau_descent = v;
	continue;
      }

      /* 	integer(c_int):: box_max_ntrfail */
      if(strcmp(key_name, "box_max_ntrfail")==0) {
	long v = PyInt_AsLong(value);
	if(v==-1 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['box_max_ntrfail'] must be an integer.");
	  return false;
	}
	options->box_max_ntrfail = (int) v;
	continue;
      }

      /* 	integer(c_int):: box_quad_match */
      if(strcmp(key_name, "box_quad_match")==0) {
	long v = PyInt_AsLong(value);
	if(v==-1 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['box_quad_match'] must be an integer.");
	  return false;
	}
	options->box_quad_match = (int) v;
	continue;
      }

      /* 	real(wp) :: box_alpha_scale */
      if(strcmp(key_name, "box_alpha_scale")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['box_alpha_scale'] must be a float.");
	  return false;
         }
	options->box_alpha_scale = v;
	continue;
      }

      /* 	real(wp) :: box_Delta_scale */
      if(strcmp(key_name, "box_Delta_scale")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['box_Delta_scale'] must be a float.");
            return false;
	}
	options->box_Delta_scale = v;
	continue;
      }

      /* 	real(wp) :: box_tau_min */
      if(strcmp(key_name, "box_tau_min")==0) {
	double v = PyFloat_AsDouble(value);
	if(v==-1.0 && PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['box_tau_min'] must be a float.");
            return false;
	}
         options->box_tau_min = v;
         continue;
      }

      /* 	integer(c_int):: box_ls_step_maxit */
      if(strcmp(key_name, "box_ls_step_maxit")==0) {
	long v = PyInt_AsLong(value);
	if(v==-1 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['box_ls_step_maxit'] must be an integer.");
	  return false;
	}
	options->box_ls_step_maxit = (int) v;
	continue;
      }

      /* 	integer(c_int):: box_linesearch_type */
      if(strcmp(key_name, "box_linesearch_type")==0) {
	long v = PyInt_AsLong(value);
	if(v==-1 && PyErr_Occurred()) {
	  PyErr_SetString(PyExc_RuntimeError, "options['box_linesearch_type'] must be an integer.");
	  return false;
	}
	options->box_linesearch_type = (int) v;
	continue;
      }

      /* 	real(wp) :: fd_step */
      if(strcmp(key_name, "fd_step")==0) {
	     double v = PyFloat_AsDouble(value);
	     if(v<=0 || PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['fd_step'] must be a positve float.");
            return false;
	     }
         options->fd_step = v;
         continue;
      }

     /* Integer       :: check_derivatives = 0 */
      if(strcmp(key_name, "check_derivatives")==0) {
	    long v = PyInt_AsLong(value);
	       if(v==-1 && PyErr_Occurred()) {
	          PyErr_SetString(PyExc_RuntimeError, "options['check_derivatives'] must be an integer.");
	          return false;
	       }
	    options->check_derivatives = (int) v;
	    continue;
      }

     /* Real(Kind=wp) :: derivative_test_tol = 1.0e-4_wp */
      if(strcmp(key_name, "derivative_test_tol")==0) {
	     double v = PyFloat_AsDouble(value);
	     if(v<=0 || PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, "options['derivative_test_tol'] must be a positve float.");
            return false;
	     }
         options->derivative_test_tol = v;
         continue;
      }

      // If we reach this point, unrecognised option
      char errmsg[200];
      snprintf(errmsg, 200, "Bad key options['%s']\n", key_name);
      PyErr_SetString(PyExc_RuntimeError, errmsg);
      return false;
   }

   return true; // success
}

//
// Take the info struct from C, and turn it into
// a python dictionary
//

static PyObject*
make_info_dict(const struct ral_nlls_inform *inform) {
   PyObject *pyinfo = PyDict_New();

   PyDict_SetItemString(pyinfo, "iter", PyInt_FromLong(inform->iter));
   PyDict_SetItemString(pyinfo, "inner_iter", PyInt_FromLong(inform->inner_iter));
   // bool:: inner_iter_success?
   PyDict_SetItemString(pyinfo, "f_eval", PyInt_FromLong(inform->f_eval));
   PyDict_SetItemString(pyinfo, "g_eval", PyInt_FromLong(inform->g_eval));
   PyDict_SetItemString(pyinfo, "h_eval", PyInt_FromLong(inform->h_eval));
   PyDict_SetItemString(pyinfo, "convergence_normf",PyInt_FromLong(inform->convergence_normf));
   PyDict_SetItemString(pyinfo, "convergence_normg",PyInt_FromLong(inform->convergence_normg));
   PyDict_SetItemString(pyinfo, "convergence_norms",PyInt_FromLong(inform->convergence_norms));
   //   PyDict_SetItemString(pyinfo, "resinf", PyFloat_FromDouble(inform->resinf));
   //   PyDict_SetItemString(pyinfo, "gradinf", PyFloat_FromDouble(inform->gradinf));
   PyDict_SetItemString(pyinfo, "obj", PyFloat_FromDouble(inform->obj));
   PyDict_SetItemString(pyinfo, "norm_g", PyFloat_FromDouble(inform->norm_g));
   PyDict_SetItemString(pyinfo, "scaled_g",PyFloat_FromDouble(inform->scaled_g));
   //   PyDict_SetItemString(pyinfo, "step",PyFloat_FromDouble(inform->step));
   PyDict_SetItemString(pyinfo, "ls_step_iter", PyInt_FromLong(inform->ls_step_iter));
   PyDict_SetItemString(pyinfo, "f_eval_ls", PyInt_FromLong(inform->f_eval_ls));
   PyDict_SetItemString(pyinfo, "g_eval_ls", PyInt_FromLong(inform->g_eval_ls));
   PyDict_SetItemString(pyinfo, "pg_step_iter", PyInt_FromLong(inform->pg_step_iter));
   PyDict_SetItemString(pyinfo, "f_eval_pg", PyInt_FromLong(inform->f_eval_pg));
   PyDict_SetItemString(pyinfo, "g_eval_pg", PyInt_FromLong(inform->g_eval_pg));
   PyDict_SetItemString(pyinfo, "fd_f_eval", PyInt_FromLong(inform->fd_f_eval));

   return pyinfo;
}


///
/// call the solve routine
///


/*
 * x = ral_nlls.solve(x0, f, J=None, Hr=None, params=(), options={}
 *                    weights=None, Hp=None,
 *                    lower_bounds=None, upper_bounds=None)
 */
static PyObject*
ral_nlls_solve(PyObject* self, PyObject* args, PyObject* keywds)
{
  PyObject *x0ptr=NULL, *options_ptr=NULL, *weightsptr=NULL;
  PyObject *lower_boundsptr=NULL, *upper_boundsptr=NULL;
  PyObject *arglist=NULL, *result=NULL;
  PyArrayObject *x0=NULL, *f=NULL, *x=NULL;
  PyArrayObject *weights=NULL, *lower_bounds=NULL, *upper_bounds=NULL;

  double *weights_val, *lower_bounds_val, *upper_bounds_val;

  struct callback_data data;
  data.f = NULL;
  data.J = NULL;
  data.Hr = NULL;
  data.params = NULL;
  data.Hp = NULL;

  static char *kwlist[] = {"x0",
			   "r",
			   "J",
			   "Hr",
			   "params",
			   "options",
			   "weights",
			   "Hp",
			   "lower_bounds",
			   "upper_bounds",
			    NULL};
  if(!PyArg_ParseTupleAndKeywords(args,
				  keywds,
				  "OOO|OOO!OOOO",
				  kwlist,
				  &x0ptr,
				  &data.f,
				  &data.J,
				  &data.Hr,
				  &data.params,
				  &PyDict_Type,
				  &options_ptr,
				  &weightsptr,
				  &data.Hp,
				  &lower_boundsptr,
				  &upper_boundsptr
				  )
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

  /* weights */
  if (weightsptr == NULL) {
    weights_val = NULL;
  }
  else{
    weights = (PyArrayObject*) PyArray_FROM_OTF(weightsptr,
						NPY_FLOAT64,
						NPY_ARRAY_IN_ARRAY);
    if(PyArray_NDIM(weights) != 1) {
      PyErr_SetString(PyExc_RuntimeError, "weights must be a rank-1 array");
      goto fail;
    }
    weights_val = (double*) PyArray_DATA(weights);
  }


  /* lower_bounds */
  if (lower_boundsptr == NULL) {
    lower_bounds_val = NULL;
  }
  else{
    lower_bounds = (PyArrayObject*) PyArray_FROM_OTF(lower_boundsptr,
						NPY_FLOAT64,
						NPY_ARRAY_IN_ARRAY);
    if(PyArray_NDIM(lower_bounds) != 1) {
      PyErr_SetString(PyExc_RuntimeError, "lower_bounds must be a rank-1 array");
      goto fail;
    }
    lower_bounds_val = (double*) PyArray_DATA(lower_bounds);
  }

  /* upper_bounds */
  if (upper_boundsptr == NULL) {
    upper_bounds_val = NULL;
  }
  else{
    upper_bounds = (PyArrayObject*) PyArray_FROM_OTF(upper_boundsptr,
						NPY_FLOAT64,
						NPY_ARRAY_IN_ARRAY);
    if(PyArray_NDIM(upper_bounds) != 1) {
      PyErr_SetString(PyExc_RuntimeError, "upper_bounds must be a rank-1 array");
      goto fail;
    }
    upper_bounds_val = (double*) PyArray_DATA(upper_bounds);
  }

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
   if(data.Hp == NULL){
     nlls_solve_d(n, m, xval, eval_f, eval_J, eval_Hr, &data, &options, &inform,
		  weights_val, NULL, lower_bounds_val, upper_bounds_val);
   }
   else{
   nlls_solve_d(n, m, xval, eval_f, eval_J, eval_Hr, &data, &options, &inform,
   		  weights_val, eval_Hp, lower_bounds_val, upper_bounds_val);
   }
   switch(inform.status) {
   case 0: // Clean exit
     break;
   case -1: // Exceeded max itr
     PyErr_WarnEx(PyExc_RuntimeWarning,
		  "Exceeded maximum number of iterations",2);
     goto output;
   case -8: // No progress being made
     PyErr_WarnEx(PyExc_RuntimeWarning,
		  "No progress being made towards the solution",2);
     goto output;
   default:
     PyErr_SetString(PyExc_RuntimeError,
		     inform.error_message); // print out the error message passed
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
    "   x = ral_nlls.solve(x0, r, J, Hr=None, params=None, options={}, \n"
    "                      weights=None, Hp=None,\n"
    "                      lower_bounds=None, upper_bounds=None)"
   },
   {NULL, NULL, 0, NULL} /* Sentinel */
};


#if PY_MAJOR_VERSION >= 3

static int ral_nlls_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int ral_nlls_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef ral_nllsdef = {
        PyModuleDef_HEAD_INIT,
        "ral_nlls",
        NULL,
        sizeof(struct ral_nlls_state),
        RalNllsMethods,
        NULL,
        ral_nlls_traverse,
        ral_nlls_clear,
        NULL
};

#define INITERROR return NULL

PyMODINIT_FUNC
PyInit_ral_nlls(void)

#else
#define INITERROR return

void
initral_nlls(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
PyObject *ral_nlls = PyModule_Create(&ral_nllsdef);
#else
PyObject *ral_nlls = Py_InitModule("ral_nlls", RalNllsMethods);
#endif
if (ral_nlls == NULL)
  INITERROR;
struct ral_nlls_state *st = GETSTATE(ral_nlls);
st->error = PyErr_NewException("ral_nlls.Error", NULL, NULL);
if (st->error == NULL) {
   Py_DECREF(ral_nlls);
   INITERROR;
}
import_array();
#if PY_MAJOR_VERSION >= 3
return ral_nlls;
#endif
}
