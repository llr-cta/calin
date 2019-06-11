/*

   calin/calin_typemaps_vectorxd.i -- Stephen Fegan -- 2017-01-19

   SWIG interface file for common calin typemaps : Eigen::VectorXd

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, CNRS, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

// =============================================================================
//
// Typemaps for using Eigen::VectorXd - these require data be copied once on
// input and again on output (for non-const references)
//
// =============================================================================

%fragment("Calin_Python_to_EigenVec",
          "header",
          fragment="NumPy_Array_Requirements",
          fragment="NumPy_Backward_Compatibility",
          fragment="NumPy_Macros",
          fragment="NumPy_Utilities")
{
  static bool calin_python_to_eigen_vec(PyObject* input, Eigen::VectorXd* vec)
  {
    const int typecode = NPY_DOUBLE;

    if(!_swig_numpy_is_array(input))
      {
        const char* desired_type = typecode_string(typecode);
        const char* actual_type  = pytype_string(input);
        PyErr_Format(PyExc_TypeError,
                     "Array of type '%s' required.  A '%s' was given",
                     desired_type,
                     actual_type);
        return false;
      }

    PyArrayObject* in_array = (PyArrayObject*) input;

    if(_swig_numpy_array_numdims(in_array) > 1)
    {
      PyErr_Format(PyExc_TypeError,
                   "Array must have 1 dimension. "
                   "Given array has %d dimensions",
                   _swig_numpy_array_numdims(in_array));
      return false;
    }

    if(_swig_numpy_array_numdims(in_array)==0 or _swig_numpy_array_size(in_array, 0)==0)
    {
      *vec = Eigen::VectorXd();
      return true;
    }

    npy_intp size[1] = { _swig_numpy_array_size(in_array, 0) };
    vec->resize(size[0]);

    PyArrayObject* out_array = (PyArrayObject*)
        PyArray_SimpleNewFromData(1, size, typecode, vec->data());
    if(out_array == nullptr)return false;

    if(PyArray_CopyInto(out_array, in_array) != 0)
      {
        Py_DECREF(out_array);
        return false;
      }

    Py_DECREF(out_array);
    return true;
  }

} // fragment("Calin_Python_to_EigenVec"

%fragment("Calin_EigenVec_to_Python",
          "header",
          fragment="NumPy_Array_Requirements",
          fragment="NumPy_Backward_Compatibility",
          fragment="NumPy_Macros",
          fragment="NumPy_Utilities")
{
  static bool calin_eigen_vec_to_python(Eigen::VectorXd* vec,
                                        PyObject* output)
  {
    const int typecode = NPY_DOUBLE;

    if(!_swig_numpy_is_array(output))
      {
        const char* desired_type = typecode_string(typecode);
        const char* actual_type  = pytype_string(output);
        PyErr_Format(PyExc_TypeError,
                     "Array of type '%s' required.  A '%s' was given",
                     desired_type,
                     actual_type);
        return false;
      }

    npy_intp size[1] = { vec->size() };
    PyArrayObject* in_array = (PyArrayObject*)
        PyArray_SimpleNewFromData(1, size, typecode, vec->data());
    if(in_array == nullptr)
      {
        return false;
      }

    PyArrayObject* out_array = (PyArrayObject*) output;

    PyArray_Dims dims = { size, 1 };
    if(PyArray_Resize(out_array, &dims, 0, NPY_ANYORDER) == nullptr)
      {
        // Do we need to call Py_DECREF on returned val??
        Py_DECREF(in_array);
        return false;
      }

    if(PyArray_CopyInto(out_array, in_array) != 0)
      {
        Py_DECREF(in_array);
        return false;
      }

    Py_DECREF(in_array);
    return true;
  }

} // fragment("Calin_EigenVec_to_Python"

// *************************** const Eigen::VectorXd& **************************

%typemap(in, fragment="Calin_Python_to_EigenVec")
     const Eigen::VectorXd& (Eigen::VectorXd temp)
{
  // typemap(in) const Eigen::VectorXd& -- calin_typemap_vectorxd.i
  $1 = &temp;
  if(!calin_python_to_eigen_vec($input, $1))SWIG_fail;
}

%typemap(out, fragment="Calin_EigenVec_to_Python") const Eigen::VectorXd&
{
  // typemap(out) const Eigen::VectorXd& -- calin_typemap_vectorxd.i
  npy_intp size[1] { $1->size() };
  $result = PyArray_EMPTY(1, size, NPY_DOUBLE, 0);
  if(!$result)SWIG_fail;
  if(!calin_eigen_vec_to_python($1, $result))SWIG_fail;
}

%typemap(argout) const Eigen::VectorXd&
{
  // typemap(argout) const Eigen::VectorXd& -- calin_typemap_vectorxd.i
  // nothing to see here
}

%typemap(typecheck, precedence=5000) const Eigen::VectorXd&
{
  // typemap(typecheck) const Eigen::VectorXd& -- calin_typemap_vectorxd.i
  $1 = _swig_numpy_is_array($input) ? 1 : 0;
}

// ****************************** Eigen::VectorXd& *****************************

%typemap(in, fragment="Calin_Python_to_EigenVec")
     Eigen::VectorXd& (Eigen::VectorXd temp)
{
  // typemap(in) Eigen::VectorXd& -- calin_typemap_vectorxd.i
  $1 = &temp;
  if(!calin_python_to_eigen_vec($input, $1))SWIG_fail;
}

%typemap(argout, fragment="Calin_EigenVec_to_Python") Eigen::VectorXd&
{
  // typemap(argout) Eigen::VectorXd& -- calin_typemap_vectorxd.i
  if(!calin_eigen_vec_to_python($1, $input))SWIG_fail;
}

%typemap(typecheck, precedence=5000) Eigen::VectorXd&
{
  // typemap(typecheck) Eigen::VectorXd& -- calin_typemap_vectorxd.i
  $1 = _swig_numpy_is_array($input) ? 1 : 0;
}

// ************************** Eigen::VectorXd &OUTPUT **************************

%typemap(in, numinputs=0) Eigen::VectorXd &OUTPUT (Eigen::VectorXd temp)
{
  // typemap(in) Eigen::VectorXd &OUTPUT -- calin_typemap_vectorxd.i
  $1 = &temp;
}

%typemap(argout, fragment="Calin_EigenVec_to_Python") Eigen::VectorXd &OUTPUT
{
  // typemap(argout) Eigen::VectorXd &OUTPUT -- calin_typemap_vectorxd.i
  npy_intp size[1] { $1->size() };
  PyObject* temp_array = PyArray_EMPTY(1, size, NPY_DOUBLE, 0);
  if(!temp_array)SWIG_fail;
  if(!calin_eigen_vec_to_python($1, temp_array))
  {
    Py_DECREF(temp_array);
    SWIG_fail;
  }
  $result = SWIG_Python_AppendOutput($result, temp_array);
}

// ************************** Eigen::VectorXd &INOUT ***************************

%typemap(in, fragment="Calin_Python_to_EigenVec")
  Eigen::VectorXd &INOUT (Eigen::VectorXd temp)
{
  // typemap(in) const Eigen::VectorXd &INOUT -- calin_typemap_vectorxd.i
  $1 = &temp;
  if(!calin_python_to_eigen_vec($input, $1))SWIG_fail;
}

%typemap(argout, fragment="Calin_EigenVec_to_Python") Eigen::VectorXd &INOUT
{
  // typemap(argout) Eigen::Vector3d &INOUT -- calin_typemap_vectorxd.i
  npy_intp size[1] { $1->size() };
  PyObject* temp_array = PyArray_EMPTY(1, size, NPY_DOUBLE, 0);
  if(!temp_array)SWIG_fail;
  if(!calin_eigen_vec_to_python($1, temp_array))
  {
    Py_DECREF(temp_array);
    SWIG_fail;
  }
  $result = SWIG_Python_AppendOutput($result, temp_array);
}

// ****************************** Eigen::VectorXd ******************************

%typemap(out, fragment="Calin_EigenVec_to_Python") Eigen::VectorXd
{
  // typemap(out) Eigen::VectorXd -- calin_typemap_vectorxd.i
  npy_intp size[1] { $1.size() };
  $result = PyArray_EMPTY(1, size, NPY_DOUBLE, 0);
  if(!$result)SWIG_fail;
  if(!calin_eigen_vec_to_python(&$1, $result))SWIG_fail;
}

%typemap(typecheck, precedence=5000) Eigen::VectorXd
{
  // typemap(typecheck) Eigen::VectorXd -- calin_typemap_vectorxd.i
  $1 = _swig_numpy_is_array($input) ? 1 : 0;
}
