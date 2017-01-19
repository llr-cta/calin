/*

   calin/calin_typemaps_matrix3d.i -- Stephen Fegan -- 2017-01-19

   SWIG interface file for common calin typemaps : Eigen::Matrix3d

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay

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
// Typemaps for using Eigen::Matrix3d - these require data be copied once on
// input and again on output (for non-const references)
//
// =============================================================================

%fragment("Calin_Python_to_EigenMat3",
          "header",
          fragment="NumPy_Array_Requirements",
          fragment="NumPy_Backward_Compatibility",
          fragment="NumPy_Macros",
          fragment="NumPy_Utilities")
{
  static bool calin_python_to_eigen_mat3(PyObject* input, Eigen::Matrix3d* mat)
  {
    const int typecode = NPY_DOUBLE;

    if(!_swig_numpy_is_array(input))
      {
        const char* desired_type = typecode_string(typecode);
        const char* actual_type  = pytype_string(input);
        PyErr_Format(PyExc_TypeError,
                     "Array of type '%s' required. A '%s' was given",
                     desired_type,
                     actual_type);
        return false;
      }

    PyArrayObject* in_array = (PyArrayObject*) input;

    if(_swig_numpy_array_numdims(in_array) != 2)
    {
      PyErr_Format(PyExc_TypeError,
                   "Array must have exactly 2 dimensions. "
                   "Given array has %d dimensions",
                   _swig_numpy_array_numdims(in_array));
      return false;
    }

    if(_swig_numpy_array_size(in_array, 0)!=3 or
       _swig_numpy_array_size(in_array, 1)!=3)
    {
      PyErr_Format(PyExc_TypeError,
                   "Matrix must have 3x3 elements. "
                   "Given matrix has %dx%d elements",
                   _swig_numpy_array_size(in_array, 0),
                   _swig_numpy_array_size(in_array, 1));
      return false;
    }

    npy_intp size[2] = { 3, 3 };

    PyArrayObject* out_array = (PyArrayObject*)
        PyArray_New(&PyArray_Type, _swig_numpy_array_numdims(in_array), size, typecode,
                    NULL, mat->data(), 0, NPY_ARRAY_FARRAY, NULL);

    if(out_array == nullptr)return false;

    if(PyArray_CopyInto(out_array, in_array) != 0)
      {
        Py_DECREF(out_array);
        return false;
      }

    Py_DECREF(out_array);
    return true;
  }

} // %fragment("Calin_Python_to_EigenMat3"

%fragment("Calin_EigenMat3_to_Python",
          "header",
          fragment="NumPy_Array_Requirements",
          fragment="NumPy_Backward_Compatibility",
          fragment="NumPy_Macros",
          fragment="NumPy_Utilities")
{

  static bool calin_eigen_mat3_to_python(Eigen::Matrix3d* mat,
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

    npy_intp size[2] = { 3, 3 };
    PyArrayObject* in_array = (PyArrayObject*)
        PyArray_New(&PyArray_Type, 2, size, typecode,
                    NULL, mat->data(), 0, NPY_ARRAY_FARRAY, NULL);
    if(in_array == nullptr)
    {
      return false;
    }

    PyArrayObject* out_array = (PyArrayObject*) output;

    PyArray_Dims dims = { size, 2 };
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

} // %fragment("Calin_EigenMat3_to_Python",

// *************************** const Eigen::Matrix3d& **************************

%typemap(in, fragment="Calin_Python_to_EigenMat3")
     const Eigen::Matrix3d& (Eigen::Matrix3d temp)
{
  // typemap(in) const Eigen::Matrix3d& -- calin_typemaps.i
  $1 = &temp;
  if(!calin_python_to_eigen_mat3($input, $1))SWIG_fail;
}

%typemap(out, fragment="Calin_EigenMat3_to_Python") const Eigen::Matrix3d&
{
  // typemap(out) const Eigen::Matrix3d& -- calin_typemaps.i
  npy_intp size[1] { 3 };
  $result = PyArray_EMPTY(1, size, NPY_DOUBLE, 0);
  if(!$result)SWIG_fail;
  if(!calin_eigen_mat3_to_python($1, $result))SWIG_fail;
}

%typemap(argout) const Eigen::Matrix3d&
{
  // typemap(argout) const Eigen::Matrix3d& -- calin_typemaps.i
  // nothing to see here
}

%typemap(typecheck, precedence=5000) const Eigen::Matrix3d&
{
  // typemap(typecheck) const Eigen::Matrix3d& -- calin_typemaps.i
  $1 = _swig_numpy_is_array($input) ? 1 : 0;
}

// ****************************** Eigen::Matrix3d& *****************************

%typemap(in, fragment="Calin_Python_to_EigenMat3")
     Eigen::Matrix3d& (Eigen::Matrix3d temp)
{
  // typemap(in) Eigen::Matrix3d& -- calin_typemaps.i
  $1 = &temp;
  if(!calin_python_to_eigen_mat3($input, $1))SWIG_fail;
}

%typemap(argout, fragment="Calin_EigenMat3_to_Python") Eigen::Matrix3d&
{
  // typemap(argout) Eigen::Matrix3d& -- calin_typemaps.i
  if(!calin_eigen_mat3_to_python($1, $input))SWIG_fail;
}

%typemap(typecheck, precedence=5000) Eigen::Matrix3d&
{
  // typemap(typecheck) Eigen::Matrix3d& -- calin_typemaps.i
  $1 = _swig_numpy_is_array($input) ? 1 : 0;
}

// ************************** Eigen::Matrix3d &OUTPUT **************************

%typemap(in, numinputs=0) Eigen::Matrix3d &OUTPUT (Eigen::Matrix3d temp)
{
  // typemap(in) Eigen::Matrix3d &OUTPUT -- calin_typemaps.i
  $1 = &temp;
}

%typemap(argout, fragment="Calin_EigenMat3_to_Python") Eigen::Matrix3d &OUTPUT
{
  // typemap(argout) Eigen::Matrix3d &OUTPUT -- calin_typemaps.i
  npy_intp size[1] { $1->size() };
  PyObject* temp_array = PyArray_EMPTY(1, size, NPY_DOUBLE, 0);
  if(!temp_array)SWIG_fail;
  if(!calin_eigen_mat3_to_python($1, temp_array))
  {
    Py_DECREF(temp_array);
    SWIG_fail;
  }
  $result = SWIG_Python_AppendOutput($result, temp_array);
}

// ************************** Eigen::Matrix3d &INOUT ***************************

%typemap(in, fragment="Calin_Python_to_EigenMat3")
  Eigen::Matrix3d &INOUT (Eigen::Matrix3d temp)
{
  // typemap(in) const Eigen::Matrix3d &INOUT -- calin_typemaps.i
  $1 = &temp;
  if(!calin_python_to_eigen_mat3($input, $1))SWIG_fail;
}

%typemap(argout, fragment="Calin_EigenMat3_to_Python") Eigen::Matrix3d &INOUT
{
  // typemap(argout) Eigen::Matrix3d &INOUT -- calin_typemaps.i
  npy_intp size[1] { $1->size() };
  PyObject* temp_array = PyArray_EMPTY(1, size, NPY_DOUBLE, 0);
  if(!temp_array)SWIG_fail;
  if(!calin_eigen_mat3_to_python($1, temp_array))
  {
    Py_DECREF(temp_array);
    SWIG_fail;
  }
  $result = SWIG_Python_AppendOutput($result, temp_array);
}

// ****************************** Eigen::Matrix3d ******************************

%typemap(out, fragment="Calin_EigenMat3_to_Python") Eigen::Matrix3d
{
  // typemap(out) Eigen::Matrix3d -- calin_typemaps.i
  npy_intp size[1] { $1.size() };
  $result = PyArray_EMPTY(1, size, NPY_DOUBLE, 0);
  if(!$result)SWIG_fail;
  if(!calin_eigen_mat3_to_python(&$1, $result))SWIG_fail;
}

%typemap(typecheck, precedence=5000) Eigen::Matrix3d
{
  // typemap(typecheck) Eigen::Matrix3d -- calin_typemaps.i
  $1 = _swig_numpy_is_array($input) ? 1 : 0;
}
