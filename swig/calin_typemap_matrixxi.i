/*

   calin/calin_typemaps_matrixxi.i -- Stephen Fegan -- 2022-10-14

   SWIG interface file for common calin typemaps : Eigen::MatrixXi

   Copyright 2022, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

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
// Typemaps for using Eigen::MatrixXi - these require data be copied once on
// input and again on output (for non-const references)
//
// =============================================================================

%fragment("Calin_Python_to_EigenIntMat",
          "header",
          fragment="NumPy_Array_Requirements",
          fragment="NumPy_Backward_Compatibility",
          fragment="NumPy_Macros",
          fragment="NumPy_Utilities")
{
  static bool calin_python_to_eigen_int_mat(PyObject* input, Eigen::MatrixXi* mat)
  {
    const int typecode = NPY_INT;

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

    if(_swig_numpy_array_numdims(in_array) > 2)
    {
      PyErr_Format(PyExc_TypeError,
                   "Array must have at most 2 dimensions. "
                   "Given array has %d dimensions",
                   _swig_numpy_array_numdims(in_array));
      return false;
    }

    if(_swig_numpy_array_numdims(in_array)==0 or _swig_numpy_array_size(in_array, 0)==0 or
       (_swig_numpy_array_numdims(in_array)==2 and _swig_numpy_array_size(in_array, 1)==0))
    {
      *mat = Eigen::MatrixXi();
      return true;
    }

    npy_intp size[2] = { _swig_numpy_array_size(in_array, 0), 1 };
    if(_swig_numpy_array_numdims(in_array)==2)
      size[1] = _swig_numpy_array_size(in_array,1);

    mat->resize(size[0], size[1]);

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

} // %fragment("Calin_Python_to_EigenIntMat"

%fragment("Calin_EigenIntMat_to_Python",
          "header",
          fragment="NumPy_Array_Requirements",
          fragment="NumPy_Backward_Compatibility",
          fragment="NumPy_Macros",
          fragment="NumPy_Utilities")
{

  static bool calin_eigen_int_mat_to_python(Eigen::MatrixXi* mat,
                                        PyObject* output)
  {
    const int typecode = NPY_INT;

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

    npy_intp size[2] = { mat->rows(), mat->cols() };
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

} // %fragment("Calin_EigenIntMat_to_Python",

// *************************** const Eigen::MatrixXi& **************************

%typemap(in, fragment="Calin_Python_to_EigenIntMat")
     const Eigen::MatrixXi& (Eigen::MatrixXi temp)
{
  // typemap(in) const Eigen::MatrixXi& -- calin_typemap_matrixxd.i
  $1 = &temp;
  if(!calin_python_to_eigen_int_mat($input, $1))SWIG_fail;
}

%typemap(out, fragment="Calin_EigenIntMat_to_Python") const Eigen::MatrixXi&
{
  // typemap(out) const Eigen::MatrixXi& -- calin_typemap_matrixxd.i
  npy_intp size[1] { $1->size() };
  $result = PyArray_EMPTY(1, size, NPY_DOUBLE, 0);
  if(!$result)SWIG_fail;
  if(!calin_eigen_int_mat_to_python($1, $result))SWIG_fail;
}

%typemap(argout) const Eigen::MatrixXi&
{
  // typemap(argout) const Eigen::MatrixXi& -- calin_typemap_matrixxd.i
  // nothing to see here
}

%typemap(typecheck, precedence=5000) const Eigen::MatrixXi&
{
  // typemap(typecheck) const Eigen::VectorXd& -- calin_typemap_vectorxd.i
  $1 = _swig_numpy_is_array($input) ? 1 : 0;
}

// ****************************** Eigen::MatrixXi& *****************************

%typemap(in, fragment="Calin_Python_to_EigenIntMat")
     Eigen::MatrixXi& (Eigen::MatrixXi temp)
{
  // typemap(in) Eigen::MatrixXi& -- calin_typemap_matrixxd.i
  $1 = &temp;
  if(!calin_python_to_eigen_int_mat($input, $1))SWIG_fail;
}

%typemap(argout, fragment="Calin_EigenIntMat_to_Python") Eigen::MatrixXi&
{
  // typemap(argout) Eigen::MatrixXi& -- calin_typemap_matrixxd.i
  if(!calin_eigen_int_mat_to_python($1, $input))SWIG_fail;
}

// ************************** Eigen::MatrixXi &OUTPUT **************************

%typemap(in, numinputs=0) Eigen::MatrixXi &OUTPUT (Eigen::MatrixXi temp)
{
  // typemap(in) Eigen::MatrixXi &OUTPUT -- calin_typemap_matrixxd.i
  $1 = &temp;
}

%typemap(argout, fragment="Calin_EigenIntMat_to_Python") Eigen::MatrixXi &OUTPUT
{
  // typemap(argout) Eigen::MatrixXi &OUTPUT -- calin_typemap_matrixxd.i
  npy_intp size[1] { $1->size() };
  PyObject* temp_array = PyArray_EMPTY(1, size, NPY_DOUBLE, 0);
  if(!temp_array)SWIG_fail;
  if(!calin_eigen_int_mat_to_python($1, temp_array))
  {
    Py_DECREF(temp_array);
    SWIG_fail;
  }
  $result = SWIG_Python_AppendOutput($result, temp_array);
}

// ************************** Eigen::MatrixXi &INOUT ***************************

%typemap(in, fragment="Calin_Python_to_EigenIntMat")
  Eigen::MatrixXi &INOUT (Eigen::MatrixXi temp)
{
  // typemap(in) const Eigen::MatrixXi &INOUT -- calin_typemap_matrixXd.i
  $1 = &temp;
  if(!calin_python_to_eigen_int_mat($input, $1))SWIG_fail;
}

%typemap(argout, fragment="Calin_EigenIntMat_to_Python") Eigen::MatrixXi &INOUT
{
  // typemap(argout) Eigen::MatrixXi &INOUT -- calin_typemap_matrixXd.i
  npy_intp size[1] { $1->size() };
  PyObject* temp_array = PyArray_EMPTY(1, size, NPY_DOUBLE, 0);
  if(!temp_array)SWIG_fail;
  if(!calin_eigen_int_mat_to_python($1, temp_array))
  {
    Py_DECREF(temp_array);
    SWIG_fail;
  }
  $result = SWIG_Python_AppendOutput($result, temp_array);
}

// ****************************** Eigen::MatrixXi ******************************

%typemap(out, fragment="Calin_EigenIntMat_to_Python") Eigen::MatrixXi
{
  // typemap(out) Eigen::MatrixXi -- calin_typemap_matrixxd.i
  npy_intp size[1] { $1.size() };
  $result = PyArray_EMPTY(1, size, NPY_DOUBLE, 0);
  if(!$result)SWIG_fail;
  if(!calin_eigen_int_mat_to_python(&$1, $result))SWIG_fail;
}

%typemap(typecheck, precedence=5000) Eigen::MatrixXi
{
  // typemap(typecheck) Eigen::MatrixXi -- calin_typemap_matrixxd.i
  $1 = _swig_numpy_is_array($input) ? 1 : 0;
}
