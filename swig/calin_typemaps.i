/*

   calin/calin_typemaps.i -- Stephen Fegan -- 2015-12-15

   SWIG interface file for common calin typemaps

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
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

//#define SWIGWORDSIZE64

%include "calin_stdint.i"
%include "calin_numpy.i"
%include "std_vector.i"
%include "std_string.i"
%include "typemaps.i"
%include "exception.i"

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%include "calin_typemap_vectorxd.i"
%include "calin_typemap_vectorxi.i"
%include "calin_typemap_vector3d.i"
%include "calin_typemap_matrixxd.i"
%include "calin_typemap_matrix3d.i"

// =============================================================================
//
// Typemaps for enums (why are these not part of SWIG?)
//
// =============================================================================

%typemap(in,noblock=1,numinputs=0)
  enum SWIGTYPE *CALIN_INT_OUTPUT ($basetype temp) {
  // typemap(in) enum SWIGTYPE *CALIN_INT_OUTPUT -- calin_typemaps.i
  temp = %static_cast(0,$basetype);
  $1 = &temp;
}
%typemap(argout,fragment=SWIG_From_frag(int),noblock=1)
  enum SWIGTYPE *CALIN_INT_OUTPUT {
  // typemap(argout) enum SWIGTYPE *CALIN_INT_OUTPUT -- calin_typemaps.i
  {
    PyObject* res_int = SWIG_From(int)(%static_cast(*$1, int));
    if(!res_int)SWIG_fail;
    $result = SWIG_Python_AppendOutput($result, res_int);
  }
}

// =============================================================================
//
// Typemaps for bytes arrays
//
// =============================================================================

%typemap(in) const std::string& CALIN_BYTES_IN (std::string temp, char* bytes, Py_ssize_t len) %{
  // typemap(in) const std::string& CALIN_BYTES_IN -- calin_typemaps.i
  if(PyBytes_AsStringAndSize($input,&bytes,&len) == -1)SWIG_fail;
  temp.assign(bytes, len);
  $1 = &temp;
%}
%typemap(argout) const std::string& CALIN_BYTES_IN %{
  // typemap(argout) const std::string& CALIN_BYTES_IN -- calin_typemaps.i
  // nothing to see here
%}
%typemap(freearg) const std::string& CALIN_BYTES_IN %{
  // typemap(freearg) const std::string& CALIN_BYTES_IN -- calin_typemaps.i
  // nothing to see here
%}
%typemap(typecheck) const std::string& CALIN_BYTES_IN %{
  // typemap(typecheck) const std::string& CALIN_BYTES_IN -- calin_typemaps.i
  $1 = PyBytes_Check($input) ? 1 : 0;
%}

%typemap(in, numinputs=0) std::string & CALIN_BYTES_OUT (std::string temp) %{
  // typemap(in) std::string & CALIN_BYTES_OUT -- calin_typemaps.i
  $1 = &temp;
%}
%typemap(argout) std::string & CALIN_BYTES_OUT %{
  // typemap(argout) std::string & CALIN_BYTES_OUT -- calin_typemaps.i
  {
    PyObject* temp_bytes = PyBytes_FromStringAndSize(&$1->front(), $1->size());
    if(!temp_bytes)SWIG_fail;
    $result = SWIG_Python_AppendOutput($result, temp_bytes);
  }
%}
%typemap(freearg) const std::string & CALIN_BYTES_OUT %{
  // typemap(freearg) std::string& CALIN_BYTES_OUT -- calin_typemaps.i
  // nothing to see here
%}

// =============================================================================
//
// Typemaps for using Eigen::Ref - these avoid copies if data types match
//
// =============================================================================

%typemap(in,
         fragment="NumPy_Fragments")
         const Eigen::Ref<const Eigen::VectorXd>&
         (PyArrayObject* array=NULL, int is_new_object=0,
          Eigen::Map<const Eigen::VectorXd>* eigenmap=0)
{
  npy_intp size[1] = { -1 };
  array = obj_to_array_contiguous_allow_conversion($input,
                                                   NPY_DOUBLE,
                                                   &is_new_object);
  if (!array || !require_dimensions(array, 1) ||
      !require_size(array, size, 1)) SWIG_fail;
  eigenmap =
      new Eigen::Map<const Eigen::VectorXd>((const double*)array_data(array),
                                            _swig_numpy_array_size(array,0));
  $1 = new $*1_ltype(*eigenmap);
}

%typemap(freearg) const Eigen::Ref<const Eigen::VectorXd>&
{
  delete arg$argnum;
  delete eigenmap$argnum;
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}
