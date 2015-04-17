/* 

   calin/math/function.i -- Stephen Fegan -- 2015-04-15

*/
//-*-mode:swig;-*-

%module (package="calin.math") function

%{
#include "package_wide_definitions.hpp"
#include "math/function.hpp"
  %}

%include "numpy.i"

%include "package_wide_definitions.hpp"

%typemap(in,
         fragment="NumPy_Fragments")
   calin::ConstVecRef
(PyArrayObject* array=NULL, int is_new_object=0)
{
  npy_intp size[1] = { -1 };
  array = obj_to_array_contiguous_allow_conversion($input,
                                                   NPY_DOUBLE,
                                                   &is_new_object);
  if (!array || !require_dimensions(array, 1) ||
      !require_size(array, size, 1)) SWIG_fail;
  $1 = new Eigen::VectorXd();
  *$1 = Eigen::Map<const Eigen::VectorXd>((const double*)array_data(array),
                                          array_size(array,0));
}
%typemap(freearg) calin::ConstVecRef
{
  delete arg$argnum;
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}


%include "math/function.hpp"
