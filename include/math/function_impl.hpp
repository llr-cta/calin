/*

   calin/math/function_impl.hpp -- Stephen Fegan -- 2017-04-05

   Implementation of template functions for function.hpp

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

namespace calin { namespace math { namespace function {

template<typename ParameterizableBaseType>
ReducedSpaceParameterizable<ParameterizableBaseType>::~ReducedSpaceParameterizable()
{
  // nothing to see here
}

template<typename ParameterizableBaseType>
unsigned ReducedSpaceParameterizable<ParameterizableBaseType>::num_parameters()
{
  return subspace_params_.size();
}

template<typename ParameterizableBaseType> std::vector<ParameterAxis>
ReducedSpaceParameterizable<ParameterizableBaseType>::parameters()
{
  std::vector<ParameterAxis> return_parameters(subspace_params_.size());
  std::vector<ParameterAxis> axes = ParameterizableBaseType::parameters();
  for(unsigned iparam=0;iparam<subspace_params_.size();iparam++)
    return_parameters[iparam] = axes[subspace_params_[iparam]];
  return return_parameters;
}

template<typename ParameterizableBaseType> Eigen::VectorXd
ReducedSpaceParameterizable<ParameterizableBaseType>::parameter_values()
{
  return original_param_vec_to_subspace(
    ParameterizableBaseType::parametersparameter_values());
}

template<typename ParameterizableBaseType> void
ReducedSpaceParameterizable<ParameterizableBaseType>::set_parameter_values(ConstVecRef values)
{
  ParameterizableBaseType::
    set_parameter_values(subspace_param_vec_to_original(values));
}

template<typename ParameterizableBaseType> bool
ReducedSpaceParameterizable<ParameterizableBaseType>::
remove_parameter_from_subspace(unsigned iparam, double value)
{
  removed_param_values_(iparam) = value;
  auto index = std::lower_bound(subspace_params_.begin(), subspace_params_.end(), iparam);
  if(index == subspace_params_.end() or *index != iparam)return false;
  subspace_params_.erase(index);
  index = std::lower_bound(removed_params_.begin(), removed_params_.end(), iparam);
  removed_params_.insert(index, iparam);
  return true;
}

template<typename ParameterizableBaseType> bool
ReducedSpaceParameterizable<ParameterizableBaseType>::replace_parameter(unsigned iparam)
{
  auto index =
      std::lower_bound(removed_params_.begin(), removed_params_.end(), iparam);
  if(index == removed_params_.end() or *index != iparam)return false;
  removed_params_.erase(index);
  index = std::lower_bound(subspace_params_.begin(), subspace_params_.end(), iparam);
  subspace_params_.insert(index, iparam);
  return true;
}

template<typename ParameterizableBaseType> Eigen::VectorXd
ReducedSpaceParameterizable<ParameterizableBaseType>::subspace_param_vec_to_original(ConstVecRef values)
{
  assert(values.size() == subspace_params_.size());
  Eigen::VectorXd o_vec = removed_param_values_;
  for(unsigned iparam=0; iparam!=subspace_params_.size(); iparam++)
    o_vec[subspace_params_[iparam]] = values[iparam];
  return o_vec;
}

template<typename ParameterizableBaseType> Eigen::VectorXd
ReducedSpaceParameterizable<ParameterizableBaseType>::original_param_vec_to_subspace(ConstVecRef values)
{
  assert(values.size() == subspace_params_.size()+removed_params_.size());
  Eigen::VectorXd m_vec(subspace_params_.size());
  for(unsigned iparam=0; iparam!=subspace_params_.size(); iparam++)
    m_vec[iparam] = values[subspace_params_[iparam]];
  return m_vec;
}

template<typename ParameterizableBaseType> Eigen::VectorXd
ReducedSpaceParameterizable<ParameterizableBaseType>::original_param_grad_to_subspace(ConstVecRef grad)
{
  assert(grad.size() == subspace_params_.size()+removed_params_.size());
  Eigen::VectorXd subspace_grad(subspace_params_.size());
  for(unsigned iparam=0; iparam!=subspace_params_.size(); iparam++)
    subspace_grad[iparam] = grad[subspace_params_[iparam]];
  return subspace_grad;
}

template<typename ParameterizableBaseType> Eigen::MatrixXd
ReducedSpaceParameterizable<ParameterizableBaseType>::original_param_hess_to_subspace(ConstMatRef hess)
{
  assert(hess.rows() == subspace_params_.size()+removed_params_.size());
  assert(hess.cols() == subspace_params_.size()+removed_params_.size());
  Eigen::MatrixXd subspace_hess(subspace_params_.size(),this->subspace_params_.size());
  for(unsigned iparam=0; iparam!=subspace_params_.size(); iparam++)
    for(unsigned jparam=0; jparam!=subspace_params_.size(); jparam++)
      subspace_hess(iparam,jparam) = hess(subspace_params_[iparam],subspace_params_[jparam]);
  return subspace_hess;
}

template<typename ParameterizableBaseType>
ReducedSpaceParameterizableMultiAxisFunction<ParameterizableBaseType>::
~ReducedSpaceParameterizableMultiAxisFunction()
{
  // nothing to see here
}

template<typename ParameterizableBaseType>
double ReducedSpaceParameterizableMultiAxisFunction<ParameterizableBaseType>::
value_and_parameter_gradient(ConstVecRef x, VecRef gradient)
{
  Eigen::VectorXd orig_grad(this->subspace_params_.size() + this->removed_params_.size());
  double val = ParameterizableBaseType::value_and_parameter_gradient(x, orig_grad);
  gradient = this->original_param_grad_to_subspace(orig_grad);
  return val;
}

template<typename ParameterizableBaseType>
double ReducedSpaceParameterizableMultiAxisFunction<ParameterizableBaseType>::
value_parameter_gradient_and_hessian(ConstVecRef x, VecRef gradient, MatRef hessian)
{
  Eigen::VectorXd orig_grad(this->subspace_params_.size() + this->removed_params_.size());
  Eigen::MatrixXd orig_hess(this->subspace_params_.size() + this->removed_params_.size(),
    this->subspace_params_.size() + this->removed_params_.size());
  double val = ParameterizableBaseType::
    value_and_parameter_gradient_and_hessian(x, orig_grad, orig_hess);
  gradient = this->original_param_grad_to_subspace(orig_grad);
  hessian = this->original_param_hess_to_subspace(orig_hess);
  return val;
}

template<typename ParameterizableBaseType>
ReducedSpaceParameterizableSingleAxisFunction<ParameterizableBaseType>::
~ReducedSpaceParameterizableSingleAxisFunction()
{
  // nothing to see here
}

template<typename ParameterizableBaseType> double
ReducedSpaceParameterizableSingleAxisFunction<ParameterizableBaseType>::
value_and_parameter_gradient_1d(double x, VecRef gradient)
{
  Eigen::VectorXd orig_grad(this->subspace_params_.size() + this->removed_params_.size());
  double val = ParameterizableBaseType::value_and_parameter_gradient_1d(x, orig_grad);
  gradient = this->original_param_grad_to_subspace(orig_grad);
  return val;
}

template<typename ParameterizableBaseType> double
ReducedSpaceParameterizableSingleAxisFunction<ParameterizableBaseType>::
value_parameter_gradient_and_hessian_1d(double x, VecRef gradient, MatRef hessian)
{
  Eigen::VectorXd orig_grad(this->subspace_params_.size() + this->removed_params_.size());
  Eigen::MatrixXd orig_hess(this->subspace_params_.size() + this->removed_params_.size(),
    this->subspace_params_.size() + this->removed_params_.size());
  double val = ParameterizableBaseType::
    value_and_parameter_gradient_and_hessian_1d(x, orig_grad, orig_hess);
  gradient = this->original_param_grad_to_subspace(orig_grad);
  hessian = this->original_param_hess_to_subspace(orig_hess);
  return val;
}

} } } // namespace calin::math::function
