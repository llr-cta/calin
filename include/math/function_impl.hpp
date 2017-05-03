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

// #include<io/log.hpp>

namespace calin { namespace math { namespace function {

template<typename T>
BasicReducedSpaceParameterizable<T>::~BasicReducedSpaceParameterizable()
{
  // nothing to see here
}

template<typename T>
unsigned BasicReducedSpaceParameterizable<T>::num_parameters()
{
  return subspace_params_.size();
}

template<typename T> std::vector<ParameterAxis>
BasicReducedSpaceParameterizable<T>::parameters()
{
  std::vector<ParameterAxis> return_parameters(subspace_params_.size());
  std::vector<ParameterAxis> axes = this->delegate_->parameters();
  for(unsigned iparam=0;iparam<subspace_params_.size();iparam++)
    return_parameters[iparam] = axes[subspace_params_[iparam]];
  return return_parameters;
}

template<typename T> Eigen::VectorXd
BasicReducedSpaceParameterizable<T>::parameter_values()
{
  return original_param_vec_to_subspace(this->delegate_->parameter_values());
}

template<typename T> void
BasicReducedSpaceParameterizable<T>::set_parameter_values(ConstVecRef values)
{
  if(values.size() != subspace_params_.size())
    throw std::invalid_argument(
      std::string("Number of parameter values does not match subspace size: "
      + std::to_string(values.size()) + " != "
      + std::to_string(subspace_params_.size())));
  // calin::io::log::LOG(calin::io::log::INFO) << "A: " << values.transpose();
  // calin::io::log::LOG(calin::io::log::INFO) << "B: " << removed_param_values_.transpose();
  // calin::io::log::LOG(calin::io::log::INFO) << "C: " << subspace_param_vec_to_original(values).transpose();
  this->delegate_->set_parameter_values(subspace_param_vec_to_original(values));
}

template<typename T> bool BasicReducedSpaceParameterizable<T>::
remove_parameter_from_subspace(unsigned unmodified_index)
{
  auto index = std::lower_bound(subspace_params_.begin(), subspace_params_.end(), unmodified_index);
  if(index == subspace_params_.end() or *index != unmodified_index)return false;
  removed_param_values_(unmodified_index) = this->delegate_->parameter_values()[unmodified_index];
  subspace_params_.erase(index);
  index = std::lower_bound(removed_params_.begin(), removed_params_.end(), unmodified_index);
  removed_params_.insert(index, unmodified_index);
  return true;
}

template<typename T> bool BasicReducedSpaceParameterizable<T>::
remove_parameter_from_subspace(const std::string& name)
{
  int index = this->delegate_->parameter_name_to_index(name);
  if(index<0)return false;
  return remove_parameter_from_subspace(unsigned(index));
}

template<typename T> bool BasicReducedSpaceParameterizable<T>::
remove_parameter_from_subspace(unsigned unmodified_index, double value)
{
  removed_param_values_(unmodified_index) = value;
  auto index = std::lower_bound(subspace_params_.begin(), subspace_params_.end(), unmodified_index);
  if(index == subspace_params_.end() or *index != unmodified_index)return false;
  subspace_params_.erase(index);
  index = std::lower_bound(removed_params_.begin(), removed_params_.end(), unmodified_index);
  removed_params_.insert(index, unmodified_index);
  return true;
}

template<typename T> bool BasicReducedSpaceParameterizable<T>::
remove_parameter_from_subspace(const std::string& name, double value)
{
  int index = this->delegate_->parameter_name_to_index(name);
  if(index<0)return false;
  return remove_parameter_from_subspace(unsigned(index), value);
}

template<typename T> void BasicReducedSpaceParameterizable<T>::
remove_all_remaining_parameters_from_subspace()
{
  while(not subspace_params_.empty())
    remove_parameter_from_subspace(subspace_params_.front());
}

template<typename T> bool
BasicReducedSpaceParameterizable<T>::replace_parameter_in_subspace(unsigned unmodified_index)
{
  auto index =
      std::lower_bound(removed_params_.begin(), removed_params_.end(), unmodified_index);
  if(index == removed_params_.end() or *index != unmodified_index)return false;
  removed_params_.erase(index);
  index = std::lower_bound(subspace_params_.begin(), subspace_params_.end(), unmodified_index);
  subspace_params_.insert(index, unmodified_index);
  return true;
}

template<typename T> bool BasicReducedSpaceParameterizable<T>::
replace_parameter_in_subspace(const std::string& name)
{
  int index = this->delegate_->parameter_name_to_index(name);
  if(index<0)return false;
  return replace_parameter_in_subspace(unsigned(index));
}

template<typename T> void BasicReducedSpaceParameterizable<T>::
replace_all_parameters_in_subspace()
{
  while(not removed_params_.empty())
    replace_parameter_in_subspace(removed_params_.front());
}

template<typename T> Eigen::VectorXd
BasicReducedSpaceParameterizable<T>::subspace_param_vec_to_original(ConstVecRef values)
{
  assert(values.size() == subspace_params_.size());
  Eigen::VectorXd o_vec = removed_param_values_;
  for(unsigned iparam=0; iparam!=subspace_params_.size(); iparam++)
    o_vec[subspace_params_[iparam]] = values[iparam];
  return o_vec;
}

template<typename T> Eigen::VectorXd
BasicReducedSpaceParameterizable<T>::original_param_vec_to_subspace(ConstVecRef values)
{
  assert(values.size() == subspace_params_.size()+removed_params_.size());
  Eigen::VectorXd m_vec(subspace_params_.size());
  for(unsigned iparam=0; iparam!=subspace_params_.size(); iparam++)
    m_vec[iparam] = values[subspace_params_[iparam]];
  return m_vec;
}

template<typename T> Eigen::VectorXd
BasicReducedSpaceParameterizable<T>::original_param_grad_to_subspace(ConstVecRef grad)
{
  assert(grad.size() == subspace_params_.size()+removed_params_.size());
  Eigen::VectorXd subspace_grad(subspace_params_.size());
  for(unsigned iparam=0; iparam!=subspace_params_.size(); iparam++)
    subspace_grad[iparam] = grad[subspace_params_[iparam]];
  return subspace_grad;
}

template<typename T> Eigen::MatrixXd
BasicReducedSpaceParameterizable<T>::original_param_hess_to_subspace(ConstMatRef hess)
{
  assert(hess.rows() == subspace_params_.size()+removed_params_.size());
  assert(hess.cols() == subspace_params_.size()+removed_params_.size());
  Eigen::MatrixXd subspace_hess(subspace_params_.size(),this->subspace_params_.size());
  for(unsigned iparam=0; iparam!=subspace_params_.size(); iparam++)
    for(unsigned jparam=0; jparam!=subspace_params_.size(); jparam++)
      subspace_hess(iparam,jparam) = hess(subspace_params_[iparam],subspace_params_[jparam]);
  return subspace_hess;
}

template<typename T>
BasicReducedSpaceParameterizableSingleAxisFunction<T>::
~BasicReducedSpaceParameterizableSingleAxisFunction()
{
  // nothing to see here
}

template<typename T>
double BasicReducedSpaceParameterizableSingleAxisFunction<T>::
value_and_parameter_gradient_1d(double x, VecRef gradient)
{
  Eigen::VectorXd orig_grad(this->subspace_params_.size() + this->removed_params_.size());
  double val = this->delegate_->value_and_parameter_gradient_1d(x, orig_grad);
  gradient = this->original_param_grad_to_subspace(orig_grad);
  return val;
}

template<typename T>
double BasicReducedSpaceParameterizableSingleAxisFunction<T>::
value_parameter_gradient_and_hessian_1d(double x, VecRef gradient, MatRef hessian)
{
  unsigned norig = this->subspace_params_.size() + this->removed_params_.size();
  Eigen::VectorXd orig_grad(norig);
  Eigen::MatrixXd orig_hess(norig, norig);
  double val = this->delegate_->value_parameter_gradient_and_hessian_1d(x, orig_grad, orig_hess);
  gradient = this->original_param_grad_to_subspace(orig_grad);
  hessian = this->original_param_hess_to_subspace(orig_hess);
  return val;
}

} } } // namespace calin::math::function
