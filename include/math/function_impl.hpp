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

#include <numeric>

namespace calin { namespace math { namespace function {

template<typename ParameterizableBaseType>
FrozenParameterizable<ParameterizableBaseType>::
FrozenParameterizable(ParameterizableBaseType* p, bool adopt_p):
  ParameterizableBaseType(), p_(p), adopt_p_(adopt_p),
  free_params_(p->num_parameters()), frozen_params_(),
  values_frozen_(p->num_parameters())
{
  std::iota(free_params_.begin(), free_params_.end(), 0);
}

template<typename ParameterizableBaseType>
FrozenParameterizable<ParameterizableBaseType>::~FrozenParameterizable()
{
  if(adopt_p_)delete p_;
}

template<typename ParameterizableBaseType>
unsigned FrozenParameterizable<ParameterizableBaseType>::num_parameters()
{
  return free_params_.size();
}

template<typename ParameterizableBaseType> std::vector<ParameterAxis>
FrozenParameterizable<ParameterizableBaseType>::parameters()
{
  std::vector<ParameterAxis> return_parameters(free_params_.size());
  std::vector<ParameterAxis> axes = p_->parameters();
  for(unsigned iparam=0;iparam<free_params_.size();iparam++)
    return_parameters[iparam] = axes[free_params_[iparam]];
  return return_parameters;

}

template<typename ParameterizableBaseType> Eigen::VectorXd
FrozenParameterizable<ParameterizableBaseType>::parameter_values()
{
  return original_param_vec_to_modified(p_->parameter_values());
}

template<typename ParameterizableBaseType> void
FrozenParameterizable<ParameterizableBaseType>::set_parameter_values(ConstVecRef values)
{
  p_->set_parameter_values(modified_param_vec_to_original(values));
}

template<typename ParameterizableBaseType> bool
FrozenParameterizable<ParameterizableBaseType>::can_calculate_parameter_gradient()
{
  return p_->can_calculate_parameter_gradient();
}

template<typename ParameterizableBaseType> bool
FrozenParameterizable<ParameterizableBaseType>::can_calculate_parameter_hessian()
{
  return p_->can_calculate_parameter_hessian();
}

template<typename ParameterizableBaseType> bool
FrozenParameterizable<ParameterizableBaseType>::freeze_parameter(unsigned iparam, double value)
{
  values_frozen_(iparam) = value;
  auto index = std::lower_bound(free_params_.begin(), free_params_.end(), iparam);
  if(index == free_params_.end() or *index != iparam)return false;
  free_params_.erase(index);
  index = std::lower_bound(frozen_params_.begin(), frozen_params_.end(), iparam);
  frozen_params_.insert(index, iparam);
  return true;
}

template<typename ParameterizableBaseType> bool
FrozenParameterizable<ParameterizableBaseType>::thaw_parameter(unsigned iparam)
{
  auto index =
      std::lower_bound(frozen_params_.begin(), frozen_params_.end(), iparam);
  if(index == frozen_params_.end() or *index != iparam)return false;
  frozen_params_.erase(index);
  index = std::lower_bound(free_params_.begin(), free_params_.end(), iparam);
  free_params_.insert(index, iparam);
  return true;
}

template<typename ParameterizableBaseType> Eigen::VectorXd
FrozenParameterizable<ParameterizableBaseType>::modified_param_vec_to_original(ConstVecRef values)
{
  assert(values.size() == free_params_.size());
  Eigen::VectorXd o_vec = values_frozen_;
  for(unsigned iparam=0; iparam!=free_params_.size(); iparam++)
    o_vec[free_params_[iparam]] = values[iparam];
  return o_vec;
}

template<typename ParameterizableBaseType> Eigen::VectorXd
FrozenParameterizable<ParameterizableBaseType>::original_param_vec_to_modified(ConstVecRef values)
{
  assert(values.size() == free_params_.size()+frozen_params_.size());
  Eigen::VectorXd m_vec(free_params_.size());
  for(unsigned iparam=0; iparam!=free_params_.size(); iparam++)
    m_vec[iparam] = values[free_params_[iparam]];
  return m_vec;
}

} } } // namespace calin::math::function
