/*

   calin/pattern/delegation.hpp -- Stephen Fegan -- 2017-04-25

   Base class for objects that delegate tasks

   Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
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

#pragma once

namespace calin { namespace pattern { namespace delegation {

template<typename T> class Delegator
{
public:
  Delegator(T* delegate, bool adopt_delegate = false):
    delegate_(delegate), adopt_delegate_(adopt_delegate) { }
  virtual ~Delegator() { if(adopt_delegate_)delete delegate_; }
  T* delegate() { return delegate_; }
  void set_delegate(T* delegate, bool adopt_delegate = false) {
    if(adopt_delegate_)delete delegate_;
    delegate_ = delegate;
    adopt_delegate_ = adopt_delegate;
  }
protected:
  T* delegate_ = nullptr;
  bool adopt_delegate_ = false;
};

} } } // namespace calin::pattern::delegation
