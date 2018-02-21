/*

   calin/math/simd_fft.hpp -- Stephen Fegan -- 2018-02-21

   SIMD FFT functions using codelets from FFTW/genfft

   Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
   LLR, Ecole Polytechnique, CNRS/IN2P3

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

#if defined(__AVX__)
#include <immintrin.h>
#endif // defined(__AVX__)

#if defined(__AVX__)

namespace calin { namespace math { namespace simd_fft { namespace m256 {
using float_type = __m256;
#include "genfft_codelets/declare_all_dft_codelets.hpp_incl"
} } } } // namespace calin::math::simd_fft::m256

namespace calin { namespace math { namespace simd_fft { namespace m256d {
using float_type = __m256d;
#include "genfft_codelets/declare_all_dft_codelets.hpp_incl"
} } } } // namespace calin::math::simd_fft::m256d

#endif // defined(__AVX__)
