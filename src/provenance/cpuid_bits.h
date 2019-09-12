/*

   calin/provenance/cpuid_bits.h -- Stephen Fegan -- 2017-10-27

   CPUID bits that may be missing from certain compiler versions

   Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   Based on original, copyright 2006, Stephen Fegan, see notice below

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

#ifndef bit_AVX2
#define bit_AVX2        (1 << 5)
#endif

#ifndef bit_BMI
#define bit_BMI         (1 << 3)
#endif

#ifndef bit_BMI2
#define bit_BMI2        (1 << 8)
#endif

#ifndef bit_ADX
#define bit_ADX         (1 << 19)
#endif

#ifndef bit_AVX512F
#define bit_AVX512F     (1 << 16)
#endif

#ifndef bit_AVX512DQ
#define bit_AVX512DQ    (1 << 17)
#endif

#ifndef bit_AVX512IFMA
#define bit_AVX512IFMA  (1 << 21)
#endif

#ifndef bit_AVX512PF
#define bit_AVX512PF    (1 << 26)
#endif

#ifndef bit_AVX512ER
#define bit_AVX512ER    (1 << 27)
#endif

#ifndef bit_AVX512CD
#define bit_AVX512CD    (1 << 28)
#endif

#ifndef bit_AVX512BW
#define bit_AVX512BW    (1 << 30)
#endif

#ifndef bit_AVX512VL
#define bit_AVX512VL    (1u << 31)
#endif

#ifndef bit_AVX512VBMI
#define bit_AVX512VBMI  (1 << 1)
#endif

#ifndef bit_AVX512VBMI2
#define bit_AVX512VBMI2 (1 << 6)
#endif

#ifndef bit_AVX512VNNI
#define bit_AVX512VNNI  (1 << 11)
#endif

#ifndef bit_AVX512BITALG
#define bit_AVX512BITALG (1 << 12)
#endif

#ifndef bit_AVX512VPOPCNTDQ
#define bit_AVX512VPOPCNTDQ (1 << 14)
#endif

#ifndef bit_AVX512_4VNNIW
#define bit_AVX512_4VNNIW   (1 << 2)
#endif

#ifndef bit_AVX512_4FMAPS
#define bit_AVX512_4FMAPS   (1 << 3)
#endif

#ifndef bit_FMA4
#define bit_FMA4        (1 << 16)
#endif
