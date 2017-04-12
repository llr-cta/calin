# calin/python/iact_data/ipno.py -- Stephen Fegan -- 2017-04-07
#
# Functions for manipulating data from IPNO test bench
#
# Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
# LLR, Ecole polytechnique, CNRS/IN2P3, Universite Paris-Saclay
#
# This file is part of "calin"
#
# "calin" is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License version 2 or later, as published by
# the Free Software Foundation.
#
# "calin" is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

import struct
import numpy
import calin.math.histogram

def read_ipno_comment(f):
    data = open(f,'rb').read()
    i64 = numpy.asarray(list(map(lambda x: struct.unpack_from('>I', x)[0], \
                [data[i:i+4] for i in range(0, len(data), 4)])))
    return ''.join(map(chr, i64[2:i64[1]-1]))

def read_ipno_adc_data(f, dtype = 'float'):
    data = open(f,'rb').read()
    i64 = numpy.asarray(list(map(lambda x: struct.unpack_from('>I', x)[0], \
                [data[i:i+4] for i in range(0, len(data), 4)])))
    index = numpy.where(i64 == 0xffffffff)[0][1:]
    return numpy.asarray(i64[index-2], dtype = dtype)

def make_ipno_adc_hist(f, value_min=-numpy.inf, value_max=numpy.inf):
    vals = read_ipno_adc_data(f)
    hist = calin.math.histogram.SimpleHist(1)
    hist.insert_vec(vals[numpy.bitwise_and(vals>=value_min, vals<=value_max)]);
    return hist
