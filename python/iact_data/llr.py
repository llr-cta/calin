# calin/python/iact_data/llr.py -- Stephen Fegan -- 2017-05-18
#
# Functions for manipulating data from LLR test bench
#
# Copyright 2017, Stephen Fegan <sfegan@llr.in2p3.fr>
# LLR, Ecole Polytechnique, CNRS/IN2P3
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

def make_lecroy_adc_hist(f, scale=1.0):
    data = numpy.loadtxt(f, delimiter=',', skiprows=5)
    dx = numpy.median(data[1:,0] - data[0:-1,0])/scale
    hist = calin.math.histogram.SimpleHist(dx)
    hist.insert_two_vec(data[:,0]/scale, data[:,1])
    return hist
