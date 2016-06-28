#!/usr/bin/env python

# calin/scripts/make_proto_deps.py -- Stephen Fegan - 2016-06-28
#
# Make CMake depenency list by parsing proto file
#
# Copyright 2016, Stephen Fegan <sfegan@llr.in2p3.fr>
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

import sys

print(sys.argv[1])
f = open(sys.argv[1])
#if(f == None):
