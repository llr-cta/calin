#!/usr/bin/env python3

# calin/scripts/compute_diagnostics.py -- Stephen Fegan - 2018-11-20
#
# Serve a CTA R1 ZFits via a ZMQ endpoint. This is not fast as it must pass
# through the JSON interface to the underlying protobufs. It could easily be
# improved.
#
# Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
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

import calin.io.zmq_inproc_push_pull
import calin.iact_data.raw_actl_l0_event_data_source
import calin.iact_data.raw_actl_r1_event_data_source
import sys
import json
import base64

endpoint = 'tcp://0.0.0.0:20000'

if len(sys.argv) == 1:
    print("Usage:",sys.argv[0],"filename [endpoint="+endpoint+"]")
    exit(1)

filename = sys.argv[1]
src = calin.iact_data.raw_actl_r1_event_data_source.ZFITSACTL_R1_CameraEventDataSource(filename)

if len(sys.argv) >= 3:
    endpoint = sys.argv[2]

ctx = calin.io.zmq_inproc_push_pull.new_zmq_ctx()
pusher = calin.io.zmq_inproc_push_pull.ZMQPusher(ctx, endpoint)

cta_msg = calin.iact_data.raw_actl_l0_event_data_source.CTAMessage()
msg = dict()
msg['sourceName'] = filename
msg['messageCount'] = 0

rh = src.get_run_header()
if rh is not None:
    msg['payloadType'] = ['R1_CAMERA_CONFIG']
    msg['payloadData'] = [base64.b64encode(rh.SerializeAsString()).decode("utf-8")];
    cta_msg.ParseFromJSON(json.dumps(msg))
    pusher.push(cta_msg.SerializeAsString())

msg['payloadType'] = ['R1_CAMERA_EVENT']
nevent = 0
seq_id, event = src.get_next()
while event:
    msg['messageCount'] += 1
    msg['payloadData'] = [base64.b64encode(event.SerializeAsString()).decode("utf-8")];
    cta_msg.ParseFromJSON(json.dumps(msg))
    pusher.push(cta_msg.SerializeAsString())
    nevent += 1
    if(nevent % 10000 == 0):
        print("Pushed",nevent,"events")
    seq_id, event = src.get_next()
print("Pushed",nevent,"events.. finished")

msg['payloadType'] = ['END_OF_STREAM']
msg['payloadData'] = ''
cta_msg.ParseFromJSON(json.dumps(msg))
pusher.push(cta_msg.SerializeAsString())
