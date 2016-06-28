#!/usr/bin/env python

# calin/scripts/compute_diagnostics.py -- Stephen Fegan - 2016-06-10
#
# Compute diagnostics from ZFits file, saving them to SQLITE3
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
import numpy

import calin.iact_data.telescope_data_source
import calin.iact_data.event_dispatcher
import calin.diagnostics.waveform
import calin.diagnostics.functional
import calin.diagnostics.value_capture
import calin.diagnostics.module
import calin.diagnostics.event_number
import calin.diagnostics.delta_t
import calin.io.sql_transceiver
import calin.io.log

py_log = calin.io.log.PythonLogger()
py_log.this.disown()
calin.io.log.default_logger().add_logger(py_log,True)
proto_log = calin.io.log.ProtobufLogger()
proto_log.this.disown()
calin.io.log.default_logger().add_logger(proto_log,True)

if(len(sys.argv) == 1):
    raise Exception('No filename supplied')
zfits_file = sys.argv[1]

sql_file = 'diagnostics.sqlite'
if(len(sys.argv) > 2):
    sql_file = sys.argv[2]

bkg_window_start = 0
if(len(sys.argv) > 3):
    bkg_window_start = int(sys.argv[3])

sig_window_start = 44
if(len(sys.argv) > 4):
    sig_window_start = int(sys.argv[4])

# Open the data source
cfg = calin.iact_data.telescope_data_source.\
    NectarCamZFITSDataSource.default_config()
#cfg.set_max_file_fragments(8)
src = calin.iact_data.telescope_data_source.\
    NectarCamZFITSDataSource(zfits_file, cfg)
src.set_next_index(1)

# Get the run info
run_info = src.get_run_configuration()

# Create the dispatcher
dispatcher = calin.iact_data.event_dispatcher.TelescopeEventDispatcher()

# Background window functional
bkg_window_sum_cfg = calin.iact_data.functional_event_visitor.\
    FixedWindowSumFunctionalTelescopeEventVisitor.default_config()
bkg_window_sum_cfg.set_integration_0(bkg_window_start)
bkg_window_sum_cfg.set_integration_n(16)
bkg_window_sum_visitor = calin.iact_data.functional_event_visitor.\
    FixedWindowSumFunctionalTelescopeEventVisitor(bkg_window_sum_cfg)
dispatcher.add_visitor(bkg_window_sum_visitor, \
    calin.iact_data.event_dispatcher.EXECUTE_SEQUENTIAL_AND_PARALLEL)

bkg_capture_adapter = [None] * run_info.configured_channel_id_size();
bkg_capture = [None] * run_info.configured_channel_id_size();
for ichan in range(run_info.configured_channel_id_size()):
    # Background capture adapter - select channel
    bkg_capture_adapter[ichan] = calin.diagnostics.functional.\
        SingleFunctionalValueSupplierVisitor(bkg_window_sum_visitor,ichan)
    dispatcher.add_visitor(bkg_capture_adapter[ichan])

    # Background capture
    bkg_capture[ichan] = calin.diagnostics.value_capture.\
        IntSequentialValueCaptureVisitor(bkg_capture_adapter[ichan],0x7FFFFFFF)
    dispatcher.add_visitor(bkg_capture[ichan])

# Signal window functional
sig_window_sum_cfg = calin.iact_data.functional_event_visitor.\
    FixedWindowSumFunctionalTelescopeEventVisitor.default_config()
#sig_window_sum_cfg.set_integration_0(30)
sig_window_sum_cfg.set_integration_0(sig_window_start)
sig_window_sum_cfg.set_integration_n(16)
sig_window_sum_visitor = calin.iact_data.functional_event_visitor.\
    FixedWindowSumFunctionalTelescopeEventVisitor(sig_window_sum_cfg)
dispatcher.add_visitor(sig_window_sum_visitor, \
    calin.iact_data.event_dispatcher.EXECUTE_SEQUENTIAL_AND_PARALLEL)

# Signal minus background functional
sig_bkg_diff_visitor = calin.iact_data.functional_event_visitor.\
    DifferencingFunctionalTelescopeEventVisitor(sig_window_sum_visitor, \
        bkg_window_sum_visitor)
dispatcher.add_visitor(sig_bkg_diff_visitor, \
    calin.iact_data.event_dispatcher.EXECUTE_SEQUENTIAL_AND_PARALLEL)

# Signal minus background stats
sig_bkg_stats_visitor = calin.diagnostics.functional.\
    FunctionalIntStatsVisitor(sig_bkg_diff_visitor)
dispatcher.add_visitor(sig_bkg_stats_visitor)

# Signal minus background capture adapter
sig_bkg_capture_adapter = calin.diagnostics.functional.\
    SingleFunctionalValueSupplierVisitor(sig_bkg_diff_visitor,0)
dispatcher.add_visitor(sig_bkg_capture_adapter)

# Signal minus background capture
sig_bkg_capture = calin.diagnostics.value_capture.\
    IntSequentialValueCaptureVisitor(sig_bkg_capture_adapter,0x7FFFFFFF)
dispatcher.add_visitor(sig_bkg_capture)

# Waveform stats
waveform_visitor = calin.diagnostics.waveform.WaveformStatsVisitor()
dispatcher.add_visitor(waveform_visitor)

# Waveform PSD stats
psd_visitor = calin.diagnostics.waveform.WaveformPSDVisitor()
dispatcher.add_visitor(psd_visitor)

# Background window stats
bkg_window_stats_visitor = calin.diagnostics.functional.\
    FunctionalIntStatsVisitor(bkg_window_sum_visitor)
dispatcher.add_visitor(bkg_window_stats_visitor)

# Glitch detection visitor
glitch_visitor = \
    calin.diagnostics.event_number.CountersEventNumberGlitchDetector()
dispatcher.add_visitor(glitch_visitor)

# Bunch event number detection visitor
bunch_event_glitch_visitor = \
    calin.diagnostics.event_number.CountersEventNumberGlitchDetector(2)
dispatcher.add_visitor(bunch_event_glitch_visitor)

# Module present visitor
mod_present_visitor = \
    calin.diagnostics.module.ModulePresentVisitor()
dispatcher.add_visitor(mod_present_visitor)

# Delta-T calculation visitor
delta_t_calc_visitor = \
    calin.diagnostics.delta_t.OneModuleCounterDeltaTVisitor(0,3,8)
dispatcher.add_visitor(delta_t_calc_visitor)

# Delta-T capture
delta_t_capture = calin.diagnostics.value_capture.\
    DoubleSequentialValueCaptureVisitor(delta_t_calc_visitor,numpy.nan)
dispatcher.add_visitor(delta_t_capture)

# T0 rise time functional
t0_calc = calin.iact_data.functional_event_visitor.\
    NoPedestalTimingFunctionalTelescopeEventVisitor()
dispatcher.add_visitor(t0_calc)

# T0 rise time stats
t0_stats_cfg = calin.diagnostics.functional.\
    FunctionalDoubleStatsVisitor.default_config()
t0_stats_cfg.hist_config().set_dxval(0.1)
t0_stats_cfg.hist_config().set_xval_units('samples')
t0_stats = calin.diagnostics.functional.\
    FunctionalDoubleStatsVisitor(t0_calc, t0_stats_cfg)
dispatcher.add_visitor(t0_stats)

# Run all the visitors
dispatcher.process_run(src,100000,3)

# Get the results
psd = psd_visitor.results()
wfs = waveform_visitor.results()
bkg = bkg_window_stats_visitor.results()
sig = sig_bkg_stats_visitor.results()
glitch = glitch_visitor.glitch_data()
bunch_event_glitch = bunch_event_glitch_visitor.glitch_data()
mod_present = mod_present_visitor.module_data()
sig_values = sig_bkg_capture.results()
bkg_values = bkg_capture[0].results()
delta_t_values = delta_t_capture.results()
t0 = t0_stats.results()

# Write the results
sql = calin.io.sql_transceiver.SQLite3Transceiver(sql_file,
    calin.io.sql_transceiver.SQLite3Transceiver.TRUNCATE_RW)
sql.create_tables_and_insert("log", proto_log.log_messages())
sql.create_tables_and_insert("run_config", run_info)
sql.create_tables_and_insert("psd", psd)
sql.create_tables_and_insert("wfs", wfs)
sql.create_tables_and_insert("bkg", bkg)
sql.create_tables_and_insert("sig", sig)
sql.create_tables_and_insert("glitch_event", glitch)
sql.create_tables_and_insert("glitch_bunch_event", bunch_event_glitch)
sql.create_tables_and_insert("mod_present", mod_present)
sql.create_tables_and_insert("sig_values", sig_values)
sql.create_tables_and_insert("bkg_values", bkg_values)
sql.create_tables_and_insert("delta_t_values", delta_t_values)
sql.create_tables_and_insert("t0", t0)

for ichan in range(1,len(bkg_capture)):
    sql.insert("bkg_values", bkg_capture[ichan].results())
