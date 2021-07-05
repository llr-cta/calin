#!/usr/bin/env python3

# calin/scripts/cta_diagnostics.py -- Stephen Fegan - 2018-11-28
#
# Compute diagnostics from ZFits file or ZMQ streams, saving them to SQLITE3
#
# Copyright 2018, Stephen Fegan <sfegan@llr.in2p3.fr>
# Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris
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

import calin.iact_data.raw_actl_r1_event_data_source
import calin.ix.io.zmq_data_source
import calin.diagnostics.run_info
import calin.iact_data.event_dispatcher
import calin.io.sql_transceiver
import calin.util.log
import calin.provenance.anthology
import calin.util.options_processor
import calin.ix.scripts.cta_diagnostics

py_log = calin.util.log.PythonLogger()
py_log.this.disown()
calin.util.log.default_logger().add_logger(calin.util.log.default_protobuf_logger(),False)
calin.util.log.default_logger().add_logger(py_log,True)

cfg = calin.iact_data.event_dispatcher.ParallelEventDispatcher.default_config()

opt = calin.ix.scripts.cta_diagnostics.CommandLineOptions()
opt.set_run_number(0)
opt.set_o('diagnostics.sqlite')
opt.set_db_stage1_table_name('diagnostics_results')
opt.set_log_frequency(cfg.log_frequency())
opt.set_nthread(1)
opt.mutable_decoder().CopyFrom(cfg.decoder())
opt.mutable_zfits().CopyFrom(cfg.zfits())
opt.mutable_zmq().CopyFrom(cfg.zmq())
opt.mutable_run_info().CopyFrom(calin.diagnostics.run_info.RunInfoDiagnosticsParallelEventVisitor.default_config())

opt_proc = calin.util.options_processor.OptionsProcessor(opt, True);
opt_proc.process_arguments(sys.argv)

if(opt_proc.help_requested()):
    print('Usage:',opt_proc.program_name(),'[options] zmq_endpoint [zmq_endpoint...]')
    print('or:   ',opt_proc.program_name(),'[options] zfits_file\n')
    print('Options:\n')
    print(opt_proc.usage())
    exit(0)

if(len(opt_proc.arguments()) < 1):
    print('No endpoint supplied! Use "-help" option to get usage information.')
    exit(1)

if(len(opt_proc.unknown_options()) != 0):
    print('Unknown options given. Use "-help" option to get usage information.\n')
    for o in opt_proc.unknown_options():
        print("  \"%s\""%o)
    exit(1)

if(len(opt_proc.problem_options()) != 0):
    print('Problems with option values (unexpected, missing, incorrect type, etc.).')
    print('Use "-help" option to get usage information.\n')
    for o in opt_proc.problem_options():
        print("  \"%s\""%o)
    exit(1)

endpoints          = opt_proc.arguments()
sql_file           = opt.o();

cfg.set_run_number(opt.run_number())
cfg.set_log_frequency(opt.log_frequency())
cfg.set_nthread(opt.nthread())
cfg.mutable_decoder().CopyFrom(opt.decoder())
cfg.mutable_zfits().CopyFrom(opt.zfits())
cfg.mutable_zmq().CopyFrom(opt.zmq())

# Create the dispatcher
dispatcher = calin.iact_data.event_dispatcher.ParallelEventDispatcher()

# Create the run info visitor
ri = calin.diagnostics.run_info.RunInfoDiagnosticsParallelEventVisitor(opt.run_info())
dispatcher.add_visitor(ri)

# Run all the visitors
if(len(endpoints)>1 or endpoints[0].startswith('tcp://')
        or endpoints[0].startswith('ipc://') or endpoints[0].startswith('pgm://')
        or endpoints[0].startswith('pgme://')):
    dispatcher.process_cta_zmq_run(endpoints, cfg)
else:
    dispatcher.process_cta_zfits_run(endpoints[0], cfg)

# Open SQL file
sql = calin.io.sql_transceiver.SQLite3Transceiver(sql_file,
    calin.io.sql_transceiver.SQLite3Transceiver.TRUNCATE_RW)

# Get the results
results = calin.ix.scripts.cta_diagnostics.Results()
results.mutable_command_line_arguments().CopyFrom(opt_proc.command_line_arguments())
results.mutable_command_line_options().CopyFrom(opt)
calin.provenance.anthology.get_current_anthology(results.mutable_provenance())
results.mutable_run_config().CopyFrom(ri.run_config())
results.mutable_run_info().CopyFrom(ri.run_info())

# Write the results
sql.create_tables(opt.db_stage1_table_name(), results.descriptor())
sql.insert(opt.db_stage1_table_name(), results)
