#!/usr/bin/env python3

# calin/scripts/cta_stage1.py -- Stephen Fegan - 2020-07-21
#
# Compute stage1 data from ZFits file or ZMQ streams, saving them to SQLITE3
#
# Copyright 2020, Stephen Fegan <sfegan@llr.in2p3.fr>
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
import traceback
import numpy

import calin.iact_data.raw_actl_r1_event_data_source
import calin.ix.io.zmq_data_source
import calin.iact_data.event_dispatcher
import calin.io.sql_serializer
import calin.util.log
import calin.util.options_processor
import calin.diagnostics.stage1
import calin.ix.scripts.cta_stage1
import calin.provenance.chronicle

py_log = calin.util.log.PythonLogger()
py_log.this.disown()
calin.util.log.default_logger().add_logger(calin.util.log.default_protobuf_logger(),False)
calin.util.log.default_logger().add_logger(py_log,True)

cfg = calin.iact_data.event_dispatcher.ParallelEventDispatcher.default_config()

opt = calin.ix.scripts.cta_stage1.CommandLineOptions()
opt.set_run_number(0)
opt.set_o('calin_stage1.sqlite')
opt.set_db_results_table_name('stage1')
opt.set_log_frequency(cfg.log_frequency())
opt.set_nthread(1)
opt.mutable_decoder().CopyFrom(cfg.decoder())
opt.mutable_zfits().CopyFrom(cfg.zfits())
opt.mutable_zmq().CopyFrom(cfg.zmq())
opt.mutable_stage1().CopyFrom(calin.diagnostics.stage1.Stage1ParallelEventVisitor.default_config())

opt_proc = calin.util.options_processor.OptionsProcessor(opt, True);
opt_proc.process_arguments(sys.argv)

if(opt_proc.help_requested()):
    print('Usage:',opt_proc.program_name(),'[options] zmq_endpoint [zmq_endpoint...]')
    print('or:   ',opt_proc.program_name(),'[options] zfits_file [zfits_file...]\n')
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

# Create the stage1 visitor
s1pev = calin.diagnostics.stage1.Stage1ParallelEventVisitor(opt.stage1())
dispatcher.add_visitor(s1pev)

# Open SQL file
sql_mode = calin.io.sql_serializer.SQLite3Serializer.EXISTING_OR_NEW_RW
if opt.truncate_db():
    sql_mode = calin.io.sql_serializer.SQLite3Serializer.TRUNCATE_RW
sql = calin.io.sql_serializer.SQLite3Serializer(sql_file, sql_mode)
sql.create_or_extend_tables(opt.db_results_table_name(),
    calin.ix.diagnostics.stage1.Stage1.descriptor())

# Run all the visitors
if(endpoints[0].startswith('tcp://') or endpoints[0].startswith('ipc://')
        or endpoints[0].startswith('pgm://') or endpoints[0].startswith('pgme://')):
    dispatcher.process_cta_zmq_run(endpoints, cfg)
    sql.insert(opt.db_results_table_name(), visitor.stage1_results())
else:
    first_file = True
    failed_files = []
    nsuccess = 0;
    nskip = 0;
    for ifile, filename in enumerate(endpoints[opt.start_file_index():]):
        ifile += opt.start_file_index()
        if not first_file:
            print("-"*80)
        first_file = False
        if(opt.skip_existing() or opt.replace_existing()):
            selector = calin.ix.diagnostics.stage1.SelectByFilename()
            selector.set_filename(filename)
            oids = sql.select_oids_matching(opt.db_results_table_name(), selector)
            if oids and opt.skip_existing():
                print("#%d: skipping %s"%(ifile,filename))
                nskip += 1
                continue
            if oids and opt.replace_existing():
                for oid in oids:
                    sql.delete_by_oid(opt.db_results_table_name(), oid)
        print("#%d: processing %s"%(ifile,filename))
        try:
            calin.util.log.prune_default_protobuf_log()
            calin.provenance.chronicle.prune_the_chronicle()
            dispatcher.process_cta_zfits_run(filename, cfg)
            good, oid = sql.insert(opt.db_results_table_name(), s1pev.stage1_results())
            if(good):
                print("Inserted stage1 results into database with OID:",oid)
                nsuccess += 1
            else:
                print("Failed to insert stage1 results into database")
                failed_files.append(filename)
        except Exception as x:
            traceback.print_exception(*sys.exc_info())
            failed_files.append(filename)
            pass

    print("")
    print("="*80)
    if(nsuccess > 0):
        print("Successfully processed",nsuccess,"runs." if nsuccess!=1 else "run.")
    if(nskip > 0):
        print("Skipped",nskip,"runs" if nsuccess!=1 else "run","that were alreadin in database.")
    if(len(failed_files) > 0):
        print("Processing failed for",len(failed_files),"runs." if len(failed_files)!=1 else "run.")
        if(len(failed_files) > 1):
            print("Failed files :")
            for filename in failed_files:
                print("--", filename)
        else:
            print("Failed file :",failed_files[0])
    print("="*80)

# The end
