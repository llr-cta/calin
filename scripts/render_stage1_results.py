#!/usr/bin/env python3

# calin/scripts/render_stage1_results.py -- Stephen Fegan - 2021-07-05
#
# Render stage1 results from SQLITE3 to filesystem or Google drive
#
# Copyright 2021, Stephen Fegan <sfegan@llr.in2p3.fr>
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

import calin.io.sql_serializer
import calin.util.log
import calin.util.options_processor
import calin.ix.diagnostics.stage1
import calin.ix.scripts.render_stage1_results
import calin.provenance.chronicle
import calin.provenance.anthology
import calin.diagnostics.stage1_plotting

py_log = calin.util.log.PythonLogger()
py_log.this.disown()
calin.util.log.default_logger().add_logger(calin.util.log.default_protobuf_logger(),False)
calin.util.log.default_logger().add_logger(py_log,True)

opt = calin.ix.scripts.render_stage1_results.CommandLineOptions()
opt.set_db('calin_stage1.sqlite')
opt.set_db_stage1_table_name('stage1')
opt.set_base_directory('.')
opt.set_summary_csv('stage1_summary.csv')

opt_proc = calin.util.options_processor.OptionsProcessor(opt, True);
opt_proc.process_arguments(sys.argv)

if(opt_proc.help_requested()):
    print('Usage:',opt_proc.program_name(),'[options]')
    print('Options:\n')
    print(opt_proc.usage())
    exit(0)

if(len(opt_proc.arguments()) > 0):
    print('Unnecessary arguments supplied. Use "-help" option to get usage information.')
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
sql_file           = opt.db();

# Open SQL file
sql_mode = calin.io.sql_serializer.SQLite3Serializer.READ_ONLY
sql = calin.io.sql_serializer.SQLite3Serializer(sql_file, sql_mode)



# The end
