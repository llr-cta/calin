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
import time
import shutil
import os
import fcntl
import concurrent.futures 

import calin.ix.io.zmq_data_source
import calin.iact_data.event_dispatcher
import calin.io.sql_serializer
import calin.util.log
import calin.io.options_processor
import calin.diagnostics.stage1
import calin.ix.scripts.cta_stage1
import calin.provenance.chronicle


def init(local_opt, local_nfile):
    global py_log, opt, cfg, sql, dispatcher, s1cfg, s1pev, nfile

    py_log = calin.util.log.PythonLogger()
    py_log.this.disown()
    calin.util.log.default_logger().add_logger(calin.util.log.default_protobuf_logger(),False)
    calin.util.log.default_logger().add_logger(py_log,True)

    opt = local_opt

    cfg = calin.iact_data.event_dispatcher.ParallelEventDispatcher.default_config()
    cfg.set_run_number(opt.run_number())
    cfg.set_log_frequency(opt.log_frequency())
    cfg.set_nthread(opt.nthread())
    cfg.mutable_decoder().CopyFrom(opt.decoder())
    cfg.mutable_zfits().CopyFrom(opt.zfits())
    cfg.mutable_zmq().CopyFrom(opt.zmq())

    # Create the dispatcher
    dispatcher = calin.iact_data.event_dispatcher.ParallelEventDispatcher()

    # Create the stage1 visitor
    s1cfg = opt.const_stage1().Clone()
    if(opt.const_stage1().ancillary_database_directory() == ""):
        s1cfg.set_ancillary_database_directory(opt.copy_ancillary_db())
    s1pev = calin.diagnostics.stage1.Stage1ParallelEventVisitor(s1cfg)
    dispatcher.add_visitor(s1pev)

    # Open SQL file
    sql_file = opt.o();
    sql = calin.io.sql_serializer.SQLite3Serializer(sql_file, 
        calin.io.sql_serializer.SQLite3Serializer.EXISTING_RW)
    
    nfile = local_nfile

def run_file_analysis(ifile, filename):
    global sql, nfile, dispatcher, dispatcher, s1pev, s1cfg
    if(dispatcher is None):
        dispatcher, s1pev, s1cfg = make_dispatcher()

    if(sql is None):
        print("THERE")

    copied_ancillary_db = ''
    if(opt.skip_existing() or opt.replace_existing()):
        selector = calin.ix.diagnostics.stage1.SelectByFilename()
        selector.set_filename(filename)
        oids = sql.select_oids_matching(opt.db_stage1_table_name(), selector)
        if oids and opt.skip_existing():
            print("#%d / %d: skipping %s"%(ifile+1,nfile,filename))
            return [ 'skipped', copied_ancillary_db ]
        if oids and opt.replace_existing():
            for oid in oids:
                print(f"Deleting old stage1 results from database, OID : {oid}")
                sql.delete_by_oid(opt.db_stage1_table_name(), oid)

    print("#%d / %d: processing %s"%(ifile+1,nfile,filename))

    if(opt.copy_ancillary_db() != ''):
        lockfile = open(s1cfg.ancillary_database_directory()+"/.lockfile",'ab')
        fcntl.lockf(lockfile, fcntl.LOCK_EX)
        src_ancillary_db = calin.diagnostics.stage1.Stage1ParallelEventVisitor.nectarcam_ancillary_database_filename(filename,0,
            opt.const_stage1().ancillary_database(), opt.const_stage1().ancillary_database_directory())
        dst_ancillary_db = calin.diagnostics.stage1.Stage1ParallelEventVisitor.nectarcam_ancillary_database_filename(filename,0,
            s1cfg.ancillary_database(), s1cfg.ancillary_database_directory())
        if(os.path.isfile(src_ancillary_db)):
            try:
                if(not os.path.isfile(dst_ancillary_db)):
                    print("Copying %s -> %s"%(src_ancillary_db,dst_ancillary_db))
                    shutil.copyfile(src_ancillary_db, dst_ancillary_db)
                    copied_ancillary_db = dst_ancillary_db
            except Exception as x:
                traceback.print_exception(*sys.exc_info())
        else:
            print("*** Could not copy ancillary DB %s ***"%(src_ancillary_db))
        fcntl.lockf(lockfile, fcntl.LOCK_UN)
        lockfile.close()

    try:
        calin.util.log.prune_default_protobuf_log()
        calin.provenance.chronicle.prune_the_chronicle()
        dispatcher.process_cta_zfits_run(filename, cfg)
        s1res = s1pev.stage1_results()
        print(f"Inserting stage1 results into database, size: {s1res.SpaceUsedLong()/1024**2:,.1f} MB")
        start_time = time.time()
        good, oid = sql.insert(opt.db_stage1_table_name(), s1res)
        if(good):
            print(f"Inserted results into database in {time.time()-start_time:,.3f} sec, OID : {oid}")
            return [ 'success', copied_ancillary_db ]
        else:
            print("Failed to insert stage1 results into database")
            return [ 'failed', copied_ancillary_db ]
    except Exception as x:
        traceback.print_exception(*sys.exc_info())
        return [ 'failed', copied_ancillary_db ]    

def proc_file(ifile_filename):
    ifile, filename = ifile_filename
    return run_file_analysis(ifile, filename)

# Entry point for the program
if __name__ == '__main__':
    # Process command line options
    default_cfg = calin.iact_data.event_dispatcher.ParallelEventDispatcher.default_config()

    opt = calin.ix.scripts.cta_stage1.CommandLineOptions()
    opt.set_run_number(0)
    opt.set_o('calin_stage1.sqlite')
    opt.set_db_stage1_table_name('stage1')
    opt.set_log_frequency(default_cfg.log_frequency())
    opt.set_nthread(1)
    opt.mutable_decoder().CopyFrom(default_cfg.decoder())
    opt.mutable_zfits().CopyFrom(default_cfg.zfits())
    opt.mutable_zmq().CopyFrom(default_cfg.zmq())
    opt.mutable_stage1().CopyFrom(calin.diagnostics.stage1.Stage1ParallelEventVisitor.default_config())

    opt_proc = calin.io.options_processor.OptionsProcessor(opt, True);
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

    endpoints = opt_proc.arguments()

    # Create or extend SQL file
    sql_mode = calin.io.sql_serializer.SQLite3Serializer.EXISTING_OR_NEW_RW
    if opt.truncate_db():
        sql_mode = calin.io.sql_serializer.SQLite3Serializer.TRUNCATE_RW
    sql = calin.io.sql_serializer.SQLite3Serializer(opt.o(), sql_mode)
    sql.create_or_extend_tables(opt.db_stage1_table_name(),
        calin.ix.diagnostics.stage1.Stage1.descriptor())
    del sql

    # Run all the visitors
    if(endpoints[0].startswith('tcp://') or endpoints[0].startswith('ipc://')
            or endpoints[0].startswith('pgm://') or endpoints[0].startswith('pgme://')):
        init(opt,1)
        dispatcher.process_cta_zmq_run(endpoints, cfg)
        sql.insert(opt.db_stage1_table_name(), visitor.stage1_results())
    elif(opt.process_pool() >= 2):
        filelist = [ (ifile, filename) for ifile, filename in enumerate(endpoints[opt.start_file_index():]) ]
        nfile = len(filelist)
        with concurrent.futures.ProcessPoolExecutor(opt.process_pool(), initializer=init, initargs=(opt,nfile)) as executor:
            results = executor.map(proc_file, filelist)
        failed_files = []
        nsuccess = 0;
        nskip = 0;    
        for ifile_filename,result in zip(filelist,results):
            if(result[0] == 'success'):
                nsuccess += 1
            elif(result[0] == 'skipped'):
                nskip += 1
            elif(result[0] == 'failed'):
                failed_files.append(ifile_filename[1])
            if(result[1] != ''):
                print("Deleting %s"%(result[1]))
                os.unlink(result[1])
        print("="*80)
        if(nsuccess > 0):
            print("Successfully processed",nsuccess,"runs." if nsuccess!=1 else "run.")
        if(nskip > 0):
            print("Skipped",nskip,"runs that were" if nskip!=1 else "run that was","already in in database.")
        if(len(failed_files) > 0):
            print("Processing failed for",len(failed_files),"runs." if len(failed_files)!=1 else "run.")
            if(len(failed_files) > 1):
                print("Failed files :")
                for filename in failed_files:
                    print("--", filename)
            else:
                print("Failed file :",failed_files[0])
        print("="*80)    
    else:
        nfile = len(endpoints[opt.start_file_index():])
        init(opt, nfile)
        copied_ancillary_db = ''
        first_file = True
        failed_files = []
        nsuccess = 0;
        nskip = 0;
        for ifile, filename in enumerate(endpoints[opt.start_file_index():]):
            if not first_file:
                print("-"*80)
            first_file = False
            status, newly_copied_ancillary_db = run_file_analysis(ifile, filename)
            if(newly_copied_ancillary_db != ''):
                print("Deleting %s"%(copied_ancillary_db))
                os.unlink(copied_ancillary_db)
                copied_ancillary_db = newly_copied_ancillary_db
            if(status == 'success'):
                nsuccess += 1
            elif(status == 'skipped'):
                nskip += 1
            elif(status == 'failed'):
                failed_files.append(filename)

        if(copied_ancillary_db != ''):
            print("Deleting %s"%(copied_ancillary_db))
            os.unlink(copied_ancillary_db)

        print("")
        print("="*80)
        if(nsuccess > 0):
            print("Successfully processed",nsuccess,"runs." if nsuccess!=1 else "run.")
        if(nskip > 0):
            print("Skipped",nskip,"runs that were" if nsuccess!=1 else "run that was","already in in database.")
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