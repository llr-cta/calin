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

def good_tag(tag):
    return "\x1b[37;42;97;38;5;15;1m" + tag + "\x1b[0m"

def bad_tag(tag):
    return "\x1b[37;41;97;101;1m" + tag + "\x1b[0m"

def info_tag(tag):
    return "\x1b[37;46;97;1m" + tag + "\x1b[0m"

def init_dispatcher(local_opt, local_nfile):
    global py_log, opt, cfg, dispatcher, s1cfg, s1pev, nfile

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

    nfile = local_nfile

def init_writer(local_opt):
    global sql, opt

    opt = local_opt

    # Open SQL file
    sql = calin.io.sql_serializer.SQLite3Serializer(opt.o(), 
        calin.io.sql_serializer.SQLite3Serializer.EXISTING_RW)

def copy_ancillary_db(src_ancillary_db, dst_ancillary_db):
    if(not os.path.isfile(src_ancillary_db)):
        # Nothing we can do except print a warning if the destination doesn't exist
        if(not os.path.isfile(dst_ancillary_db)):
            print("*** Could not copy ancillary DB %s ***"%(src_ancillary_db))
        return False
       
    src_ancillary_db_size = os.path.getsize(src_ancillary_db)
    if(os.path.isfile(dst_ancillary_db) and
            os.path.getsize(dst_ancillary_db) == src_ancillary_db_size):
        # Destination file already exists and is the same size as the source. We outie.
        return False
    
    with open(dst_ancillary_db, 'ab') as dst_file:
        fcntl.lockf(dst_file, fcntl.LOCK_EX)

        if(os.path.getsize(dst_ancillary_db) == src_ancillary_db_size):
            # Destination file is the same size as the source. No need to copy.
            fcntl.lockf(dst_file, fcntl.LOCK_UN)
            return False
        
        print(f'Copying {src_ancillary_db} -> {dst_ancillary_db} ({src_ancillary_db_size/1024**2:,.1f} MB)')

        dst_file.truncate() # If file sizes don't match then re-copy the file
        with open(src_ancillary_db, 'rb') as src_file:
            shutil.copyfileobj(src_file, dst_file)
        fcntl.lockf(dst_file, fcntl.LOCK_UN)
        return True

def run_file_analysis(ifile, filename):
    global nfile, dispatcher, dispatcher, s1pev, s1cfg, opt
    if(dispatcher is None):
        dispatcher, s1pev, s1cfg = make_dispatcher()

    print("#%d / %d: processing %s"%(ifile+1,nfile,filename))

    copied_ancillary_db = ''
    if(opt.copy_ancillary_db() != ''):
        src_ancillary_db = calin.diagnostics.stage1.Stage1ParallelEventVisitor.nectarcam_ancillary_database_filename(filename,0,
            opt.const_stage1().ancillary_database(), opt.const_stage1().ancillary_database_directory())
        dst_ancillary_db = calin.diagnostics.stage1.Stage1ParallelEventVisitor.nectarcam_ancillary_database_filename(filename,0,
            s1cfg.ancillary_database(), s1cfg.ancillary_database_directory())
        if(copy_ancillary_db(src_ancillary_db, dst_ancillary_db)):
            copied_ancillary_db = dst_ancillary_db

    try:
        calin.util.log.prune_default_protobuf_log()
        calin.provenance.chronicle.prune_the_chronicle()
        dispatcher.process_cta_zfits_run(filename, cfg)
        s1res = s1pev.stage1_results()
        return filename, 'success', s1res, copied_ancillary_db
    except Exception as x:
        traceback.print_exception(*sys.exc_info())
        return filename, 'failed', None, copied_ancillary_db

def insert_stage1_results(s1res, filename):
    global sql, opt

    try:
        run_number_str = str(s1res.run_number()) if s1res.run_number() != 0 else '????'
        if(opt.replace_existing()):
            selector = calin.ix.diagnostics.stage1.SelectByFilename()
            selector.set_filename(filename)
            oids = sql.select_oids_matching(opt.db_stage1_table_name(), selector)
            if oids:
                for oid in oids:
                    print(f'{info_tag(" DATABASE ")} Run {run_number_str} deleting old stage1 results from database, OID : {oid}')
                    sql.delete_by_oid(opt.db_stage1_table_name(), oid)

        data_size = s1res.SpaceUsedLong()/1024**2
        print(f'{info_tag(" DATABASE ")} Run {run_number_str} inserting stage1 results into database, size: {data_size:,.1f} MB')
        start_time = time.time()
        good, oid = sql.insert(opt.db_stage1_table_name(), s1res)
        if(good):
            print(f'{good_tag(" DATABASE ")} Run {run_number_str} inserted {data_size:,.1f} MB into database in {time.time()-start_time:,.3f} sec, OID : {oid}')
        else:
            print(f'{bad_tag(" DATABASE ")} Run {run_number_str} failed to insert stage1 results into database')
        return filename, good
    except Exception as x:
        traceback.print_exception(*sys.exc_info())
        return filename, False

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
    endpoints = endpoints[opt.start_file_index():]

    # Create or extend SQL file
    sql_mode = calin.io.sql_serializer.SQLite3Serializer.EXISTING_OR_NEW_RW
    if opt.truncate_db():
        sql_mode = calin.io.sql_serializer.SQLite3Serializer.TRUNCATE_RW
    sql = calin.io.sql_serializer.SQLite3Serializer(opt.o(), sql_mode)
    sql.create_or_extend_tables(opt.db_stage1_table_name(),
        calin.ix.diagnostics.stage1.Stage1.descriptor())

    if(len(endpoints) == 0):
        print("No files to process.")
        exit(0)

    # Check if the files are already in the database
    if(opt.skip_existing()):
        filtered_endpoints = []
        for endpoint in endpoints:
            selector = calin.ix.diagnostics.stage1.SelectByFilename()
            selector.set_filename(endpoint)
            oids = sql.select_oids_matching(opt.db_stage1_table_name(), selector)
            if oids:
                print("Skipping %s"%(endpoint))
            elif endpoint in filtered_endpoints:
                print("Duplicate %s"%(endpoint))
            else:
                filtered_endpoints.append(endpoint)
        endpoints = filtered_endpoints

    all_copied_ancillary_db = []
    failed_files = []
    nsuccess = 0;
    nskip = 0;

    # Run all the visitors
    if(endpoints[0].startswith('tcp://') or endpoints[0].startswith('ipc://')
            or endpoints[0].startswith('pgm://') or endpoints[0].startswith('pgme://')):
        init_dispatcher(opt,1)
        dispatcher.process_cta_zmq_run(endpoints, cfg)
        s1res = s1pev.stage1_results()
        sql.insert(opt.db_stage1_table_name(), s1res)
    elif(opt.process_pool() >= 1):
        filelist = [ (ifile, filename) for ifile, filename in enumerate(endpoints) ]
        nfile = len(filelist)
        with concurrent.futures.ProcessPoolExecutor(1, initializer=init_writer, initargs=(opt,)) as writer_executor:
            writer_futures = []
            with concurrent.futures.ProcessPoolExecutor(opt.process_pool(), initializer=init_dispatcher, initargs=(opt,nfile)) as dispatcher_executor:
                dispatcher_futures = [ dispatcher_executor.submit(run_file_analysis, *ifile_filename) for ifile_filename in filelist ]
                for future in concurrent.futures.as_completed(dispatcher_futures):
                    filename, status, s1res, copied_ancillary_db = future.result()
                    if(status == 'success' and s1res.ByteSize() > 0):
                        writer_futures.append(writer_executor.submit(insert_stage1_results, s1res, filename))
                    elif(status == 'success'):
                        print(f'bad_tag(" PROCESS POOL ERROR ") No stage1 results received for {filename}')
                        failed_files.append(filename)
                    else:
                        failed_files.append(filename)
                    
                    if(copied_ancillary_db != ''):
                        all_copied_ancillary_db.append(copied_ancillary_db)

            for future in concurrent.futures.as_completed(writer_futures):
                filename, good = future.result()
                if(good):
                    nsuccess += 1
                else:
                    failed_files.append(filename)
    else:
        init_dispatcher(opt, len(endpoints))
        copied_ancillary_db = ''
        first_file = True
        for ifile, filename in enumerate(endpoints):
            if not first_file:
                print("-"*80)
            first_file = False
            _, status, s1res, copied_ancillary_db = run_file_analysis(ifile, filename)
            if(status == 'success'):
                _, good = insert_stage1_results(s1res, filename)
                if(good):
                    nsuccess += 1
                else:
                    failed_files.append(filename)
            else:
                failed_files.append(filename)

            if(copied_ancillary_db != ''):
                all_copied_ancillary_db.append(copied_ancillary_db)

    # for copied_ancillary_db in all_copied_ancillary_db:
    #     print("Deleting %s"%copied_ancillary_db)
    #     os.unlink(copied_ancillary_db)

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