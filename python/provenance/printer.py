# calin/python/diagnostics/printer.py -- Stephen Fegan -- 2021-08-18
#
# Print the provenance
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


import numpy
import contextlib
import io

import calin.provenance.chronicle
import calin.ix.provenance.anthology
import calin.ix.provenance.chronicle
import calin.provenance.system_info

def print_provenance(writer, anthology):
    with contextlib.redirect_stdout(writer):
        c = anthology.const_chronicle()

        if(c.command_line_record_size()>0):
            cl = c.const_command_line_record(0)
            print('='*80)
            print('Command line')
            print('='*80)
            print('')
            for i in range(cl.arguments_size()):
                if(i==0 and cl.arguments_size==1):
                    print(cl.arguments(0))
                elif(i==0):
                    print(cl.arguments(0),'\\')
                elif(i != cl.arguments_size()-1):
                    print("  ",cl.arguments(i),'\\')
                else:
                    print("  ",cl.arguments(i))
            print('')

        l = anthology.const_default_log()
        print('='*80)
        print('Log')
        print('='*80)
        print('')
        for i in range(l.message_size()):
            m = l.const_message(i)
            for s in m.message().split('\n'):
                print(f'{calin.ix.util.log.LogMessage_Level_Name(m.level())} : {s}')
        print('')

        if(c.file_io_record_size()):
            print('='*80)
            print('Files opened')
            print('='*80)
            print('')
            for i in range(c.file_io_record_size()):
                f = c.const_file_io_record(i)
                print(f'{i} : {f.file_name()} ({f.const_open_timestamp().printable()})')
            print('')

        for ip in range(c.processing_record_size()):
            print('='*80)
            p = c.const_processing_record(ip)
            if(c.processing_record_size()>1):
                print(f'Processing {ip+1} ({p.type()})')
            else:
                print(f'Processing ({p.type()})')
            print('='*80)
            print('')
            dur = numpy.floor((p.close_timestamp().nsec()-p.open_timestamp().nsec())/1000000)/1000
            if(p.comment()):
                print(f'{p.description()} - {p.comment()}, (duration {dur} sec)')
            else:
                print(f'{p.description()} (duration {dur} sec)')
            for isp in range(p.subprocessing_record_size()):
                sp = p.const_subprocessing_record(isp)
                if(sp.comment()):
                    print(f'|-> {sp.description()} - {sp.comment()}')
                else:
                    print(f'|-> {sp.description()}')
            print('')

        print('='*80)
        print('Host system')
        print('='*80)
        print('')

        si = anthology.const_host_info()
        print("Host name :",si.host_name()," Host user:",si.user_name())
        print("System :",calin.provenance.system_info.system_info_string(si))
        print("Working directory :",si.current_working_directory())
        print('')

        print('='*80)
        print('Build information')
        print('='*80)
        print('')

        bi = anthology.const_build_info()
        print(calin.provenance.system_info.build_info_string(bi))
        print('')

        if(c.command_line_record_size()>0):
            cl = c.const_command_line_record(-1)
            for iclpo in range(cl.processed_options_size()):
                clpo = cl.const_processed_options(iclpo)
                print('='*80)
                if(clpo.type()):
                    print('Configuration ('+clpo.type()+')')
                else:
                    print('Configuration')
                print('='*80)
                print('')
                print(clpo.json())
                print('')
