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
import calin.iact_data.nectarcam_layout

import matplotlib
import matplotlib.figure
import matplotlib.backends.backend_agg
import io
import os

class FilesystemUploader:
    def __init__(self, root_directory):
        self.root_directory = os.path.normpath(os.path.expanduser(root_directory)) if root_directory else '.'
        if(not os.path.isdir(self.root_directory)):
            raise RuntimeError('Base path os not directory : '+self.root_directory)

    def make_path(self, rel_path):
        if(not rel_path):
            return ''
        rel_path = os.path.normpath(rel_path)
        abs_path = os.path.normpath(os.path.join(self.root_directory, rel_path))
        if((self.root_directory == '.' and (abs_path.startswith('../') or abs_path.startswith('/')))
                or (self.root_directory != '.' and not abs_path.startswith(self.root_directory))):
            raise RuntimeError('Cannot make path outside of base : '+rel_path)
        if(not os.path.isdir(abs_path)):
            (head, tail) = os.path.split(rel_path)
            self.make_path(head)
            # print("mkdir",abs_path)
            os.mkdir(abs_path)
        return abs_path

    def upload_from_io(self, rel_filepath, iostream):
        (rel_path, filename) = os.path.split(rel_filepath)
        abs_path = os.path.join(self.make_path(rel_path), filename)
        mode = 'wb' if iostream is io.StringIO else 'w'
        with open(abs_path, mode) as f:
            f.write(iostream.getvalue())

    def upload_png_from_figure(self, rel_filepath, figure):
        (rel_path, filename) = os.path.split(rel_filepath)
        abs_path = os.path.join(self.make_path(rel_path), filename)
        matplotlib.backends.backend_agg.FigureCanvasAgg(figure).print_png(abs_path)

def cast_to_nectarcam_61_camera(stage1):
    if(stage1.const_run_config().configured_module_id_size() != 61):
        raise RuntimeError('cast_to_nectarcam_61_camera : run must have 61 modules')

    chan_index = stage1.const_run_config().configured_channel_index()
    mod_index = stage1.const_run_config().configured_module_index()

    ad = stage1.nectarcam().ancillary_data().Clone()
    stage1.mutable_nectarcam().mutable_ancillary_data().clear_feb_temperature()
    for imod in ad.feb_temperature_keys():
        ofm = ad.feb_temperature(imod)
        nfm = stage1.mutable_nectarcam().mutable_ancillary_data().mutable_feb_temperature(int(mod_index[imod]))
        for im in range(ofm.measurement_size()):
            msmt = nfm.add_measurement()
            msmt.CopyFrom(ofm.measurement(im))
            msmt.set_drawer(int(mod_index[msmt.drawer()]))

    stage1.mutable_nectarcam().mutable_ancillary_data().clear_hvpa_voltage()
    for ichan in ad.hvpa_voltage_keys():
        ovm = ad.hvpa_voltage(ichan)
        nvm = stage1.mutable_nectarcam().mutable_ancillary_data().mutable_hvpa_voltage(int(chan_index[ichan]))
        for im in range(ovm.measurement_size()):
            msmt = nvm.add_measurement()
            msmt.CopyFrom(ovm.measurement(im))
            msmt.set_drawer(int(mod_index[msmt.drawer()]))

    stage1.mutable_nectarcam().mutable_ancillary_data().clear_hvpa_current()
    for ichan in ad.hvpa_current_keys():
        ocm = ad.hvpa_current(ichan)
        ncm = stage1.mutable_nectarcam().mutable_ancillary_data().mutable_hvpa_current(int(chan_index[ichan]))
        for im in range(ocm.measurement_size()):
            msmt = ncm.add_measurement()
            msmt.CopyFrom(ocm.measurement(im))
            msmt.set_drawer(int(mod_index[msmt.drawer()]))

    stage1.mutable_run_config().clear_camera_layout()
    stage1.mutable_run_config().mutable_camera_layout().CopyFrom(calin.iact_data.nectarcam_layout.nectarcam_61module_layout())
    stage1.mutable_run_config().clear_configured_channel_id()
    stage1.mutable_run_config().set_configured_channel_id(numpy.asarray(range(61*7),dtype=numpy.int32))
    stage1.mutable_run_config().clear_configured_channel_index()
    stage1.mutable_run_config().set_configured_channel_index(numpy.asarray(range(61*7),dtype=numpy.int32))

    stage1.mutable_run_config().clear_configured_module_id()
    stage1.mutable_run_config().set_configured_module_id(numpy.asarray(range(61),dtype=numpy.int32))
    stage1.mutable_run_config().clear_configured_module_index()
    stage1.mutable_run_config().set_configured_module_index(numpy.asarray(range(61),dtype=numpy.int32))

py_log = calin.util.log.PythonLogger()
py_log.this.disown()
calin.util.log.default_logger().add_logger(calin.util.log.default_protobuf_logger(),False)
calin.util.log.default_logger().add_logger(py_log,True)

opt = calin.ix.scripts.render_stage1_results.CommandLineOptions()
opt.set_db('calin_stage1.sqlite')
opt.set_db_stage1_table_name('stage1')
opt.set_base_directory('.')
opt.set_summary_csv('stage1_summary.csv')
opt.set_figure_dpi(200)

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

# Open SQL file
sql_file = opt.db();
sql_mode = calin.io.sql_serializer.SQLite3Serializer.READ_ONLY
sql = calin.io.sql_serializer.SQLite3Serializer(sql_file, sql_mode)

# Uploader
uploader = FilesystemUploader(opt.base_directory())

sql_where = ''
if(opt.run_number() > 0):
    sql_where = 'WHERE run_number=%d'%opt.run_number()
elif(opt.from_run_number() > 0 and opt.to_run_number() > 0):
    sql_where = 'WHERE run_number>=%d AND run_number<=%d'%(
        opt.from_run_number(), opt.to_run_number())
elif(opt.from_run_number() > 0):
    sql_where = 'WHERE run_number>=%d'%opt.from_run_number()
elif(opt.to_run_number() > 0):
    sql_where = 'WHERE run_number<=%d'%opt.to_run_number()

if(sql_where):
    all_oid = sql.select_oids_by_sql(opt.db_stage1_table_name(), sql_where)
else:
    all_oid = sql.select_all_oids(opt.db_stage1_table_name())


figure_dpi = max(opt.figure_dpi(), 60)

for oid in all_oid:
    stage1 = calin.ix.diagnostics.stage1.Stage1()
    sql.retrieve_by_oid(opt.db_stage1_table_name(), oid, stage1)
    runno = stage1.run_number()
    print('Run :', runno)

    if(opt.force_nectarcam_61_camera()):
        cast_to_nectarcam_61_camera(stage1)

    def upload_figure(runno, quantity, f):
        runbatchpath = 'runs%d-%d'%(int(runno/1000)*1000, (int(runno/1000)+1)*1000)
        uploader.upload_png_from_figure(
            'by run/%s/run%d/run%d_%s.png'%(runbatchpath, runno, runno, quantity), f)
        uploader.upload_png_from_figure(
            'by quantity/%s/%s/run%d_%s.png'%(quantity, runbatchpath, runno, quantity), f)

    ############################################################################
    # FIGURE : Missing components
    ############################################################################

    ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
    calin.diagnostics.stage1_plotting.draw_missing_components_fraction(stage1,
        axis=ax, mod_label_fontsize=4, aux_label_fontsize=5.5)
    ax.set_title('Missing components, run : %d'%runno)
    upload_figure(runno, 'missing_components', ax.figure)

    ############################################################################
    # FIGURE : Pedestals
    ############################################################################

    if(stage1.has_charge_stats() and stage1.const_charge_stats().has_high_gain()
            and max(stage1.const_charge_stats().const_high_gain().ped_trigger_event_count())>0):
        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        calin.diagnostics.stage1_plotting.draw_pedestal_value(stage1,all_events_ped_win=False,low_gain=False, axis=ax)
        ax.set_title('High-gain pedestal mean (ped events), run : %d'%runno)
        upload_figure(runno, 'pedestal_mean_hg_ped_evt', ax.figure)

        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        calin.diagnostics.stage1_plotting.draw_pedestal_rms(stage1,all_events_ped_win=False,low_gain=False, axis=ax)
        ax.set_title('High-gain pedestal rms (ped events), run : %d'%runno)
        upload_figure(runno, 'pedestal_rms_hg_ped_evt', ax.figure)

    if(stage1.has_charge_stats() and stage1.const_charge_stats().has_low_gain()
            and max(stage1.const_charge_stats().const_low_gain().ped_trigger_event_count())>0):
        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        calin.diagnostics.stage1_plotting.draw_pedestal_value(stage1,all_events_ped_win=False,low_gain=True, axis=ax)
        ax.set_title('Low-gain pedestal mean (ped events), run : %d'%runno)
        upload_figure(runno, 'pedestal_mean_lg_ped_evt', ax.figure)

        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        calin.diagnostics.stage1_plotting.draw_pedestal_rms(stage1,all_events_ped_win=False,low_gain=True, axis=ax)
        ax.set_title('Low-gain pedestal rms (ped events), run : %d'%runno)
        upload_figure(runno, 'pedestal_rms_lg_ped_evt', ax.figure)

    if(stage1.has_charge_stats() and stage1.const_charge_stats().has_high_gain()
            and max(stage1.const_charge_stats().const_high_gain().all_trigger_event_count())>0):
        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        calin.diagnostics.stage1_plotting.draw_pedestal_value(stage1,all_events_ped_win=True,low_gain=False, axis=ax)
        ax.set_title('High-gain pedestal mean (all events), run : %d'%runno)
        upload_figure(runno, 'pedestal_mean_hg_all_evt', ax.figure)

        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        calin.diagnostics.stage1_plotting.draw_pedestal_rms(stage1,all_events_ped_win=True,low_gain=False, axis=ax)
        ax.set_title('High-gain pedestal rms (all events), run : %d'%runno)
        upload_figure(runno, 'pedestal_rms_hg_all_evt', ax.figure)

    if(stage1.has_charge_stats() and stage1.const_charge_stats().has_low_gain()
            and max(stage1.const_charge_stats().const_low_gain().all_trigger_event_count())>0):
        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        calin.diagnostics.stage1_plotting.draw_pedestal_value(stage1,all_events_ped_win=True,low_gain=True, axis=ax)
        ax.set_title('Low-gain pedestal mean (all events), run : %d'%runno)
        upload_figure(runno, 'pedestal_mean_lg_all_evt', ax.figure)

        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        calin.diagnostics.stage1_plotting.draw_pedestal_rms(stage1,all_events_ped_win=True,low_gain=True, axis=ax)
        ax.set_title('Low-gain pedestal rms (all events), run : %d'%runno)
        upload_figure(runno, 'pedestal_rms_lg_all_evt', ax.figure)

    ############################################################################
    # FIGURE : FEB temperatures
    ############################################################################

    if(stage1.has_nectarcam() and stage1.const_nectarcam().has_ancillary_data() and
            len(stage1.const_nectarcam().const_ancillary_data().feb_temperature_keys())>0):
        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        calin.diagnostics.stage1_plotting.draw_nectarcam_feb_temperatures(stage1,1,axis=ax)
        ax.set_title('Temperature (FEB 1), run : %d'%runno)
        upload_figure(runno, 'temperature_1', ax.figure)

        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        calin.diagnostics.stage1_plotting.draw_nectarcam_feb_temperatures(stage1,2,axis=ax)
        ax.set_title('Temperature (FEB 2), run : %d'%runno)
        upload_figure(runno, 'temperature_2', ax.figure)

        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        calin.diagnostics.stage1_plotting.draw_nectarcam_feb_temperatures_minmax(stage1,1,axis=ax)
        ax.set_title('Temperature spread (FEB 1), run : %d'%runno)
        upload_figure(runno, 'temperature_spread_1', ax.figure)

        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        calin.diagnostics.stage1_plotting.draw_nectarcam_feb_temperatures_minmax(stage1,2,axis=ax)
        ax.set_title('Temperature spread (FEB 2), run : %d'%runno)
        upload_figure(runno, 'temperature_spread_2', ax.figure)

    ############################################################################
    # FIGURE : clock regression
    ############################################################################

    iclock = 0
    if(stage1.has_clock_regression() and \
            stage1.const_clock_regression().module_clock_size()>=iclock and \
            stage1.const_clock_regression().const_module_clock(iclock).modules_size() > 0):
        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        ax2 = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        ax3 = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        ax4 = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        ax5 = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)

        calin.diagnostics.stage1_plotting.draw_all_clock_regression(stage1,ax,ax2,ax3,ax4,ax5)

        ax.set_title('Clock frequency error, run : %d'%runno)
        upload_figure(runno, 'clock_frequency_error', ax.figure)

        ax2.set_title('Clock offset from UCTS, run : %d'%runno)
        upload_figure(runno, 'clock_offset', ax2.figure)

        ax3.set_title('Clock vs UCTS fit RMS residual, run : %d'%runno)
        upload_figure(runno, 'clock_residual', ax3.figure)

        ax4.set_title('Clock frequency spread, run : %d'%runno)
        upload_figure(runno, 'clock_frequency_spread', ax4.figure)

        ax5.set_title('Clock offset from UCTS spread, run : %d'%runno)
        upload_figure(runno, 'clock_offset_spread', ax5.figure)

    ############################################################################
    # FIGURE : channel and module data order
    ############################################################################

    ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
    calin.diagnostics.stage1_plotting.draw_module_dataorder(stage1, axis=ax)
    ax.set_title('Module data ordering, run : %d'%runno)
    upload_figure(runno, 'data_ordering_module', ax.figure)

    ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
    calin.diagnostics.stage1_plotting.draw_channel_dataorder(stage1, axis=ax)
    ax.set_title('Channel data ordering, run : %d'%runno)
    upload_figure(runno, 'data_ordering_channel', ax.figure)

    ############################################################################
    # FIGURE : FPM voltage and current
    ############################################################################

    if(stage1.has_nectarcam() and stage1.const_nectarcam().has_ancillary_data() and
            len(stage1.const_nectarcam().const_ancillary_data().hvpa_voltage_keys())>0):
        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        ax2 = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        ax3 = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        ax4 = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        ax5 = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        ax6 = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)

        calin.diagnostics.stage1_plotting.draw_nectarcam_fpm_measurements(stage1,ax,ax2,ax3,ax4,ax5,ax6)

        ax.set_title('HVPA mean voltage, run : %d'%runno)
        upload_figure(runno, 'hvpa_voltage_mean', ax.figure)

        ax2.set_title('HVPA Cockroft-Walton current, run : %d'%runno)
        upload_figure(runno, 'hvpa_cw_current_mean', ax2.figure)

        ax3.set_title('HVPA board current, run : %d'%runno)
        upload_figure(runno, 'hvpa_board_current_mean', ax3.figure)

        ax4.set_title('HVPA voltage spread, run : %d'%runno)
        upload_figure(runno, 'hvpa_voltage_spread', ax4.figure)

        ax5.set_title('HVPA Cockroft-Walton current spread, run : %d'%runno)
        upload_figure(runno, 'hvpa_cw_current_spread', ax5.figure)

        ax6.set_title('HVPA board current spread, run : %d'%runno)
        upload_figure(runno, 'hvpa_board_current_spread', ax6.figure)


# The end