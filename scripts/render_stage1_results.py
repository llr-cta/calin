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
import concurrent.futures

import calin.io.sql_serializer
import calin.util.log
import calin.util.options_processor
import calin.ix.diagnostics.stage1
import calin.ix.scripts.render_stage1_results
import calin.provenance.chronicle
import calin.provenance.anthology
import calin.diagnostics.stage1_plotting
import calin.iact_data.nectarcam_layout
import calin.io.uploader

import matplotlib
import matplotlib.figure
import matplotlib.backends.backend_agg

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
opt.set_overwrite(True)

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

sql_file = opt.db();
sql_mode = calin.io.sql_serializer.SQLite3Serializer.READ_ONLY

figure_dpi = max(opt.figure_dpi(), 60)

def get_oids():
    # Open SQL file
    sql = calin.io.sql_serializer.SQLite3Serializer(sql_file, sql_mode)

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

    return all_oid

def render_oid(oid):
    # Open SQL file
    sql = calin.io.sql_serializer.SQLite3Serializer(sql_file, sql_mode)

    stage1 = calin.ix.diagnostics.stage1.Stage1()
    sql.retrieve_by_oid(opt.db_stage1_table_name(), oid, stage1)
    runno = stage1.run_number()
    print('Run :', runno)

    if(opt.upload_to_google_drive()):
        uploader = calin.io.uploader.GoogleDriveUploader(opt.google().token_file(),
            opt.google().base_directory(), opt.google().credentials_file(),
            overwrite=opt.overwrite(), loud=opt.loud_upload())
    else:
        uploader = calin.io.uploader.FilesystemUploader(opt.base_directory(),
            overwrite=opt.overwrite(), loud=opt.loud_upload())

    if(opt.force_nectarcam_61_camera()):
        cast_to_nectarcam_61_camera(stage1)

    def upload_figure(runno, quantity, f):
        runbatchpath = 'runs%d-%d'%(int(runno/1000)*1000, (int(runno/1000)+1)*1000)
        filenames = [ \
            'by run/%s/run%d/run%d_%s.png'%(runbatchpath, runno, runno, quantity),
            'by quantity/%s/%s/run%d_%s.png'%(quantity, runbatchpath, runno, quantity) ]
        uploader.upload_png_from_figure(filenames, f)

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

    ############################################################################
    # FIGURE : L0 trigger frequency
    ############################################################################

    if(stage1.const_charge_stats().channel_triggered_count_size()>0):
        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        calin.diagnostics.stage1_plotting.draw_trigger_event_fraction(stage1, axis=ax)
        ax.set_title('L0 trigger bit frequency, run : %d'%runno)
        upload_figure(runno, 'trigger_l0_bit_frequency', ax.figure)

    if(stage1.const_charge_stats().phy_trigger_few_neighbor_channel_triggered_count_size()>0):
        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        calin.diagnostics.stage1_plotting.draw_nn_failed_phy_trigger_event_fraction(stage1, axis=ax)
        ax.set_title('L0 trigger bit in sub-3NN events, run : %d'%runno)
        upload_figure(runno, 'trigger_sub_3nn_l0_bit_frequency', ax.figure)

    if(stage1.const_charge_stats().const_num_channel_triggered_hist().sum_w()):
        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        calin.diagnostics.stage1_plotting.draw_num_channel_triggered_hist(stage1,axis=ax)
        ax.set_title('L0 trigger bit histogram (all events), run : %d'%runno)
        upload_figure(runno, 'trigger_l0_bit_count', ax.figure)

    if(stage1.const_charge_stats().const_phy_trigger_num_channel_triggered_hist().sum_w()):
        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        calin.diagnostics.stage1_plotting.draw_num_channel_triggered_hist(stage1,axis=ax,
            phys_trigger=True)
        ax.set_title('L0 trigger bit histogram (physics events), run : %d'%runno)
        upload_figure(runno, 'trigger_l0_bit_count_phys', ax.figure)

        ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
        calin.diagnostics.stage1_plotting.draw_num_channel_triggered_hist(stage1,axis=ax,
            phys_trigger=True,zoom=True)
        ax.set_title('L0 trigger bit histogram (physics events), run : %d'%runno)
        upload_figure(runno, 'trigger_l0_bit_count_phys_zoom', ax.figure)

    ############################################################################
    # FIGURE : mean waveforms
    ############################################################################

    if(stage1.has_mean_wf_pedestal() and stage1.has_charge_stats()):
        if(stage1.const_mean_wf_pedestal().channel_high_gain_size() \
                and stage1.const_mean_wf_pedestal().has_camera_high_gain() \
                and stage1.const_charge_stats().has_high_gain() \
                and numpy.count_nonzero(stage1.const_charge_stats().const_high_gain().ped_trigger_event_count())):
            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_mean_wf(stage1,axis=ax,dataset='pedestal')
            ax.set_title('Mean waveform (pedestal, high-gain), run : %d'%runno)
            upload_figure(runno, 'waveform_mean_pedestal_hg', ax.figure)

            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_mean_wf(stage1,axis=ax,dataset='pedestal',
                subtract_pedestal=True)
            ax.set_title('Mean waveform offset (pedestal, high-gain), run : %d'%runno)
            upload_figure(runno, 'waveform_mean_pedestal_hg_offset', ax.figure)

            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_mean_wf_deviation_from_camera_mean(stage1,
                axis=ax,dataset='pedestal')
            ax.set_title('RMS waveform offset (pedestal, high-gain), run : %d'%runno)
            upload_figure(runno, 'waveform_mean_pedestal_hg_offset_rms', ax.figure)

        if(stage1.const_mean_wf_pedestal().channel_low_gain_size() \
                and stage1.const_mean_wf_pedestal().has_camera_low_gain() \
                and stage1.const_charge_stats().has_low_gain() \
                and numpy.count_nonzero(stage1.const_charge_stats().const_low_gain().ped_trigger_event_count())):
            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_mean_wf(stage1,axis=ax,dataset='pedestal',low_gain=True)
            ax.set_title('Mean waveform (pedestal, low-gain), run : %d'%runno)
            upload_figure(runno, 'waveform_mean_pedestal_lg', ax.figure)

            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_mean_wf(stage1,axis=ax,dataset='pedestal',low_gain=True,
                subtract_pedestal=True)
            ax.set_title('Mean waveform offset (pedestal, low-gain), run : %d'%runno)
            upload_figure(runno, 'waveform_mean_pedestal_lg_offset', ax.figure)

            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_mean_wf_deviation_from_camera_mean(stage1,
                axis=ax, dataset='pedestal', low_gain=True)
            ax.set_title('RMS waveform offset (pedestal, low-gain), run : %d'%runno)
            upload_figure(runno, 'waveform_mean_pedestal_lg_offset_rms', ax.figure)

    if(stage1.has_mean_wf_physics() and stage1.has_charge_stats()):
        if(stage1.const_mean_wf_physics().channel_high_gain_size() \
                and stage1.const_mean_wf_physics().has_camera_high_gain() \
                and stage1.const_charge_stats().has_high_gain() \
                and numpy.count_nonzero(stage1.const_charge_stats().const_high_gain().phy_trigger_event_count())):
            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_mean_wf(stage1,axis=ax,dataset='physics')
            ax.set_title('Mean waveform (physics, high-gain), run : %d'%runno)
            upload_figure(runno, 'waveform_mean_physics_hg', ax.figure)
        if(stage1.const_mean_wf_physics().channel_low_gain_size() \
                and stage1.const_mean_wf_physics().has_camera_low_gain() \
                and stage1.const_charge_stats().has_low_gain() \
                and numpy.count_nonzero(stage1.const_charge_stats().const_low_gain().phy_trigger_event_count())):
            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_mean_wf(stage1,axis=ax,dataset='physics',low_gain=True)
            ax.set_title('Mean waveform (physics, low-gain), run : %d'%runno)
            upload_figure(runno, 'waveform_mean_physics_lg', ax.figure)

    if(stage1.has_mean_wf_external_flasher() and stage1.has_charge_stats()):
        if(stage1.const_mean_wf_external_flasher().channel_high_gain_size() \
                and stage1.const_mean_wf_external_flasher().has_camera_high_gain() \
                and stage1.const_charge_stats().has_high_gain() \
                and numpy.count_nonzero(stage1.const_charge_stats().const_high_gain().ext_trigger_event_count())):
            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_mean_wf(stage1,axis=ax,dataset='external_flasher')
            ax.set_title('Mean waveform (ext flasher, high-gain), run : %d'%runno)
            upload_figure(runno, 'waveform_mean_external_flasher_hg', ax.figure)
        if(stage1.const_mean_wf_external_flasher().channel_low_gain_size() \
                and stage1.const_mean_wf_external_flasher().has_camera_low_gain() \
                and stage1.const_charge_stats().has_low_gain() \
                and numpy.count_nonzero(stage1.const_charge_stats().const_low_gain().ext_trigger_event_count())):
            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_mean_wf(stage1,axis=ax,dataset='external_flasher',low_gain=True)
            ax.set_title('Mean waveform (ext flasher, low-gain), run : %d'%runno)
            upload_figure(runno, 'waveform_mean_external_flasher_lg', ax.figure)

    if(stage1.has_mean_wf_internal_flasher() and stage1.has_charge_stats()):
        if(stage1.const_mean_wf_internal_flasher().channel_high_gain_size() \
                and stage1.const_mean_wf_internal_flasher().has_camera_high_gain() \
                and stage1.const_charge_stats().has_high_gain() \
                and numpy.count_nonzero(stage1.const_charge_stats().const_high_gain().int_trigger_event_count())):
            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_mean_wf(stage1,axis=ax,dataset='internal_flasher')
            ax.set_title('Mean waveform (int flasher, high-gain), run : %d'%runno)
            upload_figure(runno, 'waveform_mean_internal_flasher_hg', ax.figure)
        if(stage1.const_mean_wf_internal_flasher().channel_low_gain_size() \
                and stage1.const_mean_wf_internal_flasher().has_camera_low_gain() \
                and stage1.const_charge_stats().has_low_gain() \
                and numpy.count_nonzero(stage1.const_charge_stats().const_low_gain().int_trigger_event_count())):
            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_mean_wf(stage1,axis=ax,dataset='internal_flasher',low_gain=True)
            ax.set_title('Mean waveform (int flasher, low-gain), run : %d'%runno)
            upload_figure(runno, 'waveform_mean_internal_flasher_lg', ax.figure)

    return True

all_oid = get_oids()
all_status = []

if(opt.nthread()>1):
    with concurrent.futures.ProcessPoolExecutor(max_workers=opt.nthread()) as executor:
        results = executor.map(render_oid, all_oid)
        for status in results:
            try:
                if(status):
                    all_status.append(True)
            except Exception as e:
                traceback.print_exception(*sys.exc_info())
                all_status.append(False)
else:
    for oid in all_oid:
        try:
            good = render_oid(oid)
            all_status.append(good)
        except Exception as e:
            traceback.print_exception(*sys.exc_info())
            all_status.append(False)

print("=================================== RESULTS ===================================")
sql = calin.io.sql_serializer.SQLite3Serializer(sql_file, sql_mode)
for oid_status in zip(all_oid, all_status):
    stage1pod = calin.ix.diagnostics.stage1.Stage1POD()
    sql.retrieve_by_oid(opt.db_stage1_table_name(), oid_status[0], stage1pod)
    print(stage1pod.run_number(), "success" if oid_status[1] else "*** FAILED ***")

# The end
