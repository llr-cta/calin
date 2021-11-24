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
import io

import calin.io.sql_serializer
import calin.util.log
import calin.io.options_processor
import calin.ix.diagnostics.stage1
import calin.ix.scripts.render_stage1_results
import calin.provenance.chronicle
import calin.provenance.anthology
import calin.provenance.printer
import calin.diagnostics.stage1_plotting
import calin.diagnostics.stage1_summary
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
opt.set_summary_sheet('stage1_summary.csv')
opt.set_run_log_sheet('logsheet.csv')
opt.set_figure_dpi(200)
opt.set_overwrite(True)

opt_proc = calin.io.options_processor.OptionsProcessor(opt, True);
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


def new_uploader():
    if(opt.upload_to_google_drive()):
        uploader = calin.io.uploader.GoogleDriveUploader(opt.const_google().token_file(),
            opt.base_directory(), opt.const_google().credentials_file(),
            overwrite=opt.overwrite(), loud=opt.loud_upload())
    else:
        uploader = calin.io.uploader.FilesystemUploader(opt.base_directory(),
            overwrite=opt.overwrite(), loud=opt.loud_upload())
    return uploader

def get_oids():
    all_oid = []
    all_filename = []
    all_runno = []

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

    stage1pod = calin.ix.diagnostics.stage1.Stage1POD()
    for oid in all_oid:
        sql.retrieve_by_oid(opt.db_stage1_table_name(), oid, stage1pod)
        all_runno.append(stage1pod.run_number())
        all_filename.append(stage1pod.filename())

    return all_oid, all_runno, all_filename

def render_oid(oid):
    try:
        sql = calin.io.sql_serializer.SQLite3Serializer(sql_file, sql_mode)
        stage1 = calin.ix.diagnostics.stage1.Stage1()
        sql.retrieve_by_oid(opt.db_stage1_table_name(), oid, stage1)
        return unprotected_render_oid(stage1)
    except:
        runno = stage1.run_number()
        print('Exception while rendering run :', runno)
        raise

def unprotected_render_oid(stage1):
    runno = stage1.run_number()
    print('Started run :', runno)

    uploader = new_uploader()
    figure_factory = calin.plotting.MatplotlibFigureFactory(dpi=figure_dpi)

    if(opt.force_nectarcam_61_camera()):
        cast_to_nectarcam_61_camera(stage1)

    def runbatch_path(runno):
        return 'runs%d-%d'%(int(runno/1000)*1000, (int(runno/1000)+1)*1000)

    def run_path(runno):
        return 'by run/%s/run%d'%(runbatch_path(runno), runno)

    def quantity_path(quantity, runno):
        return 'by quantity/%s/%s'%(quantity, runbatch_path(runno))

    def filenames(runno, quantity, extension):
        filename = '/run%d_%s.%s'%(runno, quantity, extension)
        filenames = [ \
            run_path(runno) + filename,
            quantity_path(quantity, runno) + filename ]
        return filenames

    def upload_figure(runno, quantity, f):
        uploader.upload_png_from_figure(filenames(runno, quantity, 'png'), f)

    def upload_figure_dict(runno, fig_dict):
        if(fig_dict):
            for fig_name in fig_dict.keys():
                upload_figure(runno, fig_name, fig_dict[fig_name][0])

    def upload_text(runno, quantity, iostream):
        uploader.upload_from_io(filenames(runno, quantity, 'txt'), 'text/plain', iostream)

    draw_all = True
    draw_all = False if(opt.draw_psd()) else draw_all
    draw_all = False if(opt.draw_high_low()) else draw_all
    draw_all = False if(opt.draw_charge()) else draw_all
    draw_all = False if(opt.draw_missing_components()) else draw_all
    draw_all = False if(opt.draw_pedestal()) else draw_all
    draw_all = False if(opt.draw_temperature()) else draw_all
    draw_all = False if(opt.draw_clock()) else draw_all
    draw_all = False if(opt.draw_data_ordering()) else draw_all
    draw_all = False if(opt.draw_hvpa()) else draw_all
    draw_all = False if(opt.draw_trigger()) else draw_all
    draw_all = False if(opt.draw_waveform_mean()) else draw_all
    draw_all = False if(opt.draw_event()) else draw_all

    ############################################################################
    # FIGURE : power spectra
    ############################################################################

    if(opt.draw_psd() or draw_all):
        fig_dict = calin.diagnostics.stage1_plotting.draw_psd(stage1,
            dataset='pedestal', figure_factory=figure_factory)
        upload_figure_dict(runno, fig_dict)

        fig_dict = calin.diagnostics.stage1_plotting.draw_psd(stage1,
            dataset='pedestal', low_gain=True, figure_factory=figure_factory)
        upload_figure_dict(runno, fig_dict)

        fig_dict = calin.diagnostics.stage1_plotting.draw_psd(stage1,
            dataset='all', figure_factory=figure_factory)
        upload_figure_dict(runno, fig_dict)

        fig_dict = calin.diagnostics.stage1_plotting.draw_psd(stage1,
            dataset='all', low_gain=True, figure_factory=figure_factory)
        upload_figure_dict(runno, fig_dict)

    ############################################################################
    # FIGURE : high-gain vs low-gain values
    ############################################################################

    if(opt.draw_high_low() or draw_all):
        fig_dict, _, _ = calin.diagnostics.stage1_plotting.draw_high_gain_low_gain(stage1,
            dataset='max_sample', subtract_pedestal=True, figure_factory=figure_factory)
        upload_figure_dict(runno, fig_dict)

        fig_dict, _, _ = calin.diagnostics.stage1_plotting.draw_high_gain_low_gain(stage1,
            dataset='opt_sum', subtract_pedestal=True, figure_factory=figure_factory)
        upload_figure_dict(runno, fig_dict)

    ############################################################################
    # FIGURE : charge histograms for all trigger types
    ############################################################################

    if(opt.draw_charge() or draw_all):
        for trigger_type in [ 'physics', 'pedestal', 'external_flasher', 'internal_flasher' ]:
            for low_gain in [ False, True ]:
                fig_dict = calin.diagnostics.stage1_plotting.draw_charge_spectrum(stage1,
                    dataset=trigger_type, low_gain = low_gain,
                    draw_median=True, draw_scale=True, figure_factory = figure_factory)
                upload_figure_dict(runno, fig_dict)

    ############################################################################
    # FIGURE : Missing components
    ############################################################################

    if(opt.draw_missing_components() or draw_all):
        fig_dict = calin.diagnostics.stage1_plotting.draw_missing_components_fraction(stage1,
            figure_factory = figure_factory, mod_label_fontsize=4, aux_label_fontsize=5.5)
        upload_figure_dict(runno, fig_dict)

    ############################################################################
    # FIGURE : Pedestals
    ############################################################################

    if(opt.draw_pedestal() or draw_all):
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

    if(opt.draw_temperature() or draw_all):
        if(stage1.has_nectarcam() and stage1.const_nectarcam().has_ancillary_data() and
                len(stage1.const_nectarcam().const_ancillary_data().feb_temperature_keys())>0):
            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_nectarcam_feb_temperatures(stage1,1,axis=ax)
            ax.set_title('Temperature (FEB 1), run : %d'%runno)
            upload_figure(runno, 'temperature_feb1', ax.figure)

            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_nectarcam_feb_temperatures(stage1,2,axis=ax)
            ax.set_title('Temperature (FEB 2), run : %d'%runno)
            upload_figure(runno, 'temperature_feb2', ax.figure)

            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_nectarcam_feb_temperatures_minmax(stage1,1,axis=ax)
            ax.set_title('Temperature spread (FEB 1), run : %d'%runno)
            upload_figure(runno, 'temperature_spread_feb1', ax.figure)

            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_nectarcam_feb_temperatures_minmax(stage1,2,axis=ax)
            ax.set_title('Temperature spread (FEB 2), run : %d'%runno)
            upload_figure(runno, 'temperature_spread_feb2', ax.figure)

    ############################################################################
    # FIGURE : clock regression
    ############################################################################

    if(opt.draw_clock() or draw_all):
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

    if(opt.draw_data_ordering() or draw_all):
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

    if(opt.draw_hvpa() or draw_all):
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

    if(opt.draw_trigger() or draw_all):
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

    if(opt.draw_waveform_mean() or draw_all):
        if(stage1.has_mean_wf_pedestal()):
            if(stage1.const_mean_wf_pedestal().channel_high_gain_size() \
                    and stage1.const_mean_wf_pedestal().has_camera_high_gain() \
                    and numpy.count_nonzero([ stage1.const_mean_wf_pedestal().channel_high_gain(i).num_entries() \
                        for i in range(stage1.const_mean_wf_pedestal().channel_high_gain_size()) ])):
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
                    and numpy.count_nonzero([ stage1.const_mean_wf_pedestal().channel_low_gain(i).num_entries() \
                        for i in range(stage1.const_mean_wf_pedestal().channel_low_gain_size()) ])):
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

        if(stage1.has_mean_wf_physics()):
            if(stage1.const_mean_wf_physics().channel_high_gain_size() \
                    and stage1.const_mean_wf_physics().has_camera_high_gain() \
                    and numpy.count_nonzero([ stage1.const_mean_wf_physics().channel_high_gain(i).num_entries() \
                        for i in range(stage1.const_mean_wf_physics().channel_high_gain_size()) ])):
                ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
                calin.diagnostics.stage1_plotting.draw_mean_wf(stage1,axis=ax,dataset='physics')
                ax.set_title('Mean waveform (physics, high-gain), run : %d'%runno)
                upload_figure(runno, 'waveform_mean_physics_hg', ax.figure)
            if(stage1.const_mean_wf_physics().channel_low_gain_size() \
                    and stage1.const_mean_wf_physics().has_camera_low_gain() \
                    and numpy.count_nonzero([ stage1.const_mean_wf_physics().channel_low_gain(i).num_entries() \
                        for i in range(stage1.const_mean_wf_physics().channel_low_gain_size()) ])):
                ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
                calin.diagnostics.stage1_plotting.draw_mean_wf(stage1,axis=ax,dataset='physics',low_gain=True)
                ax.set_title('Mean waveform (physics, low-gain), run : %d'%runno)
                upload_figure(runno, 'waveform_mean_physics_lg', ax.figure)

        if(stage1.has_mean_wf_external_flasher()):
            if(stage1.const_mean_wf_external_flasher().channel_high_gain_size() \
                    and stage1.const_mean_wf_external_flasher().has_camera_high_gain() \
                    and numpy.count_nonzero([ stage1.const_mean_wf_external_flasher().channel_high_gain(i).num_entries() \
                        for i in range(stage1.const_mean_wf_external_flasher().channel_high_gain_size()) ])):
                ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
                calin.diagnostics.stage1_plotting.draw_mean_wf(stage1,axis=ax,dataset='external_flasher')
                ax.set_title('Mean waveform (ext flasher, high-gain), run : %d'%runno)
                upload_figure(runno, 'waveform_mean_external_flasher_hg', ax.figure)
            if(stage1.const_mean_wf_external_flasher().channel_low_gain_size() \
                    and stage1.const_mean_wf_external_flasher().has_camera_low_gain() \
                    and numpy.count_nonzero([ stage1.const_mean_wf_external_flasher().channel_low_gain(i).num_entries() \
                        for i in range(stage1.const_mean_wf_external_flasher().channel_low_gain_size()) ])):
                ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
                calin.diagnostics.stage1_plotting.draw_mean_wf(stage1,axis=ax,dataset='external_flasher',low_gain=True)
                ax.set_title('Mean waveform (ext flasher, low-gain), run : %d'%runno)
                upload_figure(runno, 'waveform_mean_external_flasher_lg', ax.figure)

        if(stage1.has_mean_wf_internal_flasher()):
            if(stage1.const_mean_wf_internal_flasher().channel_high_gain_size() \
                    and stage1.const_mean_wf_internal_flasher().has_camera_high_gain() \
                    and numpy.count_nonzero([ stage1.const_mean_wf_internal_flasher().channel_high_gain(i).num_entries() \
                        for i in range(stage1.const_mean_wf_internal_flasher().channel_high_gain_size()) ])):
                ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
                calin.diagnostics.stage1_plotting.draw_mean_wf(stage1,axis=ax,dataset='internal_flasher')
                ax.set_title('Mean waveform (int flasher, high-gain), run : %d'%runno)
                upload_figure(runno, 'waveform_mean_internal_flasher_hg', ax.figure)
            if(stage1.const_mean_wf_internal_flasher().channel_low_gain_size() \
                    and stage1.const_mean_wf_internal_flasher().has_camera_low_gain() \
                    and numpy.count_nonzero([ stage1.const_mean_wf_internal_flasher().channel_high_gain(i).num_entries() \
                        for i in range(stage1.const_mean_wf_internal_flasher().channel_high_gain_size_size()) ])):
                ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
                calin.diagnostics.stage1_plotting.draw_mean_wf(stage1,axis=ax,dataset='internal_flasher',low_gain=True)
                ax.set_title('Mean waveform (int flasher, low-gain), run : %d'%runno)
                upload_figure(runno, 'waveform_mean_internal_flasher_lg', ax.figure)

    ############################################################################
    # FIGURE : mean event rate and on-disk fraction
    ############################################################################

    if(opt.draw_event() or draw_all):
        if(stage1.const_run_info().const_elapsed_time_histogram().sum_w()):
            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_elapsed_time_hist(stage1, axis=ax)
            ax.set_title('Event rate on disk, run : %d'%runno)
            upload_figure(runno, 'event_rate', ax.figure)

        if(stage1.const_run_info().const_log10_delta_t_histogram_all_recorded().sum_w()):
            f = matplotlib.figure.Figure(dpi=figure_dpi)
            f.subplots_adjust(top=0.85)
            ax = f.subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_log_delta_t_histogram(stage1, event_set = 'all', axis=ax)
            ax.set_title('Event $\Delta$T distribution (all events), run : %d'%runno)
            upload_figure(runno, 'event_delta_t_all', ax.figure)

        if(stage1.const_run_info().const_log10_delta_t_histogram().sum_w()):
            f = matplotlib.figure.Figure(dpi=figure_dpi)
            f.subplots_adjust(top=0.85)
            ax = f.subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_log_delta_t_histogram(stage1, event_set = 'consecutive', axis=ax)
            ax.set_title('Event $\Delta$T distribution (consecutive events), run : %d'%runno)
            upload_figure(runno, 'event_delta_t_consecutive', ax.figure)

        if(stage1.const_run_info().const_log10_delta_t_histogram_trigger_physics().sum_w()):
            f = matplotlib.figure.Figure(dpi=figure_dpi)
            f.subplots_adjust(top=0.85)
            ax = f.subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_log_delta_t_histogram(stage1, event_set = 'physics', axis=ax)
            ax.set_title('Event $\Delta$T distribution (consecutive physics events), run : %d'%runno)
            upload_figure(runno, 'event_delta_t_physics', ax.figure)

        if(stage1.const_run_info().const_event_number_histogram().sum_w()):
            ax = matplotlib.figure.Figure(dpi=figure_dpi).subplots(1,1)
            calin.diagnostics.stage1_plotting.draw_event_number_histogram(stage1, axis=ax)
            ax.set_title('Event fraction on disk, run : %d'%runno)
            upload_figure(runno, 'event_fraction', ax.figure)

    ############################################################################
    # PROVENANCE LOG
    ############################################################################

    writer = io.StringIO()
    calin.provenance.printer.print_provenance(writer, stage1.const_provenance_anthology())
    upload_text(runno, 'provenance_log_stage1', writer)

    ############################################################################
    # WRITE SUMMARY SHEET
    ############################################################################

    if(opt.summary_sheet()):
        summary_elements = calin.diagnostics.stage1_summary.stage1_summary_elements(stage1,
            logsheet.get(stage1.run_number(), ''), uploader.get_url(run_path(runno)))
        uploader.append_row_to_sheet(opt.summary_sheet(), summary_elements, row_start=3)

    ############################################################################
    # THE END
    ############################################################################

    print('Finished run :', runno)
    return True

uploader = new_uploader()

logsheet = dict()
if(opt.run_log_sheet()):
    db_rows = uploader.retrieve_sheet(opt.run_log_sheet(),row_start=1)
    logsheet = calin.diagnostics.stage1_summary.make_logsheet_dict(db_rows)

all_oid = []
all_runno = []
all_filename = []
all_status = []

if(opt.summary_sheet() and opt.skip_existing()):
    superset_all_oid, superset_all_runno, superset_all_filenames = get_oids()
    summary_rows = uploader.retrieve_sheet(opt.summary_sheet(),row_start=3)
    for oid, runno, filename in zip(superset_all_oid, superset_all_runno, superset_all_filenames):
        skip_oid = False
        for row in summary_rows:
            if(filename.endswith(row[-1])):
                skip_oid = True
                break
        if(skip_oid):
            print("Skipping :",filename)
        else:
            all_oid.append(oid)
            all_runno.append(runno)
            all_filename.append(filename)
else:
    all_oid, all_runno, all_filename = get_oids()

del uploader

if(not all_oid):
    print("No runs to process, exiting")
    sys.exit()

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


if(opt.upload_to_google_drive() and opt.summary_sheet()):
    uploader = new_uploader()
    uploader.sort_sheet(opt.summary_sheet(), 33, row_start=3)
    del uploader

print("=================================== RESULTS ===================================")
for runno_status in zip(all_runno, all_status):
    print(runno_status[0], "success" if runno_status[1] else "*** FAILED ***")

# The end
