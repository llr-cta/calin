# calin/python/diagnostics/stage1_summary.py -- Stephen Fegan -- 2021-09-02
#
# Generate stage1 summary line elements
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
import numpy
import calin.diagnostics.stage1
import calin.diagnostics.stage1_analysis
import calin.iact_data.instrument_layout

def zero_suppress(x, fmt=None, threshold=0):
    if x is None or numpy.isnan(x):
        return ''
    x_str = fmt.format(x) if fmt else str(x)
    return x_str if x>threshold else ''

def make_logsheet_dict(logsheet_db_rows):
    logsheet = dict()
    for irow, row in enumerate(logsheet_db_rows):
        try:
            run_lo = int(row[0])
            run_hi = int(row[1])
        except:
            run_lo = -1
            run_hi = -1
        url = row[3]
        if(run_lo > 0):
            for run in range(run_lo, run_hi + 1):
                if(run in logsheet):
                    print("Duplicate run : ",run,logsheet[run],row,file=sys.stderr)
                else:
                    logsheet[run] = url
    return logsheet

def stage1_summary_elements(stage1, logsheet_url='', dqm_url=''):
    run_info = stage1.const_run_info()
    run_config = stage1.const_run_config()
    charge_stats = stage1.const_charge_stats()

    camera_layout = run_config.const_camera_layout()

    run_duration = calin.diagnostics.stage1_analysis.run_duration(stage1)
    num_event_missing = int(run_info.event_numbers_found().end_index()[-1])-run_info.num_events_found()-1+int(sum(run_info.const_duplicate_event_numbers().count()))

    nsb_convert_60 = 58**2*1.19*(60-2.0308)
    nsb_convert_16 = 58**2*1.19*(16-2.0308)
    elec_var_60 = 0
    elec_var_16 = 0

    elements = [
        '%d'%run_config.run_number(),
        '%d'%run_info.num_events_found(),
        zero_suppress(run_info.num_mono_trigger()),
        zero_suppress(run_info.num_pedestal_trigger()),
        zero_suppress(run_info.num_external_calibration_trigger()),
        zero_suppress(run_info.num_internal_calibration_trigger()),

        zero_suppress(int(numpy.sum(run_info.const_duplicate_event_numbers().count()))),
        zero_suppress(num_event_missing),
        zero_suppress(run_info.num_events_missing_cdts()),
        zero_suppress(run_info.num_events_missing_tib()),
        zero_suppress(run_info.num_events_missing_tib_and_cdts()),
        zero_suppress(run_info.num_events_missing_modules()),

        '%.3f'%run_duration if run_duration is not None else '',
        '%.3f'%((run_info.num_events_found()+num_event_missing)/run_duration) if run_duration is not None else '',
        '%.3f'%(run_info.num_events_found()/run_duration) if run_duration is not None else '',
        '%.3f'%numpy.std(run_info.const_elapsed_time_histogram().bins()[1:-1]) if run_duration is not None else '',
        zero_suppress(run_info.num_mono_trigger()/run_duration, fmt='{:.3f}') if run_duration is not None else '',
        zero_suppress(run_info.num_pedestal_trigger()/run_duration, fmt='{:.3f}') if run_duration is not None else '',

        '%.3f'%((numpy.mean(charge_stats.const_high_gain().ped_trigger_full_wf_var())-elec_var_60)/nsb_convert_60) if run_info.num_pedestal_trigger() else '',
        '%.3f'%((numpy.mean(charge_stats.const_high_gain().all_trigger_ped_win_var())-elec_var_16)/nsb_convert_16) if run_info.num_external_calibration_trigger()==0 else '',
        zero_suppress(calin.diagnostics.stage1_analysis.median_voltage(stage1),fmt='{:.1f}',threshold=-1),
        zero_suppress(calin.diagnostics.stage1_analysis.median_feb_temp(stage1),fmt='{:.2f}',threshold=-40),
        zero_suppress(calin.diagnostics.stage1_analysis.spread_feb_temp(stage1),fmt='{:.2f}',threshold=-1),
        '%.2f'%(charge_stats.const_low_gain().ext_trigger_all_channel_opt_win_mean()/run_config.configured_channel_id_size()-250*16) if charge_stats.const_low_gain().ext_trigger_all_channel_count() else '',
        '%.2f'%(numpy.sqrt(charge_stats.const_low_gain().ext_trigger_all_channel_opt_win_var())/run_config.configured_channel_id_size()) if charge_stats.const_low_gain().ext_trigger_all_channel_count() else '',

        '=HYPERLINK("' + logsheet_url + '","Logbook")' if logsheet_url else '',
        '=HYPERLINK("' + dqm_url + '","DQM")' if dqm_url else '',
        '%d'%run_config.fragment_filename_size(),
        '%.3f'%(run_config.file_size()/1e9),
        '%d'%run_config.configured_module_id_size(),
        '%.1f'%calin.diagnostics.stage1_analysis.num_wf(stage1),
        '%d'%run_config.num_samples(),
        '%.1f'%(run_config.file_size()/run_info.num_events_found()/1e3) if run_info.num_events_found()>10000 else '',

        '/'.join(run_config.filename().split('/')[-3:])\

    ]

    return elements
