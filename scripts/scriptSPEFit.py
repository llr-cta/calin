#!/usr/bin/env python3

import numpy as np
import calin.calib.spe_fit
import calin.calib.pmt_ses_models
import calin.math.data_modeling
import calin.math.optimizer
import calin.iact_data.llr
import calin.SPEFit

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import calin.plotting

import csv
import pickle
import sys

def runScript():
    DataListSignal = []
    DataListPedestal = []
    #Filename = '_OD1_18_16_dark_opt_ASCII_histo'
    
    # Filename is a launch argument
    Filename = sys.argv[1]
    Path = sys.argv[2]

    #DataListSignal.append(calin.iact_data.llr.make_lecroy_adc_hist(sys.argv[1],scale=1))
    # ~ DataListSignal = calin.iact_data.llr.make_lecroy_adc_hist('/CTA/DataTestBench/PierreJeanTestBench/spectrum_spe_ch0_1540V.txt',scale=1)
    # ~ DataListSignal = calin.iact_data.llr.make_lecroy_adc_hist('/CTA/SPE_fitter/data/_OD0_flash7.8V_LED7_0.3mA_ASCII_histo.txt', scale=1)
    
    # Loading of the file and contruction of the calin histogram called DataListSignal
    # here delimiter are , but can be change to other stuffs, skiprows is the number of column that we have to skip at the begin of the ASCII file
    scale=1
    data = np.loadtxt(Path+Filename+'.txt', delimiter=',', skiprows=5)
    dx = np.median(data[1:,0] - data[0:-1,0])/scale
    DataListSignal = calin.math.histogram.SimpleHist(dx)
    DataListSignal.insert_two_vec(data[:,0]/scale, data[:,1])
    

    SimpleModel = False

    if (not SimpleModel):
        # 2 Gauss Fit with constrained model, fast mode
        # Fit_2_gauss(hsignal), if you specify Fit_2_gauss(hsignal,hped), than 
        # On the following I put hped = DataListSignal, because hped is not used
        # output : result is a table with parameters
        # aerror is a table with their uncertainty
        
        result,aerror,mes_2g,ses_2g,pedestal_2g,param_2g,delta_param_2g,gain_2g,delta_gain_2g = calin.SPEFit.Fit_2_gauss(DataListSignal)
        print(param_2g)
        print(delta_param_2g)
        
        # 2 Gauss Fit with constrained model, Robust mode, all parameters free, should only work at high voltage (HV > 1200)
        # Use it at your own risk
        # ~ result,aerror,b,c,d = calin.SPEFit.Fit_2_gauss(DataListSignal,DataListSignal,1000,False,True,True)
    else:
        result,aerror,mes_2g,c,d = calin.SPEFit.Fit_1_gauss(DataListSignal)

    print(aerror)

    # Compute the mean of the single electron spectrum, this is the gain
    #Compute also the variance and the resolution
    if (not SimpleModel):
        ses_mean = gain_2g
        ses_mean_uncertainty = delta_gain_2g
        ses_var = sum(mes_2g.all_ses_x()**2 * mes_2g.single_electron_spectrum())*DataListSignal.dxval()-ses_mean**2
        ses_res = np.sqrt(ses_var)/ses_mean
    else:
        ses_mean = sum(mes_2g.all_ses_x() * mes_2g.single_electron_spectrum())*DataListSignal.dxval()
        ses_mean_uncertainty = aerror[3]
        ses_var = sum(mes_2g.all_ses_x()**2 * mes_2g.single_electron_spectrum())*DataListSignal.dxval()-ses_mean**2
        ses_res = np.sqrt(ses_var)/ses_mean
    
    # plot the result on a pdf called "myfig"
    
    calin.plotting.plot_histogram(DataListSignal, lw=2, label='SPE data')
    plt.xlabel('Signal [DC]')
    plt.ylabel('Events per %d DC bin'%DataListSignal.dxval())

    ihist = range(0,DataListSignal.nbin())
    xhist = DataListSignal.all_xval_center()

    mes_2g.set_parameter_values(result)
    ymodel_2g = list(map(lambda x: DataListSignal.sum_w()*DataListSignal.dxval()*mes_2g.pdf_mes(x),xhist))
    matplotlib.pyplot.plot(xhist,ymodel_2g, 'r', label='Gaussian \n Gain = '+str(round(ses_mean,2))+' +/- ' + str(round(ses_mean_uncertainty,2))+'\n res = '+str(round(ses_res,3))+' +/- ' + str(round(delta_param_2g[2],3)))
    a=list(plt.axis())
    plt.axis(a)
    plt.legend()
    plt.grid()
    plt.savefig('/CTA/myfig'+Filename)

    F = open("/CTA/"+Filename+"Output.txt","w")
    print("SPE: spe_mean = %-6.2f ; spe_var = %-6.2f; spe_res = %-6.3f; "%(ses_mean,np.sqrt(ses_var),np.sqrt(ses_var)/ses_mean))
    for i in range(len(a)):
        F.write(str(result[i])+'\t')
        F.write(str(aerror[i])+'\n')
    F.write(str(ses_mean)+'\t')
    F.write(str(ses_mean_uncertainty)+'\n')
    F.close()

runScript()
