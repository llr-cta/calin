import numpy as np
import calin.calib.spe_fit
import calin.calib.pmt_ses_models
import calin.math.data_modeling
import calin.math.optimizer
import calin.iact_data.llr

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import calin.plotting
from scipy.misc import derivative

import csv
import pickle
import sys

PedestalGuess = 0

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def PMax(r):
    if (np.pi*r**2/(np.pi*r**2 + np.pi - 2*r**2 - 2) <= 1):
        return np.pi*r**2/(np.pi*r**2 + np.pi - 2*r**2 - 2)
    else:
        return 1

def ax(p,res):
    return ((2/np.pi)*p**2-p/(res**2+1))

def bx(p,mu2):
    return (np.sqrt(2/np.pi)*2*p*(1-p)*mu2)

def cx(sig2,mu2,res,p):
    return (1-p)**2*mu2**2 - (1-p)*(sig2**2+mu2**2)/(res**2+1)

def delta(p,res,sig2,mu2):
    return bx(p,mu2)*bx(p,mu2) - 4*ax(p,res)*cx(sig2,mu2,res,p)

def ParamU(p,r):
    return ((8*(1-p)**2*p**2)/np.pi - 4*(2*p**2/np.pi - p/(r**2+1))*((1-p)**2-(1-p)/(r**2+1)))

def ParamS(p,r):
    return (4*(2*p**2/np.pi - p/(r**2+1))*(1-p))/(r**2+1)

def SigMin(p,res,mu2):
    return mu2*np.sqrt((-ParamU(p,res)+(bx(p,mu2)**2/mu2**2))/(ParamS(p,res)))

def SigMax(p,res,mu2):
    return mu2*np.sqrt((-ParamU(p,res))/(ParamS(p,res)))

def sigma1(p,res,sig2,mu2):
    return (-bx(p,mu2)+np.sqrt(delta(p,res,sig2,mu2)))/(2*ax(p,res))

def sigma2(n,p,res,mu2):
    if ((-ParamU(p,res)+(bx(p,mu2)**2/mu2**2))/(ParamS(p,res)) > 0):
        return SigMin(p,res,mu2)+n*(SigMax(p,res,mu2)-SigMin(p,res,mu2))
    else:
        return n*SigMax(p,res,mu2)

def Gain(pp,mu1,res,mu2,n):
    p = pp*PMax(res)
    sig2 = sigma2(n,p,res,mu2)
    return (1-p)*mu2 + 2*p*sigma1(p,res,sig2,mu2)/np.sqrt(2*np.pi)

def partial_derivative(func, var=0, point=[]):
    args = point[:]
    def wraps(x):
        args[var] = x
        return func(*args)
    return derivative(wraps, point[var], dx = 1e-6)

def GainUncertainty(pp,mu1,res,mu2,n,CovMatrix):
    uncertainty = 0
    derivative = []
    for i in range(5):
        derivative.append(partial_derivative(Gain, i, [pp,mu1,res,mu2,n]))
    for i in range(5):
        for j in range(5):
            uncertainty = uncertainty + derivative[i]*derivative[j]*CovMatrix[i,j]
    return np.sqrt(uncertainty)

def Fit_1_gauss(h,hped = 0,UseHped=False):
    if (not UseHped):
        hped = h
	
    import time
    histTable = np.zeros(h.size())
    histCharge = np.zeros(h.size())
    histTablePedestal = np.zeros(hped.size())
    histChargePedestal = np.zeros(hped.size())
    global MinCharge
    global NSample
    
    for i in range(h.size()):
        histTable[i] = h.weight(i)
        histCharge[i] = h.xval_center(i)
    for i in range(hped.size()):
        histTablePedestal[i] = hped.weight(i)
        histChargePedestal[i] = hped.xval_center(i)
    # ~ print(h.dxval())
    # ~ print(hped.dxval())
    MinCharge = min(histCharge[0],histChargePedestal[0])
    #print(MinCharge)
    MaxCharge = max(histCharge[len(histCharge)-1],histChargePedestal[len(histChargePedestal)-1])
    #print(MaxCharge)
    if (h.dxval() != hped.dxval()):
        print("hped and hsignal does not have the same dx !")
    NSample = int((MaxCharge-MinCharge)/h.dxval())+1
    #print(NSample)
    MaxIndex = np.argmax(histTable)   
    PedestalGuess = np.sum(histCharge[MaxIndex-10:MaxIndex+10]*histTable[MaxIndex-10:MaxIndex+10])/np.sum(histTable[MaxIndex-10:MaxIndex+10])
    # ~ print(histTable)
    # ~ print("PedestalGuess ")
    # ~ print(PedestalGuess)
    diffTable = np.diff(histTable)
    x = np.linspace(int(-h.size()/2.), int(h.size()/2.), h.size())
    Fporte = np.where(abs(x)<=10, 1, 0)
    diffTable = np.convolve(diffTable, Fporte,'same')
    
    # ~ print(diffTable[MaxIndex:])
    PositiveTable = 0
    index = 0
    while (PositiveTable == 0):
        PositiveTableTemp = np.where( diffTable[MaxIndex+index:] > 0 )
        if PositiveTableTemp:
            PositiveTable = PositiveTableTemp[0][0]
            index = index+1
        else:
            PositiveTable = len(histCharge)-1
			
    # ~ print(PositiveTable)
    # ~ print("here")
    # ~ print(histTable[:MaxIndex+PositiveTable])
    MeanPedestalGuess = np.sum(histCharge[:MaxIndex+PositiveTable]*histTable[:MaxIndex+PositiveTable])/np.sum(histTable[:MaxIndex+PositiveTable])
    PedestalRMSGuess = np.sum(np.square(histCharge[:MaxIndex+PositiveTable])*histTable[:MaxIndex+PositiveTable])/np.sum(histTable[:MaxIndex+PositiveTable]) - MeanPedestalGuess**2
    PedestalRMSGuess = np.sqrt(PedestalRMSGuess)
    # ~ print(PedestalRMSGuess)
    # ~ print("Gain guess")
    # ~ print(histCharge[np.argmax(histTable[MaxIndex+PositiveTable:])+MaxIndex+PositiveTable]-histCharge[MaxIndex])
    GainGuess = histCharge[np.argmax(histTable[MaxIndex+PositiveTable:])+MaxIndex+PositiveTable]-histCharge[MaxIndex]
    GainGuess2 = (np.sum(histCharge[MaxIndex+PositiveTable:]*histTable[MaxIndex+PositiveTable:])/np.sum(histTable[MaxIndex+PositiveTable:])) - histCharge[MaxIndex]
    if (GainGuess < 0.5*GainGuess2):
        GainGuess = GainGuess2
    MaxGain = 0.2*NSample*h.dxval()
    if (GainGuess > MaxGain):
        GainGuess = 0.13*NSample*h.dxval()
    LuminosityGuess = histTable[np.argmax(histTable[MaxIndex+PositiveTable:])+MaxIndex+PositiveTable]*GainGuess*0.3/(histTable[MaxIndex]*PedestalRMSGuess)
    # ~ print(LuminosityGuess)
    # ~ print(histTable[np.argmax(histTable[MaxIndex+PositiveTable:])+MaxIndex+PositiveTable])
    # ~ print(histTable[MaxIndex])
    # ~ if (GainGuess < MeanPedestalGuess and GainGuess < PedestalGuess):
		# ~ GainGuess = 
    #PedestalRMSGuess = diffTable
    #print("PedestalRMSGuess ")
    #print(PedestalRMSGuess)
    
    ftime0 = time.time()
    K = 6.68060614323e-13
    res = 0.3
    if (LuminosityGuess > 1.5):
        LuminosityGuess = 1
    LuminosityGuess = 1
    iv = [LuminosityGuess, MeanPedestalGuess, PedestalRMSGuess, GainGuess, res*GainGuess]
    FixedParam = [-1]
    a,aerror,b,c,d,e,f,g = Fit_1_gauss_bis(h,hped,iv,True,False,60,FixedParam)
    return (a,aerror,b,f,g)
    
def Fit_2_gauss(h,hped = 0,UseHped=False,RobustMode=False,FreeMode = False, LightFixed = -1, startParameters = [], UseStartParameters = False, UseSplineMode = False):
	
    if (not UseHped):
        hped = h
	
    import time
    #global PedestalGuess
    histTable = np.zeros(h.size())
    histCharge = np.zeros(h.size())
    histTablePedestal = np.zeros(hped.size())
    histChargePedestal = np.zeros(hped.size())
    global MinCharge
    global NSample
    
    for i in range(h.size()):
        histTable[i] = h.weight(i)
        histCharge[i] = h.xval_center(i)
    for i in range(hped.size()):
        histTablePedestal[i] = hped.weight(i)
        histChargePedestal[i] = hped.xval_center(i)
    # ~ print(h.dxval())
    # ~ print(hped.dxval())
    MinCharge = min(histCharge[0],histChargePedestal[0])
    #print(MinCharge)
    MaxCharge = max(histCharge[len(histCharge)-1],histChargePedestal[len(histChargePedestal)-1])
    #print(MaxCharge)
    if (h.dxval() != hped.dxval()):
        print("hped and hsignal does not have the same dx !")
    NSample = int((MaxCharge-MinCharge)/h.dxval())+1
    #print(NSample)
    MaxIndex = np.argmax(histTable)   
    PedestalGuess = np.sum(histCharge[MaxIndex-10:MaxIndex+10]*histTable[MaxIndex-10:MaxIndex+10])/np.sum(histTable[MaxIndex-10:MaxIndex+10])
    # ~ print(histTable)
    # ~ print("PedestalGuess ")
    # ~ print(PedestalGuess)
    diffTable = np.diff(histTable)
    x = np.linspace(int(-h.size()/2.), int(h.size()/2.), h.size())
    Fporte = np.where(abs(x)<=10, 1, 0)
    diffTable = np.convolve(diffTable, Fporte,'same')
    # ~ print(diffTable[MaxIndex:])
    PositiveTable = 0
    index = 0
    while (PositiveTable == 0):
        PositiveTableTemp = np.where( diffTable[MaxIndex+index:] > 0 )
        if PositiveTableTemp:
            PositiveTable = PositiveTableTemp[0][0]
            index = index+1
        else:
            PositiveTable = len(histCharge)-1
			
    # ~ print(PositiveTable)
    # ~ print("here")
    # ~ print(histTable[:MaxIndex+PositiveTable])
    MeanPedestalGuess = np.sum(histCharge[:MaxIndex+PositiveTable]*histTable[:MaxIndex+PositiveTable])/np.sum(histTable[:MaxIndex+PositiveTable])
    PedestalRMSGuess = np.sum(np.square(histCharge[:MaxIndex+PositiveTable])*histTable[:MaxIndex+PositiveTable])/np.sum(histTable[:MaxIndex+PositiveTable]) - MeanPedestalGuess**2
    PedestalRMSGuess = np.sqrt(PedestalRMSGuess)
    # ~ print(PedestalRMSGuess)
    # ~ print("Gain guess")
    # ~ print(histCharge[np.argmax(histTable[MaxIndex+PositiveTable:])+MaxIndex+PositiveTable]-histCharge[MaxIndex])
    GainGuess = histCharge[np.argmax(histTable[MaxIndex+PositiveTable:])+MaxIndex+PositiveTable]-histCharge[MaxIndex]
    GainGuess2 = (np.sum(histCharge[MaxIndex+PositiveTable:]*histTable[MaxIndex+PositiveTable:])/np.sum(histTable[MaxIndex+PositiveTable:])) - histCharge[MaxIndex]
    if (GainGuess < 0.5*GainGuess2):
        GainGuess = GainGuess2
    MaxGain = 0.2*NSample*h.dxval()
    if (GainGuess > MaxGain):
        GainGuess = 0.13*NSample*h.dxval()
    LuminosityGuess = histTable[np.argmax(histTable[MaxIndex+PositiveTable:])+MaxIndex+PositiveTable]*GainGuess*0.3/(histTable[MaxIndex]*PedestalRMSGuess)
    # ~ print(LuminosityGuess)
    # ~ print(histTable[np.argmax(histTable[MaxIndex+PositiveTable:])+MaxIndex+PositiveTable])
    # ~ print(histTable[MaxIndex])
    # ~ if (GainGuess < MeanPedestalGuess and GainGuess < PedestalGuess):
		# ~ GainGuess = 
    #PedestalRMSGuess = diffTable
    #print("PedestalRMSGuess ")
    #print(PedestalRMSGuess)
    
    
    K = 6.68060614323e-13
    res = 0.485
    pp = 0.45
    LuminosityGuess = 1
    luminosity = LuminosityGuess
    PedestalMean = PedestalGuess
    xspline = 0
    #print(K*HV**(4.71080288507))
    nParameterPed = 2
    if RobustMode:
        iv = [luminosity, MeanPedestalGuess, PedestalRMSGuess,pp, res, GainGuess,0.715] 
        if UseStartParameters:
            iv = startParameters
        FixedParam = [3,6]
        if FreeMode:
           FixedParam = [-1]
        if (LightFixed > 0):
            iv[0] = LightFixed
            FixedParam = FixedParam+[0]
    if UseSplineMode:
        nSpline = 5
        nParameterPed = nSpline+1
        if RobustMode:
            iv = [luminosity, MeanPedestalGuess, PedestalRMSGuess,pp, res, GainGuess,0.715] 
            xspline = np.linspace(MeanPedestalGuess-5*PedestalRMSGuess,MeanPedestalGuess+5*PedestalRMSGuess,nSpline)
            spline = [0.05] + gaussian(xspline,MeanPedestalGuess,PedestalRMSGuess).tolist()
            iv = [iv[0]] + spline + iv[3:7]
        else:
            iv = [luminosity, MeanPedestalGuess, PedestalRMSGuess,pp,0, res, GainGuess,0.715]
            xspline = np.linspace(MeanPedestalGuess-5*PedestalRMSGuess,MeanPedestalGuess+5*PedestalRMSGuess,nSpline)
            spline = [0.05] + gaussian(xspline,MeanPedestalGuess,PedestalRMSGuess).tolist()
            iv = [iv[0]] + spline + iv[3:8]
        
        if UseStartParameters:
            iv = startParameters
        if RobustMode:
            FixedParam = [len(spline)+1,len(spline)+4]
        else:
            FixedParam = [len(spline)+1,len(spline)+2,len(spline)+5]
        if FreeMode:
            if RobustMode:
                FixedParam = [-1]
            else:
                FixedParam = [len(spline)+2]
        if (LightFixed > 0):
            iv[0] = LightFixed
            FixedParam = FixedParam+[0]
    if (not UseSplineMode and not RobustMode):    
        iv = [luminosity, MeanPedestalGuess, PedestalRMSGuess,pp,0, res, GainGuess,0.715] #1150V
        if UseStartParameters:
            iv = startParameters
        FixedParam = [3,4,7]
        if FreeMode:
            FixedParam = [4]
        if (LightFixed > 0):
            iv[0] = LightFixed
            FixedParam = FixedParam+[0]
    
    
    #FixedParam = [4]
    #FixedParam = [-1]
    
    end = True
    while(end and RobustMode):
        a,aerror,b,c,d,e,f,g = Fit_2_gauss_bis(h,hped,iv,False,False,2,FixedParam,UseHped,UseSplineMode,xspline)
        n=0
        iv = a
        print("results")
        print(a)
        for j in range(len(iv)):
            #print(a)
            #print(c)
            #print(d)
            #or abs(a[j]-c[j])/abs(a[j])<0.02
            if ((a[j] <= c[j]) and not(j in FixedParam)):
                iv[j] = c[j]
                n += 1
            elif ((a[j] >= d[j]) and not(j in FixedParam)):
                iv[j] = d[j]
                n += 1
            elif (not(j in FixedParam)):
                iv[j] = a[j]
        if (n == 0):
            a,aerror,b,c,d,e,f,g = Fit_2_gauss_bis(h,hped,iv,False,False,5,FixedParam,UseHped,UseSplineMode,xspline)
            iv = a
            for j in range(len(iv)):
                if ((a[j] <= c[j]) and not(j in FixedParam)):
                    iv[j] = c[j]
                    n += 1
                elif ((a[j] >= d[j]) and not(j in FixedParam)):
                    iv[j] = d[j]
                    n += 1
                elif (not(j in FixedParam)):
                    iv[j] = a[j]
            print(n)
        if (n ==0):
            a,aerror,b,c,d,e,f,g = Fit_2_gauss_bis(h,hped,iv,False,False,20,FixedParam,UseHped,UseSplineMode,xspline)
            iv = a
            for j in range(len(iv)):
                if ((a[j] <= c[j]) and not(j in FixedParam)):
                    iv[j] = c[j]
                    n += 1
                elif ((a[j] >= d[j]) and not(j in FixedParam)):
                    iv[j] = d[j]
                    n += 1
                elif (not(j in FixedParam)):
                    iv[j] = a[j]
            print(n)
        
        if (n ==0):
            end = False
            print("Fixed Param")
            print(FixedParam)
            if not UseSplineMode:
                for i in range(len(FixedParam)):
                    if (FixedParam[i] > 4):
                        FixedParam[i] += 1
                FixedParam.append(4)
            else:
                for i in range(len(FixedParam)):
                    if (FixedParam[i] > len(xspline)+3):
                        FixedParam[i] += 1
                FixedParam.append(len(xspline)+3)
                # ~ for i in range(len(xspline)+1):
                    # ~ FixedParam.append(i+1)
            print("Fixed Param")
            print(FixedParam)
            iv[:] = iv[:] #+ 0.0001*iv[:]#[iv[0]+0.0001*iv[0], iv[1]+0.0001*iv[1], iv[2]+0.0001*iv[2],iv[3]+0.0001*iv[3],0, iv[4]+0.0001*iv[4], iv[5]+0.0001*iv[5],iv[6]+0.0001*iv[6]]
            if (UseSplineMode):
                iv = np.concatenate((np.concatenate((iv[:len(xspline)+3],[0]), axis=None), iv[len(xspline)+3:]),axis=None)
            else:
                iv = [iv[0], iv[1], iv[2],iv[3],0, iv[4], iv[5],iv[6]]

    # ~ if (not RobustMode):
    a,aerror,b,c,d,e,f,g = Fit_2_gauss_bis(h,hped,iv,True,False,60,FixedParam,UseHped,UseSplineMode,xspline)
    gain = sum(b.all_ses_x() * b.single_electron_spectrum())*h.dxval()
    #var = sum(b.all_ses_x()**2 * b.single_electron_spectrum())*h.dxval()-gain**2
    CovMatrix = np.zeros((5,5))
    # ~ for i in range(nParameterPed+1):
        # ~ FixedParam = FixedParam + [i+1]
    for i in range(5):
        for j in range(5):
            if ((j+nParameterPed+1 in FixedParam) or (i+nParameterPed+1 in FixedParam)):
                CovMatrix[i,j] = 0
            else:
                ni = 0
                nj = 0
                for k in range(j+nParameterPed+1,0,-1):
                    if (k in FixedParam):
                        nj = nj+1
                for k in range(i+nParameterPed+1,0,-1):
                    if (k in FixedParam):
                        ni = ni+1
                CovMatrix[i,j] = aerror[i+nParameterPed+1-ni,j+nParameterPed+1-nj]
				
    print(CovMatrix)
    resultParam = np.zeros(5)
    for i in range(5):
        if (i+nParameterPed+1 in FixedParam):
            print(i)
            resultParam[i] = iv[i+nParameterPed+1]
        else:
            ni=0
            for k in range(i+nParameterPed+1,0,-1):
                if (k in FixedParam):
                    ni = ni+1
            #if (not RobustMode):
            #    ni = ni-1
            resultParam[i] = a[i+nParameterPed+1-ni]
    DeltaGain = GainUncertainty(resultParam[0],resultParam[1],resultParam[2],resultParam[3],resultParam[4],CovMatrix)
    #DeltaVar = np.sqrt(CovMatrix[2,2])
    UncertParam = np.zeros(5)
    for i in range(5):
        UncertParam[i] = np.sqrt(CovMatrix[i,i])
    
    print("Result Param : ",resultParam)
    print("Gain ",gain," +/- ",DeltaGain)
    
    return (a,aerror,b,f,g,resultParam,UncertParam,gain,DeltaGain)
        
        
def Fit_1_gauss_bis(h,hped,iv,Optimizer=False,Verbose=False,WallTime=2,FixedParam=[-1],UseHped=False):
    nSpline = len(iv)-7
    ped = calin.math.pdf_1d.BinnedGaussianPDF(h.dxval())
    ped.set_parameter_values(np.asarray([0, 50]))
    
    ses_2g = calin.math.pdf_1d.BinnedGaussianPDF(h.dxval())
    ses_2g.set_parameter_values(np.asarray([1,0.5]))
    #mes_2g = calin.calib.spe_fit.GeneralPoissonMES(-200, h.dxval(), 4094, ses_2g, ped) #test benche
    global MinCharge
    global NSample
    mes_2g = calin.calib.spe_fit.GeneralPoissonMES(MinCharge, h.dxval(), NSample, ses_2g, ped) #test bench
    
    if (UseHped):
        cost_2g = calin.calib.spe_fit.SPELikelihood(mes_2g, h, hped)
    else:
        cost_2g = calin.calib.spe_fit.SPELikelihood(mes_2g, h)
    
    if (Optimizer):
        opt_2g = calin.math.optimizer.NLOptOptimizer("LD_LBFGS", cost_2g)
    else:
        opt_2g = calin.math.optimizer.NLOptOptimizer("GD_MLSL", cost_2g)
    opt_2g.set_verbosity_level(calin.math.optimizer.OptimizerVerbosityLevel_ELEVATED);
    opt_2g.set_abs_tolerance(0.00001);
    opt_2g.set_max_walltime(WallTime);
    
    limitLow = np.linspace(0,1,len(iv))
    limitUp = np.linspace(0,1,len(iv))
    Stepsize = np.linspace(0,1,len(iv))
    
    for i in range(0,len(iv)):
        limitLow[i] = iv[i] - np.sign(iv[i])*iv[i]*0.95
        limitUp[i] = iv[i] + np.sign(iv[i])*iv[i]*1
        Stepsize[i] = abs(iv[i])*0.01
        if (i in FixedParam):
            print(i)
            limitLow[i] = iv[i] - np.sign(iv[i])*iv[i]*1e-15
            limitUp[i] = iv[i] + np.sign(iv[i])*iv[i]*1e-15
            Stepsize[i] = abs(iv[i])*1e-15
    
    opt_2g.set_initial_values(iv);
    #opt_2g.set_initial_values([limitLow[0]]+lowb+limitLow[3:6]);
    opt_2g.set_limits_lo(limitLow)
    opt_2g.set_limits_hi(limitUp)
    opt_2g.set_scale(Stepsize)
    #opt_2g.set_limits_lo([0.1, -10, 1, 0.1, 50,   209, 1])
    #opt_2g.set_limits_hi([5,-0.1, 15, 0.5,  150, 211, 200])
    status, xopt_2g, fval_2g = opt_2g.minimize()
    if (Optimizer):
        status, err_mat= opt_2g.calc_error_matrix()
        xerr = np.sqrt(err_mat.diagonal())
        print(xerr)
        print(fval_2g)
        print(status)
    if Verbose:
        calin.plotting.plot_histogram(h, lw=2, label='SPE data')
        xlabel('Signal [DC]')
        ylabel('Events per %d DC bin'%h.dxval())
        ihist = range(0,h.nbin());
        xhist = h.all_xval_center()
        mes_2g.set_parameter_values(xopt_2g)
        ymodel_2g = \
        list(map(lambda x: h.sum_w()*h.dxval()*mes_2g.pdf_mes(x),xhist))
        plot(xhist,ymodel_2g, 'r', label='Two Gaussian')
        a=list(axis())
    return (xopt_2g,xerr,mes_2g,limitLow,limitUp,status,ses_2g,ped)

    
def Fit_2_gauss_bis(h,hped,iv,Optimizer=False,Verbose=False,WallTime=2,FixedParam=[-1],UseHped=False, UseSplineMode=False, Xspline = 0):
    
    print(iv)
    if (UseSplineMode):
        nSpline = len(Xspline)-1
    
    if (not UseSplineMode):
        ped = calin.math.pdf_1d.BinnedGaussianPDF(h.dxval())
        ped.set_parameter_values(np.asarray([0, 50]))
    else:
        ped = calin.math.pdf_1d.LogQuadraticSpline1DPDF(Xspline,Xspline[0],Xspline[len(Xspline)-1])
        if not (0 in FixedParam):
            print("Parameter size is : ", len(iv[1:nSpline+3]))
            print("Parameter are ", iv[1:nSpline+3])
            ped.set_parameter_values(np.asarray(iv[1:nSpline+3]))
        else:
            ped.set_parameter_values(np.asarray(iv[:nSpline+2]))
    
    ses_2g = calin.calib.pmt_ses_models.TwoGaussianSESConstrained(h.dxval())
    ses_2g.set_parameter_values(np.asarray([0.2,0.5,1,0.5]))
    
    
    iv2 = []
    if (UseSplineMode):
        NPed = nSpline+2
    else:
        NPed = 2
    if (Optimizer):
        for i in range(0,len(iv)):
            if not(i in FixedParam and i >= NPed+1):
                iv2 = iv2 + [iv[i]]
            else:
                if (i >= NPed+1):
                    print(i)
                    print(NPed+1)
                    print("Remove parameter :", i-(NPed+1))
                    ses_2g.remove_parameter_from_subspace(i-(NPed+1),iv[i])
    else:
        iv2=iv
    print("iv2")
    print(iv2)
    global MinCharge
    global NSample
    #print(MinCharge)
    #print("max is")
    #print(MinCharge+NSample*h.dxval())
    mes_2g = calin.calib.spe_fit.GeneralPoissonMES(MinCharge, h.dxval(), NSample, ses_2g, ped) #test benche
    
    if (UseHped):
        cost_2g = calin.calib.spe_fit.SPELikelihood(mes_2g, h, hped)
    else:
        cost_2g = calin.calib.spe_fit.SPELikelihood(mes_2g, h)
    #cost_2g = calin.calib.spe_fit.SPELikelihood(mes_2g, h, hped)
    if (Optimizer):
        opt_2g = calin.math.optimizer.NLOptOptimizer("LD_LBFGS", cost_2g)
    else:
        opt_2g = calin.math.optimizer.NLOptOptimizer("GD_MLSL", cost_2g)
    opt_2g.set_verbosity_level(calin.math.optimizer.OptimizerVerbosityLevel_ELEVATED);
    opt_2g.set_abs_tolerance(0.0000001);
    opt_2g.set_max_walltime(WallTime);
    
    limitLow = np.linspace(0,1,len(iv2))
    limitUp = np.linspace(0,1,len(iv2))
    Stepsize = np.linspace(0,1,len(iv2))
    nRes = 0
    nWidth = 0
    for k in range(NPed+3,0,-1):
        if (k in FixedParam):
            nRes = nRes +1
    for k in range(NPed+5,0,-1):
        if (k in FixedParam):
            nWidth = nWidth +1
    if(Optimizer):
        for i in range(0,len(iv2)):
            limitLow[i] = iv2[i] - np.sign(iv2[i])*iv2[i]*(1-1e-15)
            limitUp[i] = iv2[i] + np.sign(iv2[i])*iv2[i]*100
            Stepsize[i] = abs(iv2[i])*0.0001
            if (i == NPed+3-nRes and not (5 in FixedParam)):
                print("Resolution param : ", iv2[i])
                limitLow[i] = iv2[i] - np.sign(iv2[i])*iv2[i]*0.20
                limitUp[i] = iv2[i] + np.sign(iv2[i])*iv2[i]*0.20
                Stepsize[i] = abs(iv2[i])*0.0001
            if (i == NPed+5-nWidth and not (7 in FixedParam)):
                print("Width param : ", iv2[i])
                if (iv2[i] > 0.85):
                    iv2[i] = 0.8
                limitUp[i] = 0.85
                Stepsize[i] = abs(iv2[i])*0.0001
            if ((i in FixedParam) and i < NPed+1):
                limitLow[i] = iv2[i] - np.sign(iv2[i])*iv2[i]*1e-15
                limitUp[i] = iv2[i] + np.sign(iv2[i])*iv2[i]*1e-15
                Stepsize[i] = abs(iv2[i])*1e-15
            # ~ if ((i == 1 and not(0 in FixedParam)) or (0 in FixedParam) and i == 0):
                # ~ limitLow[i] = iv2[i] - np.sign(iv2[i])*iv2[i]*0.1
                # ~ limitUp[i] = iv2[i] + np.sign(iv2[i])*iv2[i]*0.1
                # ~ Stepsize[i] = abs(iv2[i])*0.0001
            # ~ if (i == 0):
                # ~ limitLow[i] = iv2[i] - np.sign(iv2[i])*iv2[i]*0.05
                # ~ limitUp[i] = iv2[i] + np.sign(iv2[i])*iv2[i]*0.05
                # ~ Stepsize[i] = abs(iv2[i])*0.01 
    else:
        for i in range(0,len(iv2)):
            limitLow[i] = iv2[i] - np.sign(iv2[i])*iv2[i]*0.3
            limitUp[i] = iv2[i] + np.sign(iv2[i])*iv2[i]*0.3
            Stepsize[i] = abs(iv2[i])*0.01
            if (i in FixedParam):
                limitLow[i] = iv2[i] - np.sign(iv2[i])*iv2[i]*1e-15
                limitUp[i] = iv2[i] + np.sign(iv2[i])*iv2[i]*1e-15
                Stepsize[i] = abs(iv2[i])*1e-15
    if (not UseSplineMode):
        if (limitUp[1] > MinCharge+NSample*h.dxval()):
            limitUp[1] = MinCharge+NSample*h.dxval()-1
        if (limitLow[1] < MinCharge):
            limitLow[1] = MinCharge
        
    print(" Limits are : ")    
    print(limitUp)
    print(limitLow)
    opt_2g.set_initial_values(iv2);
    opt_2g.set_limits_lo(limitLow)
    opt_2g.set_limits_hi(limitUp)
    opt_2g.set_scale(Stepsize)
    status, xopt_2g, fval_2g = opt_2g.minimize()
    print(status)
    xerr = np.zeros(len(iv2))
    print(len(iv2))
    err_mat = np.zeros((len(iv2),len(iv2)))
    print(len(iv2))
    if (Optimizer):
        status, err_mat= opt_2g.calc_error_matrix()
        xerr = np.sqrt(err_mat.diagonal())
        print("Error Mat")
        print(status)
        print(err_mat)
        print(xerr)
        print(fval_2g)
        print(status)
    if Verbose:
        calin.plotting.plot_histogram(h, lw=2, label='SPE data')
        xlabel('Signal [DC]')
        ylabel('Events per %d DC bin'%h.dxval())
        ihist = range(0,h.nbin());
        xhist = h.all_xval_center()
        mes_2g.set_parameter_values(xopt_2g)
        ymodel_2g = \
        list(map(lambda x: h.sum_w()*h.dxval()*mes_2g.pdf_mes(x),xhist))
        if (not SimpleModel):
            plot(xhist,ymodel_2g, 'r', label='Two Gaussian')
        else:
            plot(xhist,ymodel_2g, 'r', label='Gaussian')
        a=list(axis())
    return (xopt_2g,err_mat,mes_2g,limitLow,limitUp,status,ses_2g,ped)
