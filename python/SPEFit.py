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

import csv
import pickle
import sys

PedestalGuess = 0

def PMax(r):
    if (np.pi*r**2/(np.pi*r**2 + np.pi - 2*r**2 - 2) <= 1):
        return np.pi*r**2/(np.pi*r**2 + np.pi - 2*r**2 - 2)
    else:
        return 1
def UncertaintySig1mu2(n,pp,res,mu1,mu2):
    return (-2*np.sqrt(2)*np.pi*pp*res**2*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*np.sqrt(1/np.pi)/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + (8*np.pi**2*mu2*pp**2*res**4*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2) - (2*mu2*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2 - (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(2*mu2*n**2*(res**2 + 1)*((-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))*((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2 - (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)/(res**2 + 1)) - 8*np.pi**2*pp**2*res**4*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))/((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))) + 2*mu2)/(res**2 + 1))*(-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))/2)/np.sqrt(8*np.pi**2*mu2**2*pp**2*res**4*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2) - (mu2**2*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2 - (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(mu2**2*n**2*(res**2 + 1)*((-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))*((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2 - (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)/(res**2 + 1)) - 8*np.pi**2*pp**2*res**4*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))/((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))) + mu2**2)/(res**2 + 1))*(-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))))/(-2*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 4*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))
def UncertaintySig1res(n,pp,res,mu1,mu2):
    return (-2*np.sqrt(2)*np.pi*mu2*pp*res**2*(-2*np.pi*res + 4*res)*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*np.sqrt(1/np.pi)/(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2 - 2*np.sqrt(2)*np.pi*mu2*pp*res**2*(-np.pi*pp*res**2*(-2*np.pi*res + 4*res)/(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2 - 2*np.pi*pp*res/(-2*res**2 + np.pi*res**2 - 2 + np.pi))*np.sqrt(1/np.pi)/(-2*res**2 + np.pi*res**2 - 2 + np.pi) - 4*np.sqrt(2)*np.pi*mu2*pp*res*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*np.sqrt(1/np.pi)/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + (4*np.pi**2*mu2**2*pp**2*res**4*(-4*np.pi*res + 8*res)*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**3) + 4*np.pi**2*mu2**2*pp**2*res**4*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(-2*np.pi*pp*res**2*(-2*np.pi*res + 4*res)/(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2 - 4*np.pi*pp*res/(-2*res**2 + np.pi*res**2 - 2 + np.pi))/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2) + 16*np.pi**2*mu2**2*pp**2*res**3*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2) - (mu2**2*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2 - (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(mu2**2*n**2*(res**2 + 1)*((-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))*((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2 - (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)/(res**2 + 1)) - 8*np.pi**2*pp**2*res**4*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))/((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))) + mu2**2)/(res**2 + 1))*(8*np.pi*pp*res**3/((res**2 + 1)**2*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) - 4*np.pi*pp*res**2*(-2*np.pi*res + 4*res)/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2) - 8*np.pi*pp*res/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4*(-4*np.pi*res + 8*res)/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**3) + 32*np.pi**2*pp**2*res**3/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))/2 - (-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))*(mu2**2*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(-2*np.pi*pp*res**2*(-2*np.pi*res + 4*res)/(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2 - 4*np.pi*pp*res/(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 2*res*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(mu2**2*n**2*(res**2 + 1)*((-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))*((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2 - (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)/(res**2 + 1)) - 8*np.pi**2*pp**2*res**4*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))/((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))) + mu2**2)/(res**2 + 1)**2 - (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(2*mu2**2*n**2*res*((-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))*((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2 - (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)/(res**2 + 1)) - 8*np.pi**2*pp**2*res**4*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))/((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))) + mu2**2*n**2*(res**2 + 1)*((-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))*((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2 - (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)/(res**2 + 1)) - 8*np.pi**2*pp**2*res**4*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))*(-8*np.pi*pp*res**3/((res**2 + 1)**2*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 4*np.pi*pp*res**2*(-2*np.pi*res + 4*res)/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2) + 8*np.pi*pp*res/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) - 8*np.pi**2*pp**2*res**4*(-4*np.pi*res + 8*res)/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**3) - 32*np.pi**2*pp**2*res**3/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))/((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))**2) + mu2**2*n**2*(res**2 + 1)*((-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))*((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2 - (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)/(res**2 + 1)) - 8*np.pi**2*pp**2*res**4*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))*(np.pi*pp*res**2*(-2*np.pi*res + 4*res)/(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2 + 2*np.pi*pp*res/(-2*res**2 + np.pi*res**2 - 2 + np.pi))/((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2*(-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))) + mu2**2*n**2*(res**2 + 1)*((-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))*(2*res*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)/(res**2 + 1)**2 + (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(-2*np.pi*pp*res**2*(-2*np.pi*res + 4*res)/(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2 - 4*np.pi*pp*res/(-2*res**2 + np.pi*res**2 - 2 + np.pi)) - (-np.pi*pp*res**2*(-2*np.pi*res + 4*res)/(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2 - 2*np.pi*pp*res/(-2*res**2 + np.pi*res**2 - 2 + np.pi))/(res**2 + 1)) + ((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2 - (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)/(res**2 + 1))*(8*np.pi*pp*res**3/((res**2 + 1)**2*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) - 4*np.pi*pp*res**2*(-2*np.pi*res + 4*res)/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2) - 8*np.pi*pp*res/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4*(-4*np.pi*res + 8*res)/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**3) + 32*np.pi**2*pp**2*res**3/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2)) - 8*np.pi**2*pp**2*res**4*(-4*np.pi*res + 8*res)*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**3) - 8*np.pi**2*pp**2*res**4*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(-2*np.pi*pp*res**2*(-2*np.pi*res + 4*res)/(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2 - 4*np.pi*pp*res/(-2*res**2 + np.pi*res**2 - 2 + np.pi))/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2) - 32*np.pi**2*pp**2*res**3*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))/((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))))/(res**2 + 1) - (-np.pi*pp*res**2*(-2*np.pi*res + 4*res)/(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2 - 2*np.pi*pp*res/(-2*res**2 + np.pi*res**2 - 2 + np.pi))*(mu2**2*n**2*(res**2 + 1)*((-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))*((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2 - (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)/(res**2 + 1)) - 8*np.pi**2*pp**2*res**4*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))/((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))) + mu2**2)/(res**2 + 1))/2)/np.sqrt(8*np.pi**2*mu2**2*pp**2*res**4*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2) - (mu2**2*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2 - (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(mu2**2*n**2*(res**2 + 1)*((-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))*((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2 - (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)/(res**2 + 1)) - 8*np.pi**2*pp**2*res**4*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))/((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))) + mu2**2)/(res**2 + 1))*(-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))))/(-2*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 4*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2)) + (-2*np.sqrt(2)*np.pi*mu2*pp*res**2*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*np.sqrt(1/np.pi)/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + np.sqrt(8*np.pi**2*mu2**2*pp**2*res**4*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2) - (mu2**2*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2 - (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(mu2**2*n**2*(res**2 + 1)*((-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))*((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2 - (-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)/(res**2 + 1)) - 8*np.pi**2*pp**2*res**4*(-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)**2/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))/((-np.pi*pp*res**2/(-2*res**2 + np.pi*res**2 - 2 + np.pi) + 1)*(-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))) + mu2**2)/(res**2 + 1))*(-4*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 8*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))))*(-4*np.pi*pp*res**3/((res**2 + 1)**2*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 2*np.pi*pp*res**2*(-2*np.pi*res + 4*res)/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2) + 4*np.pi*pp*res/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) - 4*np.pi**2*pp**2*res**4*(-4*np.pi*res + 8*res)/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**3) - 16*np.pi**2*pp**2*res**3/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))/(-2*np.pi*pp*res**2/((res**2 + 1)*(-2*res**2 + np.pi*res**2 - 2 + np.pi)) + 4*np.pi**2*pp**2*res**4/(np.pi*(-2*res**2 + np.pi*res**2 - 2 + np.pi)**2))**2

def GainUncertainty(n,pp,res,mu2,Dmu2,Dres,mu1=0):
    return np.sqrt(((1-pp*PMax(res))*Dmu2)**2 + ((2*pp*PMax(res)/np.sqrt(2*np.pi))*np.sqrt((UncertaintySig1mu2(res,pp,res,mu1,mu2)*Dmu2)**2 + (UncertaintySig1res(res,pp,res,mu1,mu2)*Dres)**2))**2)

def Fit_1_gauss(h,hped,HV):
    import time
    ftime0 = time.time()
    K = 6.68060614323e-13
    res = 0.3
    iv = [1.1, 2800, 15, 1.1*K*HV**(4.68516277625), res*K*HV**(4.68516277625)]
    FixedParam = [-1]
    a,aerror,b,c,d,e,f,g = Fit_1_gauss_bis(h,hped,iv,True,False,60,FixedParam)
    return (a,aerror,b,f,g)
    
    
def Fit_2_gauss(h,hped,HV,UseHped=False,RobustMode=False,FreeMode = False, LightFixed = -1):
    import time
    global PedestalGuess
    histTable = np.zeros(h.size())
    histCharge = np.zeros(h.size())
    global MinCharge
    global NSample
    
    for i in range(h.size()):
        histTable[i] = h.weight(i)
        histCharge[i] = h.xval_center(i)
    MinCharge = histCharge[0]
    NSample = len(histCharge)
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
    luminosity = LuminosityGuess
    PedestalMean = PedestalGuess
    #print(K*HV**(4.71080288507))
    if RobustMode:
        iv = [luminosity, MeanPedestalGuess, PedestalRMSGuess,pp, res, GainGuess,0.715] 
        FixedParam = [3,6]
        if FreeMode:
           FixedParam = [-1]
        if (LightFixed > 0):
            iv[0] = LightFixed
            FixedParam = [0]
    else:    
        iv = [luminosity, MeanPedestalGuess, PedestalRMSGuess,pp,0, res, GainGuess,0.715] #1150V
        FixedParam = [3,4,7]
        if (LightFixed > 0):
            iv[0] = LightFixed
            FixedParam = [0,3,4,7]
    
    
    #FixedParam = [4]
    #FixedParam = [-1]
    
    end = True
    while(end and RobustMode):
        a,aerror,b,c,d,e,f,g = Fit_2_gauss_bis(h,hped,iv,False,False,2,FixedParam,UseHped)
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
            a,aerror,b,c,d,e,f,g = Fit_2_gauss_bis(h,hped,iv,False,False,5,FixedParam,UseHped)
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
            a,aerror,b,c,d,e,f,g = Fit_2_gauss_bis(h,hped,iv,False,False,20,FixedParam,UseHped)
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
            for i in range(len(FixedParam)):
                if (FixedParam[i] > 4):
                    FixedParam[i] += 1
            FixedParam.append(4)
            print("Fixed Param")
            print(FixedParam)
            iv = [iv[0], iv[1], iv[2],iv[3],0, iv[4], iv[5],iv[6]]

    if (not RobustMode):
        a,aerror,b,c,d,e,f,g = Fit_2_gauss_bis(h,hped,iv,True,False,60,FixedParam,UseHped)
    return (a,aerror,b,f,g)
        
        
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
        limitUp[i] = iv[i] + np.sign(iv[i])*iv[i]*0.95
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

    
def Fit_2_gauss_bis(h,hped,iv,Optimizer=False,Verbose=False,WallTime=2,FixedParam=[-1],UseHped=False):
    nSpline = len(iv)-7
    print(iv)
    
    ped = calin.math.pdf_1d.BinnedGaussianPDF(h.dxval())
    ped.set_parameter_values(np.asarray([0, 50]))
    
    ses_2g = calin.calib.pmt_ses_models.TwoGaussianSESConstrained(h.dxval())
    ses_2g.set_parameter_values(np.asarray([0.2,0.5,1,0.5]))
    
    
    iv2 = []
    if (Optimizer):
        for i in range(0,len(iv)):
            if not(i in FixedParam and i > 2):
                iv2 = iv2 + [iv[i]]
            else:
                if (i > 2):
                    ses_2g.remove_parameter_from_subspace(i-3,iv[i])
    else:
        iv2=iv
    print("iv2")
    print(iv2)
    global MinCharge
    global NSample
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
    if(Optimizer):
        for i in range(0,len(iv2)):
            limitLow[i] = iv2[i] - np.sign(iv2[i])*iv2[i]*0.95
            limitUp[i] = iv2[i] + np.sign(iv2[i])*iv2[i]*100
            Stepsize[i] = abs(iv2[i])*0.01
            if (i == 3 and not(5 in FixedParam)):
                limitLow[i] = iv2[i] - np.sign(iv2[i])*iv2[i]*0.08
                limitUp[i] = iv2[i] + np.sign(iv2[i])*iv2[i]*0.08
                Stepsize[i] = abs(iv2[i])*0.01
            if (i == 1):
                limitLow[i] = iv2[i] - np.sign(iv2[i])*iv2[i]*0.15
                limitUp[i] = iv2[i] + np.sign(iv2[i])*iv2[i]*0.15
                Stepsize[i] = abs(iv2[i])*0.01  
            if ((i in FixedParam) and i < 3):
                limitLow[i] = iv2[i] - np.sign(iv2[i])*iv2[i]*1e-15
                limitUp[i] = iv2[i] + np.sign(iv2[i])*iv2[i]*1e-15
                Stepsize[i] = abs(iv2[i])*1e-15
            #if (i == 0):
            #    limitLow[i] = iv2[i] - np.sign(iv2[i])*iv2[i]*0.05
            #    limitUp[i] = iv2[i] + np.sign(iv2[i])*iv2[i]*0.05
            #     Stepsize[i] = abs(iv2[i])*0.01 
    else:
        for i in range(0,len(iv2)):
            limitLow[i] = iv2[i] - np.sign(iv2[i])*iv2[i]*0.1
            limitUp[i] = iv2[i] + np.sign(iv2[i])*iv2[i]*0.1
            Stepsize[i] = abs(iv2[i])*0.01
            if (i in FixedParam):
                limitLow[i] = iv2[i] - np.sign(iv2[i])*iv2[i]*1e-15
                limitUp[i] = iv2[i] + np.sign(iv2[i])*iv2[i]*1e-15
                Stepsize[i] = abs(iv2[i])*1e-15
            
    
    opt_2g.set_initial_values(iv2);
    opt_2g.set_limits_lo(limitLow)
    opt_2g.set_limits_hi(limitUp)
    opt_2g.set_scale(Stepsize)
    status, xopt_2g, fval_2g = opt_2g.minimize()
    print(status)
    xerr = np.zeros(len(iv2))
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
        if (not SimpleModel):
            plot(xhist,ymodel_2g, 'r', label='Two Gaussian')
        else:
            plot(xhist,ymodel_2g, 'r', label='Gaussian')
        a=list(axis())
    return (xopt_2g,xerr,mes_2g,limitLow,limitUp,status,ses_2g,ped)
