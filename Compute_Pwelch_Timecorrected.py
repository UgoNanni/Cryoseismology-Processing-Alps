#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 18:18:18 2018

@author: ugonanni
"""

# --------------------------------------------------------------------------
# Import packages
# --------------------------------------------------------------------------


import glob
import matplotlib.pyplot as plt
import numpy as np
import obspy
from os.path import join
from obspy.core import read
from obspy.core import UTCDateTime

from scipy.signal import welch, hanning
from scipy import signal

import os

import scipy
import scipy.fftpack

from scipy import pi
from scipy.fftpack import fft

import time

import astropy.time
import dateutil.parser

import datetime
import pandas as pd
# --------------------------------------------------------------------------
# Define path
# --------------------------------------------------------------------------

path2dat = '/media/ugonanni/2E3D-8697/DATA/ARGENTIERE/ARG_Borehole_Ugo/SACt/Corr'
path2dat = '/media/ugonanni/2E3D-8697/DATA/ARGENTIERE/ARG_DEEp_Ugo/SAC/Corr'
path2dat = '/media/ugonanni/2E3D-8697/DATA/ARGENTIERE/ARG_P7_B03/SAC/DAY/'
#path2sav = '/home/ugonanni/Share/PhD/DATA/PWelch_output'
#path2dat = '/media/ugonanni/Nanni/DATA_ARG_backup/GDA/SAC/RAW'
#path2sav = '/home/ugonanni/Share/PhD/DATA/PWelch_output'
#path2dat = '/media/ugonanni/Nanni/DATA_ARG_backup/GDA/SAC/Corrected'
#path2sav = '/media/ugonanni/Nanni/DATA_ARG_backup/GDA/PSD'
#path2dat = '/media/ugonanni/Nanni/DATA_ARG_backup/GDA/SAC/Corr'
path2sav = '/home/ugonanni/Share/PhD/DATA/PSD_0.5toNyquist_0.5dF_4sec'


filename = ''

# --------------------------------------------------------------------------
# Define parameters
# --------------------------------------------------------------------------


# Pwelch parameters
l_win = 4 # secondes
# window
window = 'hanning'
# detrend 
detrend = 'False'
# Number of segments
nb_seg = 2**19
lseg = l_win/2


# PSD vector
len_PSD = 1001


# min an max frq
fmin = [1,5,10,20]
fmax = [5,10,20,50]

# first time
time17 = astropy.time.Time(datetime.datetime(2017, 1, 1, 0, 0))
#time17 = astropy.time.Time(datetime.datetime(2000, 1, 1, 0, 0))

t17 = time17.jd-1
# --------------------------------------------------------------------------
# Compute PWelch
# --------------------------------------------------------------------------
#for ext in ('*Z*.166*', '*Z*.167*', '*Z*.168*', '*Z*.168','*Z*.169*','*Z*.17*'):

it_day = 0
#for ext in ['*01..1Z*']:
#for ext in ['*B01*D.2018.1*SAC','*B01*D.2018.2*SAC','*B01*D.2018.05*SAC','*B01*D.2018.06*SAC','*B01*D.2018.07*SAC','*B01*D.2018.08*SAC','*B01*D.2018.09*SAC']:
#for ext in ['*B02*D.2018.25*SAC','*B02*D.2018.26*SAC','*B02*D.2018.27*SAC','*B02*D.2018.28*SAC','*B02*D.2018.29*SAC']:
#for ext in ['*B01*D.2017.350*SAC']:
#for ext in ['*B02*D.2018.297*SAC','*B02*D.2018.298*SAC','*B02*D.2018.299*SAC','*B02*D.2018.3*SAC']:
#for ext in ['*B02*D.2018.343*SAC']:
#for ext in ['*B02*D.2018.34*SAC','*B02*D.2018.35*SAC','*B02*D.2018.36*SAC','*B02*D.2019*SAC']:
#for ext in ['*B03*D.2018.343*SAC']:
for ext in ['*03..1Z*']:

    nbfiles = len(glob.glob(path2dat + '/' + ext))

    
    for file in sorted(glob.glob(path2dat + '/'+ ext)):

        # read data and time
        data = read(file)
        Time_data = data[0]

 
        # define namels
        t = Time_data.stats.starttime
        name = '' + Time_data.stats.network + '.' + Time_data.stats.station + '.' + Time_data.stats.channel + '.' + str(t.year) + '.' + str(t.month) + '.' + str(t.day)
        #print(name)
        #name = '' + Time_data.stats.network + '.' + 'B02' + '.' + Time_data.stats.channel + '.' + str(t.year) + '.' + str(t.month) + '.' + str(t.day)      
        print(name)

        # Time absolute
        #pd.datetime(datime.year, datime.month, datime.day, datime.hour, datime.minute, datime.second)
        
        # File duration in secondes
        duration = Time_data.stats.endtime - Time_data.stats.starttime
        
        # Sampling frequency of the time series. 
        fs = int(np.round(Time_data.stats.sampling_rate))
        
        # number of points of the data
        nb_psd = duration/l_win
        
         # define PSD matrix
        PSD_matrix = np.empty([int(nb_psd),int(len_PSD)])
        Time_vector = []
        Hydro = np.empty([int(nb_psd),len(fmin)])  
        Hydro_std = np.empty([int(nb_psd),len(fmin)])  

        # compute PSD
        for jj in range(0,int(nb_psd)):
            # extract short time
            tmp = Time_data[jj*l_win*fs:(jj+1)*l_win*fs]
            nb_points = len(tmp)
            
            # Length of each segment. 
            nperseg = fs*lseg # here 1 secondes
            
            # Number of points to overlap between segments. 
            #noverlap = nperseg / 2 
            noverlap = nperseg / 2 

            # Time associated with Pwelch is centered on the window
            a = int(((jj*l_win))+(l_win/2)) 
            
            #Time_vector[jj] = Time_data.stats.starttime + a
            Time_vector.append(Time_data.stats.starttime + a)
          
            
            # Compute PSD with welch
            if nb_points >= nperseg:
                f, Pxx = welch(tmp, fs, window=window, nperseg=nperseg, noverlap=noverlap, detrend=False, return_onesided=True)  
                PSD_matrix[jj,:] = Pxx
            #for ll in range(0,len(fmin)-1):
                #Hydro[jj,ll] = np.mean(PSD_matrix[jj,int(fmin[ll])*lseg:int(fmax[ll])*lseg])
                #Hydro_std[jj,ll] = np.std(PSD_matrix[jj,int(fmin[ll])*lseg:int(fmax[ll])*lseg])
               
        
        Time_vector = np.asarray(Time_vector)
        
        Tdoy = np.empty(len(Time_vector))
        for jj in range(0,len(Time_vector)):
            tim = astropy.time.Time(datetime.datetime(Time_vector[jj].year, Time_vector[jj].month, Time_vector[jj].day, Time_vector[jj].hour, Time_vector[jj].minute, Time_vector[jj].second))
            Tdoy[jj] =  tim.jd
        Tdoy = Tdoy - t17
        #plt.plot(Time_data)
        #plt.show()
        
        
        print(str(it_day+1) + ' over '  + str(nbfiles))        
        it_day = it_day + 1
        # save data
        try:
            np.savetxt(path2sav + '/' + name + '_PSD_Tcorr_'  + str(l_win) +'sec',PSD_matrix)
            np.savetxt(path2sav + '/' + name + '_Time_Tcorr_' + str(l_win) +'sec',Tdoy)
            #np.savetxt(path2sav + '/' + name + '_BF_Tcorr_' + str(l_win) +'sec',Hydro)
            #np.savetxt(path2sav + '/' + name + '_BFstd_Tcorr_' + str(l_win) +'sec',Hydro_std)
        except IOError:
            print("Not enough space")
            os.system('pause')
             # We try again
            try:
                np.savetxt(path2sav + '/' + name + '_PSD_Tcorr_'  + str(l_win) +'sec',PSD_matrix)
                np.savetxt(path2sav + '/' + name + '_Time_Tcorr_' + str(l_win) +'sec',Tdoy)
                #np.savetxt(path2sav + '/' + name + '_BF_Tcorr_' + str(l_win) +'sec',Hydro)
                #np.savetxt(path2sav + '/' + name + '_BFstd_Tcorr_' + str(l_win) +'sec',Hydro_std)
            except IOError:
                print("Still not enough space")
        
        
       
        
        
    np.savetxt(path2sav + '/'+ name  + '_Frequency_Tcorr_' + str(l_win) +'sec',f)
    
    # loop on each days for concatenating
    



















