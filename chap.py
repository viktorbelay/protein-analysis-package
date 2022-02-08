#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################################################
############################################################

# Author: Viktor Belay 
# Creation date: 02/06/2022

# Purpose of script: Save CHAP data as numpy arrays and plot things in a pretty way :3
# Also, can optionally get output.json (next version)

# Version: 0.1.0

# Funtionalty:
    
    # Inputs:
        # Required: path to topology
        # Required: path to trajectory
        # Optional: plot (True or False). Will plot results with no input
    # Outputs:
        
# Notes: Currently configured only for TMEM175
# Future plans: Configure script to work with any channel

############################################################
############################################################

import numpy as np
from matplotlib import pyplot as pl
import matplotlib.mlab as mlab
import json
from scipy.stats import norm
import paramiko
import mdtraj as md
import os
import sys
import subprocess
import getpass
from scp import SCPclient, SCPException


def do_chap(topology,trajectory,output_path):
    
    def dcd_to_xtc(topology, trajectory):
        
        dcd_traj = md.load(trajectory,top=topology)
        dcd_traj.save_xtc(output_path+'/trajectory.xtc')
        
        # Under development
        
def load_data(data):
    
    with open(data) as data_file:
        chap_data=json.load(data_file)
        
    return chap_data

    
        
def get_pore_radius_profile(chap_dat,plot='True'):
    
    pass

def get_min_pore_radius(chap_data,plot_hist=True,plot_time_series=True):
    dat=[]
    min_radius = np.array(chap_data['pathwayScalarTimeSeries']['minRadius'])*10 # Min pore (nm) radians * 10 = min pore radius in A 
    t = np.array(chap_data['pathwayScalarTime Series']['t'])
    
    dat.append(t)
    dat.append(min_radius)
    
    dat=np.array(dat)
    
    return dat

    if plot_hist == True:
        
        if plot_time_series==True:
            
            pl.figure()
            
            mean,std=norm.fit(min_radius[min_radius>=0])
            pl.hist(min_radius[min_radius>=0],bins=(len(min_radius)/4),alpha=0.9)
            pl.xlabel('Minimum pore radius (A)')
            pl.ylabel('Frequency')
            pl.title('mean= '+str(mean)+' SD= '+str(std))
            
            pl.figure()
            
            pl.plot(t,min_radius)
            pl.xlabel('time (ns)')
            pl.ylabel('minimum pore radius (A)')
            
            pl.ylim(0,max(min_radius)+2)
            
        elif plot_time_series==False:
            
            pl.figure()
            
            mean,std=norm.fit(min_radius[min_radius>=0])
            pl.hist(min_radius[min_radius>=0],bins=(len(min_radius)/4),alpha=0.9)
            pl.xlabel('Minimum pore radius (A)')
            pl.ylabel('Frequency')
            pl.title('mean= '+str(mean)+' SD= '+str(std))
            
    elif plot_hist==False:
        
        if plot_time_series==True:
            
            pl.figure()
            
            pl.plot(t,min_radius)
            pl.xlabel('time (ns)')
            pl.ylabel('minimum pore radius (A)')
            
            pl.ylim(0,max(min_radius)+2)
            
        elif plot_time_series==False:
            
            pass
        
        
        
def get_pore_solvent_density(chap_data,plot=True):
    
    dat = []
    z = np.array(chap_data['pathwayProfile']['s'])*10
    mean_solvent_density = np.array(chap_data['pathwayProfile']['densityMean'])
    
    dat.append(z)
    dat.append(mean_solvent_density)
    
    dat = np.array(dat)
    
    if plot==True:
        
        pl.figure()
        pl.plot(z,mean_solvent_density)
    
    
    
    
    
        
        
    
    
    
    
    
    
    
    

