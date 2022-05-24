#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 21 13:06:12 2022

@author: belayv

# Purpose of script: do basic analyses for IP3R JD

# Note: default submission assumes that trajectory is superposed already
"""

import mdtraj as md 
from matplotlib import pyplot as pl
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import pca
import seaborn as sns
import pandas as pd
import sys


    
def rmsd_all_prot(topology_path,trajectory_path,topology_name,trajectory_name):
    
    #load = md.load(trajectory_path+trajectory_name,top=topology_path+topology_name)
    
    # All protein RMSD (time series and hist)
    
    rmsd = md.rmsd(target=load,reference=load,frame=0,atom_indices=load.topology.select('protein and name CA'))
    
    time = [x * 1000/1000 for x in range(len(rmsd))]
    
    pl.figure()
    pl.plot(time,rmsd*10 # Convert nm to A
            );
    
    pl.xlabel('time (ns)')
    pl.ylabel('RMSD (A)')
    
    pl.legend()
    pl.title('RMSD all protein')
    
    pl.savefig(trajectory_path+'all_protein_rmsd_timeseries.png', dpi=300)
    pl.close()
    
    pl.figure()
    pl.hist(rmsd*10)
    pl.xlabel('RMSD (A)')
    pl.ylabel('Frequency')
    pl.savefig(trajectory_path+'all_protein_rmsd_hist.png')
    pl.close()
    

    
def rmsd_by_chain(topology_path,trajectory_path,topology_name,trajectory_name):
    
    rmsd1 = md.rmsd(target=load,reference=load,frame=0,atom_indices=load.topology.select('chainid 0 and protein and name CA'))
    rmsd2 = md.rmsd(target=load,reference=load,frame=0,atom_indices=load.topology.select('chainid 1 and protein and name CA'))

    
    
    time = [x * 1000/1000 for x in range(len(rmsd1))]
    
    pl.figure()
    
    pl.plot(time,rmsd1*10,label='chain A (polypeptide into TM domain)')
    pl.plot(time,rmsd2*10,label='chain B (polypeptide out of TM domain)')
    pl.legend()
    
    pl.xlabel('time (ns)')
    pl.ylabel('RMSD (A)')
    pl.title('RMSD by chain')
    
    pl.savefig(trajectory_path+'protein_rmsd_by_chain_timeseries.png',dpi=300)
    
    pl.close()
    
    
    pl.figure()
    
    pl.hist(rmsd1,label='chain A (polypeptide into TM domain)')
    pl.hist(rmsd2,label='chain B (polypeptide out of TM domain)')
    pl.legend()
    pl.savefig(trajectory_path+'protein_rmsd_by_chain_histogram')
    pl.close()
    
    
def rmsf_all_protein(topology_path,trajectory_path,topology_name,trajectory_name):
    
    rmsf = md.rmsf(target=load,reference=load,frame=0,atom_indices=load.topology.select('protein and name CA'))

    pl.figure()
    
    pl.plot(rmsf*10)
    pl.xlabel('residue')
    pl.ylabel('RMSF (A)')
    pl.ylim(0,10)
    pl.title('RMSF all residues')
    
    pl.savefig(trajectory_path+'protein_rmsf.png',dpi=300)
    
    pl.close()
def pca_diagnostic(topology_path,topology_name,trajectory_path,trajectory_name):
    pc = pca.PCA(universe, select='backbone',
             align=True, mean=None,
             n_components=None).run()
    
    backbone=universe.select_atoms('backbone')
    
    transformed = pc.transform(backbone, n_components=5)
    transformed.shape
    
    df = pd.DataFrame(transformed,
                  columns=['PC{}'.format(i+1) for i in range(5)])
    df['Time (ps)'] = df.index * universe.trajectory.dt * 21.5
    df.head()
    
    
    pl.scatter(df['PC1'],df['PC2'],c=df['Time (ps)'])
    pl.colorbar().set_label('time (ns)')
    pl.xlabel('PC1')
    pl.ylabel('PC2')
    pl.title('Trajectory PCA, first two PCs')
    
    pl.savefig(trajectory_path+'pca_analysis_first_two.png',dpi=300)
    pl.close()
    
    
def diagnostic_analysis(topology_path,topology_name,trajectory_path,trajectory_name,pcaa=True):
    
    global load
    global universe
    
    load = md.load(trajectory_path+trajectory_name,top=topology_path+topology_name)
    load.superpose(reference=load,frame=0,atom_indices='protein')
    universe = mda.Universe(topology_path+topology_name,trajectory_path+trajectory_name)
    
    rmsd_all_prot(topology_path,trajectory_path,topology_name,trajectory_name)
    
    rmsd_by_chain(topology_path,trajectory_path,topology_name,trajectory_name)
    
    rmsf_all_protein(topology_path,trajectory_path,topology_name,trajectory_name)
    
    if pcaa == True:
        pca_diagnostic(topology_path,trajectory_path,topology_name,trajectory_name)
        
    
    
if __name__=='__main__':
    
    diagnostic_analysis(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
    
    
### Testing area ####
# topo_path='/Users/belayv/Downloads/'
# topo_name='test7_zn.pdb'
# traj_path='/Users/belayv/Downloads/'
# traj_nam='test7_zn.dcd'

# diagnostic_analysis(topo_path,topo_name,traj_path,traj_nam,True)
    
    
    
    
    
    