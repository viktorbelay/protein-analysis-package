#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 15:28:39 2022

@author: belayv

# Purpose: do all relevent analyses on IP3R JD as of 05/25/22
"""

import mdtraj as md 
from matplotlib import pyplot as pl
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import pca
import seaborn as sns
import pandas as pd
import sys



def load(topology_path,topology_name,trajectory_path,trajectory_name):
    
    load = md.load(trajectory_path+trajectory_name,top=topology_path+topology_name)
    
    load.superpose(reference=load,frame=0,atom_indices=load.topology.select('protein'))
    
    universe = mda.Universe(topology_path+topology_name,trajectory_path+trajectory_name)

    
    return [load,universe]


def rmsd_all_prot(load):
    
    rmsd = md.rmsd(target=load,reference=load,frame=0,atom_indices=load.topology.select('protein and name CA'))
    
    time = [x * 1000/1000 for x in range(len(rmsd))]
    
    return [rmsd,time]
    
def rmsd_by_chain(load):
    
    rmsd1 = md.rmsd(target=load,reference=load,frame=0,atom_indices=load.topology.select('chainid 0 and protein and name CA'))
    rmsd2 = md.rmsd(target=load,reference=load,frame=0,atom_indices=load.topology.select('chainid 1 and protein and name CA'))

    
    
    time = [x * 1000/1000 for x in range(len(rmsd1))]
    
    return [rmsd1,rmsd2,time]


def rmsf_all_protein(load):
    
    rmsf = md.rmsf(target=load,reference=load,frame=0,atom_indices=load.topology.select('protein and name CA'))

    return rmsf



def pca_diagnostic(universe):
    
    ### Function purpose : will write later )
    
    
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
    
def calc_com_distances(load):
    


    # com_chainA=md.compute_center_of_mass(load1,select='(residue 2176 or residue 2177 or residue 2178 or residue 2179 or residue 2180 or residue 2181) and chainid 0')

    # com_chainB=md.compute_center_of_mass(load1,select='(residue 2580 or residue 2581 or residue 2582 or residue 2583 or residue 2584 or residue 2585 or residue 2586) and chainid 1')

    com_chainA=md.compute_center_of_mass(load,select='chainid 0')

    com_chainB=md.compute_center_of_mass(load,select='chainid 1')
    ## Get all distance arrays


    com_chainA_x = []
    com_chainA_y = []
    com_chainA_z = []

    for i in list(range(0,len(com_chainA))):

        com_chainA_x.append(com_chainA[i][0])
        com_chainA_y.append(com_chainA[i][1])
        
        com_chainA_z.append(com_chainA[i][2])
        
        
    com_chainA_x = np.array(com_chainA_x)
    com_chainA_y = np.array(com_chainA_y)
    com_chainA_z = np.array(com_chainA_z)



        
        
    com_chainB_x = []
    com_chainB_y = []
    com_chainB_z = []

    for i in list(range(0,len(com_chainB))):

        com_chainB_x.append(com_chainB[i][0])
        com_chainB_y.append(com_chainA[i][1])
        
        com_chainB_z.append(com_chainB[i][2])
        
        
    com_chainB_x=np.array(com_chainB_x)
    com_chainB_y=np.array(com_chainB_y)
    com_chainB_z=np.array(com_chainB_z)
        
    distance = np.sqrt(
        
        np.square(
            
            com_chainB_x - com_chainA_x
            
            
          ) + np.square(
              
              com_chainB_y - com_chainA_y
                        
                        ) + np.square(
                            
                            com_chainB_z - com_chainA_z
                            
                            )
        )
                            
                            
    contact = md.compute_contacts(load,contacts=[[265,343]])

    contacts=[]
    for i in list(range(0,len(contact[0]))):
    
        contacts.append(contact[0][i][0])
    
    
    np.array(contacts)
    
    time = [x * 1000/1000 for x in range(len(contacts))]
                   
    return [distance,contacts,time]


        
    
def lig_rmsd(load,lig_string):
    
    rmsd = md.rmsd(target=load,reference=load,frame=0,atom_indices=load.topology.select('resn '+lig_string))
    
    time = [x * 1000/1000 for x in range(len(rmsd))]
    
    return [rmsd,time]

    

