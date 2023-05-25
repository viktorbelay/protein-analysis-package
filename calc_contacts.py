#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 11:46:14 2022

@author: belayv

# Compute distance between the center of mass of two polypeptides

Note: currently only works for IP3R JD

Note 2: in the future, calculate any arbritrary distance
"""

import mdtraj as md
import numpy as np
from matplotlib import pyplot as pl

def calc_com_distances(top_path,top_name,traj_path,traj_name):

    # Write a description of what this function does
    # Write a description of what each input is and what it should be
    
    load1=md.load(traj_path+traj_name,top=top_path+top_name)
    
    load1.superpose(reference=load1,frame=0,atom_indices=load1.topology.select('protein'))

    # com_chainA=md.compute_center_of_mass(load1,select='(residue 2176 or residue 2177 or residue 2178 or residue 2179 or residue 2180 or residue 2181) and chainid 0')

    # com_chainB=md.compute_center_of_mass(load1,select='(residue 2580 or residue 2581 or residue 2582 or residue 2583 or residue 2584 or residue 2585 or residue 2586) and chainid 1')

    com_chainA=md.compute_center_of_mass(load1,select='chainid 0')

    com_chainB=md.compute_center_of_mass(load1,select='chainid 1')
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
                            
                            
    contact = md.compute_contacts(load1,contacts=[[265,343]])

    contacts=[]
    for i in list(range(0,len(contact[0]))):
    
        contacts.append(contact[0][i][0])
    
    
    np.array(contacts)
                   
    return [distance,contacts]


topo_path='/Users/belayv/Downloads/'
topo_name='atptest5.pdb'
traj_path='/Users/belayv/Downloads/'
traj_name='atptest5.dcd'

topo_path1='/Users/belayv/Downloads/'
topo_name1='test8apo.pdb'
traj_path1='/Users/belayv/Downloads/'
traj_name1='test8apo.dcd'

dist1=calc_com_distance(topo_path, topo_name, traj_path, traj_name)
dist2=calc_com_distance(topo_path1, topo_name1, traj_path1, traj_name1)

pl.hist(dist1[0],alpha=1,bins=25)

pl.hist(dist2[0],alpha=0.5,bins=25)

pl.figure()
pl.hist(dist1[1],alpha=1,bins=25)

pl.hist(dist2[1],alpha=0.5,bins=25)

# create a list of 100 numbers

