#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 13:12:27 2022

@author: belayv

# Purpose of script: calculate fluctuation of ligand species bound with JD
"""

import mdtraj as md 
from matplotlib import pyplot as pl
import numpy as np

## Testing area ####
topo_path='/Users/belayv/Downloads/'
topo_name='camp_pdb_1.pdb'
traj_path='/Users/belayv/Downloads/'
traj_name='camp_500ns.dcd'





load=md.load(traj_path+traj_name,top=topo_path+topo_name)

load.superpose(reference=load,frame=0,atom_indices=load.topology.select('protein'))


## Calculate ATP RMSD over time for all atoms

rmsd = md.rmsd(target=load,reference=load,frame=0,atom_indices=load.topology.select('resn ATP'))

## Calculate RMSF 

rmsf = md.rmsf(target=load, reference=load,frame=0,atom_indices=load.topology.select('resn ATP'))