#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 21 12:29:28 2022

@author: belayv

# Purpose of script: IP3R munging
"""

import mdtraj as md 

def super(topology_path,trajectory_path,trajectory_name='equilibrated.dcd',topology_name='step3_input.pdb'):
    
    load = md.load(trajectory_path+trajectory_name,top=topology_path+topology_name)
    
    load.superpose(reference=load,frame=0,atom_indices='protein')
    
    load.save_dcd(trajectory_path+'super_equilibrated.dcd')
    
    
    
if __name__=='__main__':
    
    super(topology_path,trajectory_path,trajectory_name='equilibrated.dcd',topology_name='step3_input.pdb')