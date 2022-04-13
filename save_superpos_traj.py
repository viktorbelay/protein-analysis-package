#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 12:08:31 2022

@author: belayv
"""

import mdtraj as md
import sys

def save_superpos_traj(topology,trajectory):
    
    traj = md.load(trajectory,top=topology)
    pdb = md.load_pdb(topology)
    traj.superpose(reference=pdb,frame=0)
    traj.save('superposed_trajectory.dcd')
    
if __name__ == '__main__':
    
    save_superpos_traj(sys.argv[1],sys.argv[2])