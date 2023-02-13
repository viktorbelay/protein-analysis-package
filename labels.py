#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 11:45:12 2022

@author: belayv

Purpose: Quickly generate labels.txt file for use in demystifying

Use: 
    num_traj: The number of trajectories in a SINGLE class (eg, the number of trajectories from one RUN)
    label: The label you want for this particular class (eg, 0 or 1 for a two-state system)
    file: [].csv
    
    User should manually concatenate a label.txt containing multiple class labels.
"""
import csv
import os
import sys

def labels(num_traj,label,filename):
    
    #trajs=list(range(start,num_traj+1))
    labels=[label]*num_traj
    
    to_file= zip(labels)
    
    os.chdir(os.getcwd())
    with open(filename, 'w') as f:
        
        writer = csv.writer(f,delimiter='\t')
        writer.writerows(to_file)
        
        
        
if __name__=='__main__':
    labels(int(sys.argv[1]),int(sys.argv[2]),sys.argv[3])
    
    
