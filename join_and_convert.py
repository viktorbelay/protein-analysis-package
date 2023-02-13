#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 10:04:20 2022

@author: belayv

# Purpose: Join and convert FAH trajectories to XTC for demystifying


Parameters:
    
    Location: accepts 'local' or 'lilac', local being your local machine and lilac being the MSKCC HPC 
    project_id: accepts 'all' or any integer greater than 0. This determines how many trajectories get merged.
    count: total number of trajectories to be merged
    yaml_path: path to the yaml containing path to trajectories
    
"""

import mdtraj as md
import yaml
import os
import sys
#ogdir=os.getcwd()
def join_and_convert(location,project_id,stride_value=1,count='all',yaml_path='/Users/belayv/fah-analysis-tmem175/metadata/projects.yaml'):
    
    with open(yaml_path) as file: # Load YAML containing location of projects
        
        projects = yaml.load(file,Loader=yaml.FullLoader) 
        
    trajs_h5_run0=[]
    trajs_h5_run1=[]
    trajs_h5=[]
    index=0
    
    if location == 'local':
        os.chdir(projects[project_id]['local_fn'])
        
        if count == 'all':
            
            for h5 in os.listdir(projects[project_id]['local_fn']):
                
                if h5.startswith('run0') and h5.endswith('.h5'):
                    
            
                    trajs_h5_run0.append(h5)
                    
        
                    
                elif h5.startswith('run1') and h5.endswith('.h5'):
                    
                    trajs_h5_run1.append(h5)
                    
                    
                    
            for h5 in os.listdir(projects[project_id]['local_fn']): 
                if h5.endswith('.h5'):
                    
                    trajs_h5.append(h5)
                    
                        
                    

            cat_traj_run0=str(md.load(trajs_h5_run0,stride=stride_value))
            cat_traj_run1=str(md.load(trajs_h5_run1,stride=stride_value))
            
            traj_h5_loaded=md.load(trajs_h5,stride=stride_value)
            os.chdir(ogdir)

            
            traj_h5_loaded.save_xtc('joined_traj.xtc')
            print(traj_h5_loaded,'a concatenation of','run0  trajecories:',cat_traj_run0,'and run1 trajectories:',cat_traj_run1,'has been saved as joined_traj.xtc')
            
            
        else:
            
            for h5 in os.listdir(projects[project_id]['local_fn']):
                
                if h5.startswith('run0') and h5.endswith('.h5'):
                    
                    if index == count:
                        
                        break
                    
                    index += 1
                    
                    trajs_h5_run0.append(h5)
                    trajs_h5.append(h5)
            index=0       
            for h5 in os.listdir(projects[project_id]['local_fn']):
                
                if h5.startswith('run1') and h5.endswith('.h5'):
                    
                    if index == count:
                        
                        break
                    
                    index += 1
                    
                    trajs_h5_run1.append(h5)
                    trajs_h5.append(h5)
                    

            cat_traj_run0=str(md.load(trajs_h5_run0,stride=stride_value))
            cat_traj_run1=str(md.load(trajs_h5_run1,stride=stride_value))
            traj_h5_loaded=md.load(trajs_h5,stride=stride_value)
            os.chdir(ogdir)

            traj_h5_loaded.save_xtc('joined_traj.xtc')
            print(traj_h5_loaded,'a concatenation of','run0  trajecories:',cat_traj_run0,'and run1 trajectories:',cat_traj_run1,'has been saved as joined_traj.xtc')
                    
    elif location == 'lilac':
        os.chdir(projects[project_id]['raw_data'])
        
        if count == 'all':
            
            for h5 in os.listdir(projects[project_id]['raw_data']):
                
                if h5.startswith('run0') and h5.endswith('.h5'):
                    
            
                    trajs_h5_run0.append(h5)
                    
                    
                    
                elif h5.startswith('run1') and h5.endswith('.h5'):
                    
                    trajs_h5_run1.append(h5)
                    
                    
                    
            for h5 in os.listdir(projects[project_id]['raw_data']): 
                if h5.endswith('.h5'):
                    
                    trajs_h5.append(h5)
                    
                        
                    

            cat_traj_run0=str(md.load(trajs_h5_run0,stride=stride_value))
            cat_traj_run1=str(md.load(trajs_h5_run1,stride=stride_value))
            
            traj_h5_loaded=md.load(trajs_h5,stride=stride_value)
            os.chdir(ogdir)

            traj_h5_loaded.save_xtc('joined_traj.xtc')
            print(traj_h5_loaded,'a concatenation of','run0  trajecories:',cat_traj_run0,'and run1 trajectories:',cat_traj_run1,'has been saved as joined_traj.xtc')
            
            
        else:
            
            for h5 in os.listdir(projects[project_id]['raw_data']):
                
                if h5.startswith('run0') and h5.endswith('.h5'):
                    
                    if index == count:
                        
                        break
                    
                    index += 1
                    
                    trajs_h5_run0.append(h5)
                    trajs_h5.append(h5)
            index=0       
            for h5 in os.listdir(projects[project_id]['raw_data']):
                
                if h5.startswith('run1') and h5.endswith('.h5'):
                    
                    if index == count:
                        
                        break
                    
                    index += 1
                    
                    trajs_h5_run1.append(h5)
                    trajs_h5.append(h5)
                    

            cat_traj_run0=str(md.load(trajs_h5_run0,stride=stride_value))
            cat_traj_run1=str(md.load(trajs_h5_run1,stride=stride_value))
            traj_h5_loaded=md.load(trajs_h5,stride=stride_value)
            os.chdir(ogdir)
            traj_h5_loaded.save_xtc('joined_traj.xtc')
            print(traj_h5_loaded,'a concatenation of','run0  trajecories:',cat_traj_run0,'and run1 trajectories:',cat_traj_run1,'has been saved as joined_traj.xtc')
            


if __name__=='__main__':
    
    ogdir=os.getcwd()

    if sys.argv[4] == 'all':
        
        join_and_convert(sys.argv[1], int(sys.argv[2]),int(sys.argv[3]),sys.argv[4],sys.argv[5])
        
    else:
        
        join_and_convert(sys.argv[1], int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),sys.argv[5])
