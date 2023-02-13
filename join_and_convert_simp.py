"""
Created on Thu Aug 18 10:04:20 2022

@author: belayv

# Purpose: Join and convert Lilac trajectories to XTC for demystifying


Parameters:
    

    
"""


import mdtraj as md
import yaml
import os
import sys

def join_and_convert_simp(location,stride_value=1,count='all'):

    trajs=[]
    index=0

    os.chdir(location)

    if count == 'all':

        for dcd in os.listdir(location):

            if dcd.endswith('.dcd'):

                trajs.append(md.load_dcd(dcd,top='../equil/openmm/step3_input.pdb',stride=stride_value))

        

        trajs_joined = md.join(trajs)

        os.chdir(ogdir)

        trajs_joined.save_xtc('joined_traj.xtc')

    else:

        pass


if __name__=='__main__':
    
    ogdir=os.getcwd()

    if sys.argv[3] == 'all':
        
        join_and_convert_simp(sys.argv[1], int(sys.argv[2]),sys.argv[3])
        
    else:
        
        join_and_convert_simp(sys.argv[1], int(sys.argv[2]),int(sys.argv[3]))
