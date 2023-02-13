"""
Created on Thu Aug 18 10:04:20 2022

@author: belayv

# Purpose: Join and convert short Lilac trajectories into one large trajectory for featurization and/or demystifying


Parameters:

NOTE: add a feature to convert to XTC on the fly
 
    
"""

import mdtraj as md
import os
import sys

def join_lite(working_dir,stride_value=1):

    trajs = []

    os.chdir(working_dir)



    for pdb in os.listdir(working_dir):

        if pdb.endswith('.pdb'):

            pdb_path=pdb



    for dcd in os.listdir(working_dir):

        if dcd.endswith('.dcd'):

            trajs.append(md.load(dcd,top=pdb_path,stride=stride_value))


    cat_traj=md.join(trajs)

    cat_traj.save_dcd(working_dir+'/joined.dcd')

    print(cat_traj,'has been saved as joined.dcd')


if __name__=='__main__':

    join_lite(sys.argv[1],int(sys.argv[2]))



