"""
Created on Thu Aug 18 10:04:20 2022

@author: belayv

# Purpose: Strip IP3R JD trajectories of everything except protein (waters,ions,ligand) and save result


Parameters:
    
traj: path to trajectory
topo: path to related topology
    
"""

import mdtraj as md
import os
import sys

def strip(traj,topology):

    a_traj=md.load(traj,top=topology)
    a_topo=md.load(topology)
    a_topo_protein=a_topo.atom_slice(a_topo.topology.select('chainid 0 or chainid 1'))

    a_traj_protein=a_traj.atom_slice(a_traj.topology.select('chainid 0 or chainid 1'))

    os.chdir(ogdir)

    a_traj_protein.save_dcd('joined_stripped.dcd')
    a_topo_protein.save_pdb('joined_stripped.pdb')

    print(a_traj,'stripped of waters, ions, and ligand saved as',a_traj_protein,'filename joined_stripped.dcd')
    print(a_topo,'stripped of waters, ions, and ligand saved as',a_topo_protein,'filename joined_stripped.pdb')

if __name__=='__main__':

    ogdir=os.getcwd()

    strip(sys.argv[1],sys.argv[2])







