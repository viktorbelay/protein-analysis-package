# Purpose of script: remove specified residues from trajectory and PDB file


import mdtraj as md
import os
import sys

def remove_residues(top,traj,residue_start,residue_end=None):

    loaded_top=md.load(top)
    loaded_traj=md.load(traj,top=top)

    if residue_end==None:

        topo_sliced=loaded_top.atom_slice(loaded_top.topology.select('residue '+str(residue_start)))
        traj_sliced=loaded_traj.atom_slice(loaded_traj.topology.select('residue '+str(esidue_start)))

        os.chdir(ogdir)

        traj_slieced.save_xtc('joined_stripped_sliced.xtc')
        topo_sliced.save_pdb('joined_stripped_sliced.pdb')

    else:

        topo_sliced=loaded_top.atom_slice(loaded_top.topology.select('residue '+str(residue_start) + 'to '+str(residue_end)))
        traj_sliced=loaded_traj.atom_slice(loaded_traj.topology.select('residue '+str(residue_start) + 'to '+str(residue_end)))

        os.chdir(ogdir)

        traj_slieced.save_xtc('joined_stripped_sliced.xtc')
        topo_sliced.save_pdb('joined_stripped_sliced.pdb')


if __name__=='__main__':

        ogdir=os.getcwd()

        remove_residues(sys.argv[1],sys.argv[2],sys.argv[3])






