import mdtraj as md
#import os
#import sys


def dcd_to_xtc(dcd_to_convert,topology_path):

	dcd_file=md.load_dcd(dcd_to_convert,top=topology_pathf)

	print(dcd_file,'has been loaded... converting to XTC')

	md.save_xtc('traj.xtc')

	print(dcd_file,'has ben saved as traj.xtc')


if __name__=='__main__':

    dcd_to_xtc(sys.argv[1],sys.argv[2])

