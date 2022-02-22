# Author: Viktor Belay

# Purpose of script: Specify directory containing all files necessary for FAH simulation. Move all files necessary to directory in pllwskifah1

import mdtraj as md
import paramiko
from scp import SCPClient, SCPException


def send_to_fah(lilac_out_dir,fah_in_dir):

	def remove_solvent(lilac_out_dir):

		try:

			for i in os.listdir(lilac_out_dir):

				if file.endswith(('.pdb','.cif')):

					try: 
						os.chdir(lilac_out_dir)
						print('removing solvent from '+i)
						md.load(i).remove_solvent().save_pdb(lilac_out_dir+'/nosolvent_'+i)
						print('removed solvent from '+ i)
					except:

						print('something went awry :(')
				else:

					print('no topology found!')

	try:

		host = 'pllskifah1'
    	username = 'server'
    	port = 22 
    	ssh = paramiko.SSHClient()
    	ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    	ssh.connect(host,port,username=username); x=1

    except:_

    	print('something when wrong with the SSH connection to pllskifah1.')

    if x == 1:

		try:
			for i in os.listdir(lilac_out_dir):

				if file.endswith('.xml.bz2'):

					try:
						os.chdir(lilac_out_dir)
						scp=SCPClient(ssh.get_transport())
        				scp.put(i,remote_path=fah_in_dir+"/.")
        				scp.close()

        except: 

        	print('Something went wrong with transporting zipped xml files.')

        try:
        	for i in os.listdir(lilac_out_dir):

        		if file.endswith(('.pdb','.cif')):

        			try:

        				remove_solvent(lilac_out_dir)

        				os.chdir(lilac_out_dir)

        				scp=SCPClient(ssh.get_transport())
        				scp.put(i,remote_path=fah_in_dir+"/.")
        				scp.close()

        except:

        	print('Something went wrong with transporting fdude got iles.')



if name == '__main__':
	send_to_fah(sys.argv[1],sys.argv[2])
