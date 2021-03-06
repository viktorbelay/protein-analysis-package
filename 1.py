#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 11:05:03 2022

@author: belayv
"""

import mdtraj as md
import paramiko
from scp import SCPClient, SCPException
import os
import sys


def send_to_fah(lilac_out_dir,fah_in_dir):
    
    def remove_solvent(lilac_out_dir):
        
        for i in os.listdir(lilac_out_dir):

            if i.endswith('.pdb'):

                try:
                    os.chdir(lilac_out_dir)
                    print('removing solvent from '+i)
                    md.load(i).remove_solvent().save_pdb(lilac_out_dir+'/nosolvent_equilibrated.pdb')
                    print('removed solvent from '+ i)
                except:
                    print('something went awry :(')
                    
            elif i.endswith('.cif'):
                
                try:
                    os.chdir(lilac_out_dir)
                    print('removing solvent from '+i)
                    md.load(i).save_pdb(lilac_out_dir+i)
                    md.load(i).remove_solvent().save_pdb(lilac_out_dir+'/nosolvent_equilibrated.pdb')
                    print('removed solvent from '+i)
                except:
                    print('something went wrong')
            else:
                
                pass
                    
    host = 'pllwskifah1'
    username = 'server'
    port = 22 
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(host,port,username=username) 
    x=1


    if x == 1:

        try:
            for i in os.listdir(lilac_out_dir):

                if i.endswith('.xml.bz2'):

                    try:
                        os.chdir(lilac_out_dir)
                        scp=SCPClient(ssh.get_transport())
                        scp.put(i,remote_path=fah_in_dir+"/.")
                        scp.close()
                        
                    except:
                        print('something went wrong while trying to transport zipped files')
                        
        except:
            print('something went wrong!')
            
        
            
        try:
        	for i in os.listdir(lilac_out_dir):

        		if i.endswith('.pdb'):

        				remove_solvent(lilac_out_dir)

        				os.chdir(lilac_out_dir)

        				scp=SCPClient(ssh.get_transport())
        				scp.put(i,remote_path=fah_in_dir+"/.")
        				scp.close()
                                                

        except:

        	print('Something went wrong with transporting topology files.')  
            
    ssh.close()
            
            
            
            
if __name__ == '__main__':
	send_to_fah(sys.argv[1],sys.argv[2])            
            
            
            
            
            
            
            
            
            
            
