#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 12:28:01 2022

@author: belayv
"""

import mdtraj as md
import paramiko
from scp import SCPClient,SCPException
import os
import sys
<<<<<<< HEAD
=======
import openmm as mm
>>>>>>> 29a7ba72914c11d1f03d3e09f2d5429d0498fd7f


def send_to_fah(lilac_out_dir,fah_in_dir):
    
    def remove_solvent(lilac_out_dir):
        
        for i in os.listdir(lilac_out_dir):
            
            if i.endswith('.pdb'):
                
                os.chdir(lilac_out_dir)
                print('removing solvent from '+i)
                md.load(i).remove_solvent().save_pdb(lilac_out_dir+'/nosolvent_equilibrated')
                print('removed solvent from '+i)
                
<<<<<<< HEAD
            elif i.endswith('.cif'):
                
                os.chdir(lilac_out_dir)
                
                
                md.load(i).save_pdb(lilac_out_dir+'/equiliberated.pdb')
                
                md.load(i).remove_solvent().save_pdb(lilac_out_dir+'nosolvent_equiliberated.pdb')
                
                print('solvent stripped from cif.')
=======
            # elif i.endswith('.cif'):
                
            #     os.chdir(lilac_out_dir)
                
                
            #     md.load(i).save_pdb(lilac_out_dir+'/equiliberated.pdb')
                
            #     md.load(i).remove_solvent().save_pdb(lilac_out_dir+'nosolvent_equiliberated.pdb')
                
            #     print('solvent stripped from cif.')
>>>>>>> 29a7ba72914c11d1f03d3e09f2d5429d0498fd7f
                
            else:
                pass
            
            
    host=  'pllwskifah1'
    username='server'
    port = 22
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(host,port,username=username)
    x=1
    
    if x == 1:
        
        for i in os.listdir(lilac_out_dir):
            
            if i.endswith('.xml.bz2'):
                
                os.chdir(lilac_out_dir)
                scp=SCPClient(ssh.get_transport())
                scp.put(i,remote_path=fah_in_dir+'/.')
                scp.close()
                
        remove_solvent(lilac_out_dir)
        
        for i in os.listdir(lilac_out_dir):
            
            if i.endswith('.pdb'):
                
                os.chdir(lilac_out_dir)
                scp=SCPClient(ssh.get_transport())
                scp.put(i,remote_path=fah_in_dir+'/.')
                scp.close()
                
    ssh.close()
    
    
if __name__=='__main__':
    
<<<<<<< HEAD
    send_to_fah(sys.argv[1],sys.argv[2])
=======
    send_to_fah(sys.argv[1],sys.argv[2])
>>>>>>> 29a7ba72914c11d1f03d3e09f2d5429d0498fd7f
