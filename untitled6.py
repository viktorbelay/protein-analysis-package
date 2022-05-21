#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 22:35:49 2022

@author: belayv

# Purpose of script: Find residue indices with OpenMM's indexing system
"""

import openmm as mm
from openmm.app import CharmmPsfFile, PDBFile, PDBxFile, CharmmParameterSet


psf = '/Users/belayv/Downloads/apotest6.pdb'

psf_openmm = PDBFile(psf)

topology = psf_openmm.topology

atoms = [atom for atom in topology.atoms()]

for i in list(range(0,len(atoms))):
    
    if (atoms[i].residue.id == ('1909') or 
        atoms[i].residue.id == ('1910') or 
        atoms[i].residue.id == ('1911') or 
        atoms[i].residue.id == ('2214') or 
        atoms[i].residue.id == ('2215') or 
        atoms[i].residue.id == ('2216') or 
        atoms[i].residue.id == ('2549') or
        atoms[i].residue.id == ('2550') or
        atoms[i].residue.id == ('2551') or
        atoms[i].residue.id == ('2634') or
        atoms[i].residue.id == ('2635') or
        atoms[i].residue.id == ('2636')
        ) and (atoms[i].residue.chain == atoms[0].residue.chain or
               
               atoms[i].residue.chain == atoms[5000].residue.chain
               
               ) and (atoms[i].name in ('CA','C','N')
                   
                   
                   
                   ):
        
        print(atoms[i])