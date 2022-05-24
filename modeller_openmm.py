#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 22:02:26 2022

@author: belayv
"""


import openmm as mm
from openmm.app import CharmmPsfFile
from openmm.app import PDBFile
import modeller as md

psf = '/Users/belayv/Downloads/charmm-gui-5041720964/openmm/step3_input_mod.psf' # Define psf location
pdb = '/Users/belayv/Downloads/charmm-gui-5041720964/openmm/step3_input.pdb' # Define pdb location


topo= CharmmPsfFile(psf) # Load openmm topology from psf
topology= topo.topology # Load openmm positions from pdb

poso = PDBFile(pdb)
# topology= poso.topology
atoms = [atom for atom in topology.atoms()]

toDelete = [atoms for atoms in topology.atoms() if (atoms.residue.id == '2562' or atoms.residue.id == '2565') and atoms.name =='HG']

mod = mm.app.Modeller(topology,poso.positions)

mod.delete(toDelete)

new_topology=mod.getTopology()
new_positions=mod.getPositions()

atom2 = [atom2 for atom2 in new_topology.atoms()]