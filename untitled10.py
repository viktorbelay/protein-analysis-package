#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 22:31:39 2022

@author: belayv
"""

import openmm as mm
from openmm.app import CharmmPsfFile
from openmm.app import PDBFile
import modeller as md

psf = '/Users/belayv/Downloads/charmm-gui-5041720964/openmm/step3_input.psf'
pdb = '/Users/belayv/Downloads/charmm-gui-5041720964/openmm/step3_input.pdb'


topo= CharmmPsfFile(psf)
topology= topo.topology

poso = PDBFile(pdb)
# topology= poso.topology
atoms = [atom for atom in topology.atoms()]

toDelete = [r for r in topology.atoms() if r.name == 'H']

mod = mm.app.Modeller(topology,poso.positions)

mod.delete(toDelete)

new_topology=mod.getTopology()

atom2 = [atom2 for atom2 in new_topology.atoms()]

