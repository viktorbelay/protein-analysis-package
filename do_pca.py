#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 16:27:24 2022

@author: belayv
"""
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import pca, align
import seaborn as sns
from matplotlib import pyplot as pl
import sys

def do_pca(topology_path,topology_name,trajectory_path,trajectory_name):
    
    pdb_path=topology_path+topology_name
    dcd_path=trajectory_path+trajectory_name
    
    
    u = mda.Universe(pdb_path,dcd_path)

    pc = pca.PCA(u,select='backbone',align=True,mean=None,n_components=None).run()



    backbone = u.select_atoms('backbone')



    transformed = pc.transform(backbone, n_components=2)
    transformed.shape


    df = pd.DataFrame(transformed,
                      columns=['PC{}'.format(i+1) for i in range(2)])
    df['Time (ps)'] = df.index * u.trajectory.dt * 21.5
    df.head()

    pl.scatter(df['PC1'],df['PC2'],c=df['Time (ps)'])
    pl.colorbar().set_label('time (ns)')
    pl.xlabel('PC1')
    pl.ylabel('PC2')
    
    pl.title('Trajectory PCA, first two PCs')

    
    pl.savefig(trajectory_path+'pca_analysis_two_pcs.png', dpi=300)
    pl.close()
    
    
if __name__=='__main__':
    do_pca(sys.argv[1],sys.argv[2])