# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 10:02:25 2022

@author: syang
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


#
allpt=np.loadtxt('Output/Interpolation_nondipping_results.txt', skiprows=1)

# make color map
colors = [(0, 0, 1), (1, 1, 0)]  # R -> G -> B
cmap_name = 'my_list'
cmap = LinearSegmentedColormap.from_list(cmap_name, colors)
#
fig1 = plt.figure(figsize=(8,8))
ax = plt.axes(projection='3d')
ax.set_xlim(-100, 3300)
ax.set_ylim(-100, 1300)
ax.set_zlim(-80, 380)
ax.set_xlabel('X', fontsize=14)
ax.set_ylabel('Y', fontsize=14)
ax.set_zlabel('Z', fontsize=14)
ax.tick_params(axis='x',labelsize=14)
ax.tick_params(axis='y',labelsize=14)
ax.tick_params(axis='z',labelsize=14)
v3d=ax.scatter3D(allpt[:,1],allpt[:,2],allpt[:,3],c=allpt[:,4],cmap=cmap,s=10)
fig1.savefig('Plots/Nondip_lithology_view.png', dpi=500)
