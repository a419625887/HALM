# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 10:02:25 2022

@author: syang
"""


import numpy as np
import matplotlib.pyplot as plt
import voxel_plot as ve

#
model = np.loadtxt('Output/Geo-model_nondip.dat')

xx = np.unique(model[:, 1])
res = xx[1] - xx[0] # horizontal resolution of grid

# Create cells
colid=np.unique(model[:,0])
organized_model=[]
positions=[]
heights=[]
colors=[]
for i in colid:
    ind=np.where(model[:,0]==i)[0]
    col=model[ind]
    for j in col[1::]:
        positions.append([col[0][1]-res/2.0, col[0][2]-res/2.0, j[3]])
    heights += (col[0:-1, 3]-col[1::, 3]).tolist()
    for j in col[0:-1]:
        if j[4] == 1:
            colors.append('yellow')
        else:
            colors.append('blue')
    #
    organized_model.append(col.tolist())

sizes=[(res,res,p) for p in heights]
positions=np.array(positions)


#
fig1 = plt.figure(figsize=(8,8))
ax = plt.axes(projection='3d')
pc = ve.plotCubeAt(positions, sizes=sizes, colors=colors,edgecolor="none")
ax.add_collection3d(pc)
xmin, xmax = np.min(model[:, 1]), np.max(model[:, 1])
ymin, ymax = np.min(model[:, 2]), np.max(model[:, 2])
ax.set_xlim(xmin - 50, xmax + 50)
ax.set_ylim(ymin - 50, ymax + 50)
ax.set_zlim(-80, 380)
ax.set_xlabel('X', fontsize=14)
ax.set_ylabel('Y', fontsize=14)
ax.set_zlabel('Z', fontsize=14)
ax.tick_params(axis='x',labelsize=14)
ax.tick_params(axis='y',labelsize=14)
ax.tick_params(axis='z',labelsize=14)

fig1.savefig('Plots/Nondip_lithology_view.png', dpi=500)
