# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 09:31:22 2022

@author: syang
"""


import numpy as np
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

#   
with open('Output/Mapped_logs.csv', newline='') as f:
    reader = csv.reader(f,delimiter=',')
    mapped_elogs = list(reader)[1::]
# Organized elog data
log_name=[]
for i in mapped_elogs:
    if i[0] not in log_name:
        log_name.append(i[0])

organized_elogs=[]
for i in log_name:
    one_log=[]
    for j in mapped_elogs:
        if j[0]==i:
            one_log.append(j[0:3]+[round(float(j[3]),2)]+[j[4]])
    one_log[0][3]=round(one_log[0][3],0)
    one_log[-1][3]=round(one_log[-1][3],0)
    organized_elogs.append(one_log)

#
flat_top_surface=np.loadtxt('Output/restored_top_point_data.txt',skiprows=1)
flat_top_surface=flat_top_surface[flat_top_surface[:,0].argsort()]
flat_top_surface[:,-1]=300. #Set flat top surface elevation_meter
flat_bot_surface=np.loadtxt('Output/restored_bot_point_data.txt',skiprows=1)
flat_bot_surface=flat_bot_surface[flat_bot_surface[:,0].argsort()]
flat_bot_surface[:,-1]=10. #Set flat bottom surface elevation_meter

faces=np.loadtxt('Input/triangle_faces.txt',skiprows=1).astype('int32')
face_pt_top=[]
for i in faces:
    face_pt_top.append(flat_top_surface[i[1::], 1::].tolist())
face_pt_bot=[]
for i in faces:
    face_pt_bot.append(flat_bot_surface[i[1::], 1::].tolist())

#
# 3D view of horizons
fig1 = plt.figure(figsize=(8,8))
ax = plt.axes(projection='3d')
xmin, xmax = np.min(flat_bot_surface[:, 1]), np.max(flat_bot_surface[:, 1])
ymin, ymax = np.min(flat_bot_surface[:, 2]), np.max(flat_bot_surface[:, 2])
ax.set_xlim(xmin - 50, xmax + 50)
ax.set_ylim(ymin - 50, ymax + 50)
ax.set_zlim(-80, 380)
ax.set_xlabel('X', fontsize=14)
ax.set_ylabel('Y', fontsize=14)
ax.set_zlabel('Z', fontsize=14)
ax.tick_params(axis='x',labelsize=14)
ax.tick_params(axis='y',labelsize=14)
ax.tick_params(axis='z',labelsize=14)  
for i in organized_elogs:
    for j in range(1, len(i)):
        if i[j][4] == 'S':
            cl='yellow'
        else:
            cl='blue'
        ax.plot([float(i[j-1][1]), float(i[j][1])], [float(i[j-1][2]), float(i[j][2])], [float(i[j-1][3]), float(i[j][3])],
                color=cl, linewidth=5)
collection = Poly3DCollection(face_pt_top, linewidths=0.7, alpha=0.2)
collection.set_edgecolor('k')
ax.add_collection3d(collection)
collection = Poly3DCollection(face_pt_bot, linewidths=1, alpha=0.2)
collection.set_edgecolor('c')
ax.add_collection3d(collection)
fig1.savefig('Plots/Nondip_domain_view.png', dpi=500)