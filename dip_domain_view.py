# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 09:17:18 2022

@author: syang
"""


import numpy as np
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

#
logs = np.loadtxt('Output/log_intersection.txt', skiprows = 1, delimiter = ',')
logloc = logs[:, 0:2].tolist()
topbotpt = [[p[2], p[5]] for p in logs]

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

logseqs=[]
for i in range(len(organized_elogs)):
    seq=[[topbotpt[i][0], '']]
    for j in organized_elogs[i][1::]:
        r=(organized_elogs[i][0][3]-j[3])/(organized_elogs[i][0][3]-organized_elogs[i][-1][3])
        seq.append([topbotpt[i][0]+r*(topbotpt[i][1]-topbotpt[i][0]), j[4]])
    logseqs.append(seq)

#
top_surface=np.loadtxt('Input/top_surface_point.txt',skiprows=1)
bot_surface=np.loadtxt('Input/bot_surface_point.txt',skiprows=1)

flat_top_surface=np.loadtxt('Output/restored_top_point_data.txt',skiprows=1)
flat_top_surface=flat_top_surface[flat_top_surface[:,0].argsort()]
flat_top_surface[:,-1]=300. #Set flat top surface elevation_meter
flat_bot_surface=np.loadtxt('Output/restored_bot_point_data.txt',skiprows=1)
flat_bot_surface=flat_bot_surface[flat_bot_surface[:,0].argsort()]
flat_bot_surface[:,-1]=10. #Set flat bottom surface elevation_meter

faces=np.loadtxt('Input/triangle_faces.txt',skiprows=1).astype('int32')
face_pt_top=[]
for i in faces:
    face_pt_top.append(top_surface[i[1::], 1::].tolist())
face_pt_bot=[]
for i in faces:
    face_pt_bot.append(bot_surface[i[1::], 1::].tolist())

#
fig1 = plt.figure(figsize=(8,8))
ax = plt.axes(projection='3d')
ax.set_xlim(-100, 3300)
ax.set_ylim(-100, 1300)
ax.set_zlim(-80, 80)
ax.set_xlabel('X', fontsize=14)
ax.set_ylabel('Y', fontsize=14)
ax.set_zlabel('Z', fontsize=14)
ax.tick_params(axis='x',labelsize=14)
ax.tick_params(axis='y',labelsize=14)
ax.tick_params(axis='z',labelsize=14)
for i,j in zip(logloc,logseqs):
    ax.plot([i[0], i[0]], [i[1], i[1]], [60., j[0][0]], color='green', linewidth=5)
    for k in range(1, len(j)):
        if j[k][1] == 'S':
            cl='yellow'
        else:
            cl='blue'
        ax.plot([i[0], i[0]], [i[1], i[1]], [j[k-1][0], j[k][0]], color=cl, linewidth=5)
    ax.plot([i[0], i[0]], [i[1], i[1]], [j[-1][0], -60.], color='green', linewidth=5)    
collection = Poly3DCollection(face_pt_top, linewidths=0.7, alpha=0.2)
collection.set_edgecolor('k')
ax.add_collection3d(collection)
collection = Poly3DCollection(face_pt_bot, linewidths=1, alpha=0.2)
collection.set_edgecolor('c')
ax.add_collection3d(collection)
fig1.savefig('Plots/Dip_domain_view.png', dpi=500)