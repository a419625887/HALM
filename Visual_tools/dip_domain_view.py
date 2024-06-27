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
with open('Input/Logs.csv', newline='') as f:
    reader = csv.reader(f,delimiter=',')
    logdata = list(reader)

for count in range(len(logdata)):
    if '' in logdata[count]:
        logdata[count]=logdata[count][0:logdata[count].index('')]


#
top_surface=np.loadtxt('Input/top_surface_point.txt',skiprows=1)
bot_surface=np.loadtxt('Input/bot_surface_point.txt',skiprows=1)


faces=np.loadtxt('Input/triangle_faces.txt',skiprows=1).astype('int32')
face_pt_top=[]
for i in faces:
    face_pt_top.append(top_surface[i[1::], :].tolist())
face_pt_bot=[]
for i in faces:
    face_pt_bot.append(bot_surface[i[1::], :].tolist())

#
fig1 = plt.figure(figsize=(8,8))
ax = plt.axes(projection='3d')
xmin, xmax = np.min(top_surface[:, 0]), np.max(top_surface[:, 0])
ymin, ymax = np.min(top_surface[:, 1]), np.max(top_surface[:, 1])
zmin, zmax = np.min(bot_surface[:, 2]), np.max(top_surface[:, 2])
ax.set_xlim(xmin - 50, xmax + 50)
ax.set_ylim(ymin - 50, ymax + 50)
ax.set_zlim(zmin - 50, zmax + 50)
ax.set_xlabel('X', fontsize=14)
ax.set_ylabel('Y', fontsize=14)
ax.set_zlabel('Z', fontsize=14)
ax.tick_params(axis='x',labelsize=14)
ax.tick_params(axis='y',labelsize=14)
ax.tick_params(axis='z',labelsize=14)

for i in range(0, len(logdata), 2):
    logx = float(logdata[i+1][1])
    logy = float(logdata[i+1][2])
    datum = float(logdata[i+1][3])
    dep = np.array([p for p in logdata[i][3::]]).astype('float64')
    elev = (datum - dep) * 0.3048
    litho = logdata[i+1][4::]
    #
    for j in range(len(elev) - 1):
        if litho[j] == 'S':
            cl = 'yellow'
        else:
            cl = 'blue'
        #
        ax.plot([logx, logx], [logy, logy], [elev[j], elev[j+1]], color=cl, linewidth=5)

    
collection = Poly3DCollection(face_pt_top, linewidths=0.7, alpha=0.2)
collection.set_edgecolor('k')
ax.add_collection3d(collection)
collection = Poly3DCollection(face_pt_bot, linewidths=1, alpha=0.2)
collection.set_edgecolor('c')
ax.add_collection3d(collection)

fig1.savefig('Plots/Dip_domain_view.png', dpi=500)