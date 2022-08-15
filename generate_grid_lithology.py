# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 14:34:03 2020

@author: syang34
"""

import numpy as np
from joblib import Parallel, delayed
import os

# Import and process data
mapped_geo=np.loadtxt('Output/Interpolation_results_map_back.txt',skiprows=1)

horizon_cell_top=np.loadtxt('Input/domain_discretization.txt',skiprows=1, usecols = [0, 1, 2, 3])
horizon_cell_bot=np.loadtxt('Input/domain_discretization.txt',skiprows=1, usecols = [0, 1, 2, 4])

layerid=np.unique(mapped_geo[:,0])
layer_num=len(layerid)

#
# Descritilize column given top, bot, and resolution
def descritilize_col(top,bot,resol):
    col=[top.tolist()]
    t=0
    while round(col[t][-1]-resol,2) >= bot[-1]:
        next_pt=top[0:3].tolist()+[round(col[t][-1]-resol,2)]
        col.append(next_pt)
        t+=1
    #
    if col[-1][-1] > bot[-1]:
        col.append(bot.tolist())
    #
    return col

# upscale a column of geology
def upscale_column(col_geo):
    up_col=[col_geo[0]]
    for p in range(len(col_geo)-1):
        if col_geo[p][-1] != col_geo[p+1][-1]:
            up_col.append(col_geo[p+1])
    #
    if up_col[-1][3] > col_geo[-1][3]:
        up_col.append(col_geo[-1][0:4]+[-999])
    elif up_col[-1][3] == col_geo[-1][3]:
        up_col[-1][-1]=-999
    #
    return up_col

# Get column geology
def column_geology(top_data,bot_data):
#for top_data,bot_data in zip(horizon_cell_top,horizon_cell_bot):
    des_col=descritilize_col(top=top_data, bot=bot_data, resol=0.1)
    layered_col=np.linspace(des_col[0][3],des_col[-1][3],layer_num) # divide col according to layers of mapped geology
    for w in range(len(des_col)):
        v_dis=np.absolute(des_col[w][-1]-layered_col)
        tind=np.where(v_dis==min(v_dis))[0]
        des_col[w].append(tind[0]) # append target layer id on des_col
    #
    col_xy=[des_col[0][1],des_col[0][2]]
    h_dis=(mapped_geo[:,1]-col_xy[0])**2+(mapped_geo[:,2]-col_xy[1])**2
    h_dis_min=min(h_dis)
    if h_dis_min > (1.5*50)**2: # No data area is 1.5km away from the mapped back geology cloud
        col_geo=[m[0:4]+[-2] for m in des_col]
    else:
        indl=np.where(h_dis < (3*50)**2)[0]
        mapped_geo_local=mapped_geo[indl]
        col_geo=[]
        for k in des_col:
            indx=np.where(mapped_geo_local[:,0]==k[-1])[0]
            layergeo=mapped_geo_local[indx]
            #
            dist2_arr=(layergeo[:,1]-k[1])**2+(layergeo[:,2]-k[2])**2 # 2D search
            #dist2_arr=(layergeo[:,1]-k[1])**2+(layergeo[:,2]-k[2])**2+(layergeo[:,3]-k[3])**2 # 3D search
            ind=np.where(dist2_arr==min(dist2_arr))[0]
            soil=layergeo[ind[0]][-1]
            col_geo.append(k[0:4]+[soil])
    #
    col_geo_upscale=upscale_column(col_geo=col_geo)
    #
    return col_geo_upscale
    

#
grid_geology_para=Parallel(n_jobs=-1,verbose=1)(delayed(column_geology)(i,j) for i,j in zip(horizon_cell_top,horizon_cell_bot))

#
# Export
np.savetxt('Output/grid_lithology.txt',np.array([q for p in grid_geology_para for q in p]),
           header='2D_CellID X Y Z_meter Indicator',fmt='%d %.2f %.2f %.2f %d')
os.remove('Output/Interpolation_results_map_back.txt')
