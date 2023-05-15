# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 17:19:59 2020

@author: syang34
"""
# This version improves the efficiency compared to the previous version
import numpy as np
import csv
from scipy.interpolate import griddata
from scipy.spatial import Delaunay
from metpy.interpolate import geometry
from metpy.interpolate.points import natural_neighbor_point
from joblib import Parallel, delayed

#
with open('Output/Mapped_logs.csv', newline='') as f:
    reader = csv.reader(f,delimiter=',')
    mapped_elogs = list(reader)[1::]

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
top_xmax,top_xmin=max(flat_top_surface[:,1]),min(flat_top_surface[:,1])
top_ymax,top_ymin=max(flat_top_surface[:,2]),min(flat_top_surface[:,2])
bot_xmax,bot_xmin=max(flat_bot_surface[:,1]),min(flat_bot_surface[:,1])
bot_ymax,bot_ymin=max(flat_bot_surface[:,2]),min(flat_bot_surface[:,2])
xmax=max([top_xmax,bot_xmax])
xmin=min([top_xmin,bot_xmin])
ymax=max([top_ymax,bot_ymax])
ymin=min([top_ymin,bot_ymin])

# Construct grid data
cellsize = round(np.sum((flat_top_surface[0, 1:3] - flat_top_surface[1, 1:3]) ** 2) ** 0.5, 0)

h_res=cellsize / 5
x_vec=np.linspace(xmin,xmax,int((xmax-xmin)/h_res))
y_vec=np.linspace(ymin,ymax,int((ymax-ymin)/h_res))
z_vec=np.linspace(round(flat_top_surface[0][-1],2),round(flat_bot_surface[0][-1],2),300) # use 300 layers
grid_2d=np.array([[i,j] for i in x_vec for j in y_vec])
#grid_3d=[np.hstack((grid_2d,i*np.ones(len(grid_2d)).reshape((-1,1)))) for i in z_vec]
unstru_grid_2d=[]
for i in grid_2d:
    dis=(flat_top_surface[:,1]-i[0])**2+(flat_top_surface[:,2]-i[1])**2
    if min(dis) < (1.5*cellsize)**2:
        unstru_grid_2d.append(i)
        
unstru_grid_2d=np.array(unstru_grid_2d)
unstru_grid_2d=np.hstack((np.arange(len(unstru_grid_2d)).reshape((-1,1)),unstru_grid_2d)) # ColumnID, X, Y
#np.savetxt('Output/unstruct_grid2d.txt',unstru_grid_2d,fmt='%d %.2f %.2f')

#
# Get point data for each log given elevation of each layer
def get_log_data(organized_elogs,elev):
    log_data=[]
    for i in organized_elogs:
        arr=np.array([j[1:3]+[round(float(j[3]),2)] for j in i]).astype('float')
        elev_seq=arr[:,2]
        soil_seq=[j[4] for j in i]
        if elev<=elev_seq[0] and elev>=elev_seq[-1]:
            ind=np.where(elev_seq==elev)[0]
            if len(ind)>=1:
                log_point=arr[ind[0]]
                soil=soil_seq[ind[0]]
                if soil=='':
                    soil=soil_seq[ind[0]+1]
            else:
                inde=np.where(elev_seq<elev)[0]
                cut=arr[inde[0]]
                soil=soil_seq[inde[0]]
                t=(elev-arr[0][-1])/(cut[-1]-arr[0][-1])
                log_point=arr[0]+(cut-arr[0])*t
            #
            if soil=='S':
                indicator=1
            elif soil=='C':
                indicator=0
            log_data.append([i[0][0]]+log_point.tolist()+[indicator])
    #
    return log_data

# Nearest neighbor interpolation
def nearest_interpo(sample_info,query_info):
    x=np.array(sample_info)[:,1]
    y=np.array(sample_info)[:,2]
    z=np.array(sample_info)[:,4]
    #
    xq=np.array(query_info)[:,1]
    yq=np.array(query_info)[:,2]
    zq=griddata((x,y),z,(xq,yq),method='nearest')
    #
    zq_class=[]
    for i in zq:
        if i>=0.5:
            zq_class.append(1)
        else:
            zq_class.append(0)
    result=[m[0:4]+[n] for m,n in zip(query_info,zq_class)]
    return result

# NN interpolation for one layer
def nn_interpo_layer(m): # m is elevation of the layer
    sample_pt=get_log_data(organized_elogs=organized_elogs,elev=m)
    tri = Delaunay(np.array([p[1:3] for p in sample_pt]))
    members, circumcenters = geometry.find_natural_neighbors(tri, [tuple(p[1::]) for p in unstru_grid_2d.tolist()])
    #
    xp=np.array([p[1] for p in sample_pt])
    yp=np.array([p[2] for p in sample_pt])
    zp=np.array([p[4] for p in sample_pt])
    #
    members_list=[[key,value] for key, value in members.items()]
    #
    nn_layer_results=[]
    nn_nodata=[]
    for i in range(len(unstru_grid_2d)):
        if len(members_list[i][1])>0:
            val = natural_neighbor_point(xp, yp, zp, (unstru_grid_2d[i][1], unstru_grid_2d[i][2]),
                                         tri, members_list[i][1],circumcenters)
            #
            if val>=0.5:
                interpo_class=1
            else:
                interpo_class=0
            #
            nn_layer_results.append([unstru_grid_2d[i][0], unstru_grid_2d[i][1], unstru_grid_2d[i][2], m, interpo_class]) # ColumnID X Y Z indicator
            #
        else:
            nn_nodata.append([unstru_grid_2d[i][0], unstru_grid_2d[i][1], unstru_grid_2d[i][2], m, -2]) # No data area has indicator -2
    #
    nn_layer_results_noextra=nn_layer_results+nn_nodata
    # Nearest neighbor interpolation to handle extrapolation
    nn_extra=nearest_interpo(sample_info=nn_layer_results,query_info=nn_nodata)
    nn_layer_results+=nn_extra
    #
    nn_layer_results=np.array(nn_layer_results)
    nn_layer_results_noextra=np.array(nn_layer_results_noextra)
    nn_layer_results=nn_layer_results[nn_layer_results[:,0].argsort()]
    nn_layer_results_noextra=nn_layer_results_noextra[nn_layer_results_noextra[:,0].argsort()]
    #
    return nn_layer_results.tolist(),nn_layer_results_noextra.tolist()

#
# Main parallel program
nn_results_para=Parallel(n_jobs=-1,verbose=10)(delayed(nn_interpo_layer)(m) for m in z_vec)


# Export
extrapo_result=[p[0] for p in nn_results_para]
organized_extrapo_result=[]
for i in range(0,len(unstru_grid_2d)):
    column=[]
    for j in extrapo_result:
        column.append(j[i])
    organized_extrapo_result.append(column)
    
problem=[]
for column in organized_extrapo_result:
    seq=[p[3] for p in column]
    if all(m > n for m, n in zip(seq, seq[1::])) != 1:
        problem.append(column)

if len(problem)>0:
    print('problem in organize geology in flat unit')

np.savetxt('Output/Interpolation_nondipping_results.txt',np.array([q for p in organized_extrapo_result for q in p]),
           header='ColumnID X Y Z Indicator',fmt='%d %.2f %.2f %.2f %d')
