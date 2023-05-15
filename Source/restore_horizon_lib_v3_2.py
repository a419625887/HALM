# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 22:48:17 2020

@author: syang34
"""

# This module serves as library for horizon restoration
import math as math
import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

# Find adjacent cells for a given cell
def search_adj_cell(pointdata,cellid,searchrange):
    cell_info=pointdata[np.where(pointdata[:,0]==cellid)[0][0]]
    ind=np.where(((pointdata[:,2]-cell_info[2])**2+(pointdata[:,3]-cell_info[3])**2)**0.5<searchrange)[0]
    adj_cell=pointdata[ind]
    return adj_cell.tolist()

# Find leftmost cell given cell center data, and neighbors for all cells
def get_leftmost_cell(pointdata,cellsize):
    cell_adj=[]
    bound_cell=[]
    bound_cell_adj=[]
    for i in pointdata:
        adj=search_adj_cell(pointdata=pointdata,cellid=i[0],searchrange=1.6*cellsize)
        cell_adj.append(adj)
        if len(adj)<8:
            bound_cell.append(i)
            bound_cell_adj.append(np.array(adj))
    #
    left_cell=[]
    for j,k in zip(bound_cell,bound_cell_adj):
        inde=np.where((k[:,2]+0.6*cellsize<j[2]) & (k[:,3]-0.2*cellsize<j[3]) & (k[:,3]+0.2*cellsize>j[3]))[0]
        if len(inde)==0:
            left_cell.append(j.tolist())
    #
    return np.array(left_cell)[np.array(left_cell)[:,3].argsort()].tolist(),cell_adj

# Find cell ceneter in a row given a leftmost cell
def get_row_cell(cell_adj,leftcell,cellsize):
    row_cell=[leftcell]
    stop_search=0
    t=0
    #
    while stop_search==0:
        neighbor=np.array(cell_adj[int(row_cell[t][0])])
        ind=np.where((neighbor[:,2]-0.6*cellsize>row_cell[t][2]) & (neighbor[:,3]-0.2*cellsize<row_cell[t][3]) & (neighbor[:,3]+0.2*cellsize>row_cell[t][3]))[0]
        if len(ind)!=0:
            next_cell=neighbor[ind[0]].tolist()
            row_cell.append(next_cell)
        else:
            stop_search=1
        #
        t+=1
    #
    return row_cell

# Find mesh faces in a row given two row cell nodes
def get_row_face(row_node1,row_node2,face_center,face_vertex,vertex_id):
    row_face_center=[]
    row_face_pt=[]
    for i,j,k in zip(face_center,face_vertex,vertex_id):
        point = Point(i[0], i[1])
        polygon = Polygon([tuple(row_node1[0][2:4]),tuple(row_node1[-1][2:4]),
                           tuple(row_node2[-1][2:4]),tuple(row_node2[0][2:4])])
        judge=polygon.contains(point)
        #
        if judge==1:
            row_face_center.append(i.tolist())
            row_face_pt.append([[m]+n for m,n in zip(k[1::],j)])
    #
    arr=np.array(row_face_center)
    #
    return np.array(row_face_pt)[arr[:,0].argsort()].tolist()

# Rotate point given an origin
def rotate_point(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin[0:2]
    px, py = point[0:2]

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return [qx, qy, point[2]]

# Rigid translation to fit triangles with two shared nodes
def rigid_fit(datapool,shared_pt,unshared_pt):
    p1=np.array(shared_pt[0])
    p2=np.array(shared_pt[1])
    p3=np.array(unshared_pt[0])
    arr=np.array(datapool)
    dp1=arr[np.where(arr[:,0]==p1[0])[0][0]]
    dp2=arr[np.where(arr[:,0]==p2[0])[0][0]]
    #
    trans=dp1[1::]-p1[1::]
    p1_trans=p1[1::]+trans
    p2_trans=p2[1::]+trans
    #
    v1=dp2[1::]-dp1[1::]
    v2=p2_trans-p1_trans
    cos=np.dot(v1,v2)/(np.sqrt(sum(v1**2))*np.sqrt(sum(v2**2)))
    angle=math.acos(cos) #angle in radians
    #
    p2_unclock=rotate_point(origin=p1_trans,point=p2_trans,angle=angle)
    p3_unclock=rotate_point(origin=p1_trans,point=p3[1::]+trans,angle=angle)
    check1=np.dot(v1,p2_unclock-p1_trans)/(np.sqrt(sum(v1**2))*np.sqrt(sum((p2_unclock-p1_trans)**2)))
    p2_clock=rotate_point(origin=p1_trans,point=p2_trans,angle=-angle)
    p3_clock=rotate_point(origin=p1_trans,point=p3[1::]+trans,angle=-angle)
    check2=np.dot(v1,p2_clock-p1_trans)/(np.sqrt(sum(v1**2))*np.sqrt(sum((p2_clock-p1_trans)**2)))
    if abs(check1-1)<abs(check2-1):
        p2_fit=p2_unclock
        p3_fit=p3_unclock
    else:
        p2_fit=p2_clock
        p3_fit=p3_clock
    #
    if np.sqrt(sum((np.array(p2_fit)-dp2[1::])**2))<0.001:
        fit_results=[[p1[0]]+p1_trans.tolist(),[p2[0]]+p2_fit,[p3[0]]+p3_fit]
        distort=0
    else:
        #trans2=dp2[1::]-np.array(p2_fit)
        #p2_fit=np.array(p2_fit)+trans2
        #p3_fit=np.array(p3_fit)+trans2
        #fit_results=[[p1[0]]+p1_trans.tolist(),[p2[0]]+p2_fit.tolist(),[p3[0]]+p3_fit.tolist()]
        fit_results=awk_fit(datapool,shared_pt,unshared_pt)
        distort=1
    return fit_results,distort

# Fit a triangle fixed by 3 points in datapool
def three_pt_fit(datapool,shared_pt): # Should have 3 shared points
    arr=np.array(datapool)
    fit_results=[]
    for i in np.array(shared_pt):
        dp=arr[np.where(arr[:,0]==i[0])[0][0]]
        fit_results.append(dp.tolist())
    return fit_results

# Compute triangle area given three vertexs
def tri_area(tri_list):
    verte=np.array(tri_list)[:,1::]
    side1=verte[0]-verte[1]
    side2=verte[0]-verte[2]
    cos_angle=np.dot(side1,side2)/(np.sqrt(sum(side1**2))*np.sqrt(sum(side2**2)))
    area=0.5*np.sqrt(sum(side1**2))*np.sqrt(sum(side2**2))*(1-cos_angle**2)**0.5
    return area

# Strain fit for the triangle following a triangle fixed by 3 nodes in datapool
def strain_fit(datapool,shared_pt,unshared_pt,previous_fit_tri,previous_unfit_tri):
    p1=np.array(shared_pt[0])
    p2=np.array(shared_pt[1])
    p3=np.array(unshared_pt[0])
    arr=np.array(datapool)
    dp1=arr[np.where(arr[:,0]==p1[0])[0][0]]
    dp2=arr[np.where(arr[:,0]==p2[0])[0][0]]
    #
    trans1=dp1[1::]-p1[1::]
    p1_trans=p1[1::]+trans1
    p2_med=p2[1::]+trans1
    trans2=dp2[1::]-p2_med
    p2_trans=p2_med+trans2
    p3_trans=p3[1::]+trans1+trans2
    #
    tri1_fit_area=tri_area(tri_list=previous_fit_tri) # Area of a previous fitted triangle
    tri1_unfit_area=tri_area(tri_list=previous_unfit_tri) # Area of a previous triangle before fitting
    tri2_unfit_area=tri_area(tri_list=shared_pt+unshared_pt) # Area of the focused triangle before fitting
    #
    t=-np.dot(p1_trans-p3_trans,p2_trans-p1_trans)/sum((p2_trans-p1_trans)**2)
    v=p1_trans+(p2_trans-p1_trans)*t # Projection of p3_trans on the line passing p1_trans_and p2_trans
    d=2*(tri1_unfit_area+tri2_unfit_area-tri1_fit_area)/np.sqrt(sum((p2_trans-p1_trans)**2)) # Distance v and p3_fit
    D=np.sqrt(sum((p3_trans-v)**2)) # Distance between v and p3_trans
    p3_fit=v+(p3_trans-v)*(d/D)
    #
    fit_results=[[p1[0]]+p1_trans.tolist(),[p2[0]]+p2_trans.tolist(),[p3[0]]+p3_fit.tolist()]
    return fit_results

# Fit face without three shared points following a strain-fitted face
def awk_fit(datapool,shared_pt,unshared_pt):
    p1=np.array(shared_pt[0])
    p2=np.array(shared_pt[1])
    p3=np.array(unshared_pt[0])
    arr=np.array(datapool)
    dp1=arr[np.where(arr[:,0]==p1[0])[0][0]]
    dp2=arr[np.where(arr[:,0]==p2[0])[0][0]]
    #
    trans1=dp1[1::]-p1[1::]
    p1_trans=p1[1::]+trans1
    p2_med=p2[1::]+trans1
    trans2=dp2[1::]-p2_med
    p2_trans=p2_med+trans2
    p3_trans=p3[1::]+trans1+trans2
    #
    fit_results=[[p1[0]]+p1_trans.tolist(),[p2[0]]+p2_trans.tolist(),[p3[0]]+p3_trans.tolist()]
    return fit_results

# Fit out_face using rigid fit method
def fit_out_face(one_face,id_data,pt_data):
    share_pt=[]
    unshared_pt=[]
    for k in one_face:
        if k[0] in id_data:
            share_pt.append(k)
        else:
            unshared_pt.append(k)
    #
    fit_face=rigid_fit(datapool=pt_data,shared_pt=share_pt,unshared_pt=unshared_pt)[0]
    #
    return fit_face
        
# Update datapool given a fitted face
def update_datapool(id_data,pt_data,fit_face):
    for p in fit_face:
        if p[0] not in id_data:
            id_data.append(p[0])
            pt_data.append(p)
    #
    return id_data,pt_data