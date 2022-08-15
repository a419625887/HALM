# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 11:32:18 2020

@author: syang34
"""

import pyvista as pv
import numpy as np

#
xyz = np.loadtxt('Input/top_surface_point.txt')
faces = np.loadtxt('Input/triangle_faces.txt').astype('int32')
vertices=xyz[:,1::]

#
surf = pv.PolyData(vertices, faces)
surf.compute_normals(cell_normals=True, point_normals=False, inplace=True)

# Extract points for faces
face=surf.faces.reshape(-1, 4)
points=surf.points
face_pt=[]
for i in face:
    face_pt.append(points[i[1::]].tolist())

# Find face normal
norm=surf['Normals']

# Find face center
centers=surf.cell_centers().points
#
#
def normalize_v2(vec):
    ''' Normalize a numpy array of 2 component vectors '''
    lens = np.sqrt( vec[0]**2 + vec[1]**2 )
    vec[0] /= lens
    vec[1] /= lens               
    return np.array(vec.tolist()+[0])
#
#Rotate faces
rotated_face=[]
for i,j,centroid in zip(norm,face_pt,centers):
    if i[2]>0:
        n=i
    else:
        n=-i
    x=(n[2]**2/(n[1]**2+n[2]**2))**0.5
    d=np.array([x*n[0],x*n[1],-(n[0]**2+n[1]**2)**0.5])
    s=np.cross(n,d)
    #
    rotated_pt=[]
    for p in np.array(j):
        t=-np.dot(centroid-p,s)/sum(s**2)
        v=centroid+s*t
        dist=sum((p-v)**2)**0.5
        dire=normalize_v2((p-v)[0:2])
        node=v+dire*dist
        rotated_pt.append(node.tolist())
    #
    rotated_face.append(rotated_pt)

# Find triangle faces from left to right in each row
import restore_horizon_lib_v3_2 as rh
xyz=np.hstack((np.arange(len(xyz)).reshape((-1,1)),xyz))
left_cell,cell_neighbor=rh.get_leftmost_cell(pointdata=xyz,cellsize=50)
row_face_coll=[]
for i in range(len(left_cell)-1):
    row1=rh.get_row_cell(cell_adj=cell_neighbor,leftcell=left_cell[i],cellsize=50)
    row2=rh.get_row_cell(cell_adj=cell_neighbor,leftcell=left_cell[i+1],cellsize=50)
    row_face_pt=rh.get_row_face(row_node1=row1,row_node2=row2,
                                face_center=centers,face_vertex=rotated_face,vertex_id=face)
    row_face_coll.append(row_face_pt)

# Fit triangle faces one by one for the first row
strtrow=0
id_data=[i[0] for i in row_face_coll[strtrow][0]]
pt_data=row_face_coll[strtrow][0].copy()
fit_face_coll=[]
for i in row_face_coll[strtrow]:
    share_pt=[]
    unshared_pt=[]
    for j in i:
        if j[0] in id_data:
            share_pt.append(j)
        else:
            unshared_pt.append(j)
    #
    if len(share_pt)==3:
        fit_face=share_pt
    elif len(share_pt)==2:
        fit_face=rh.rigid_fit(datapool=pt_data,shared_pt=share_pt,unshared_pt=unshared_pt)[0]
    fit_face_coll.append(fit_face)
    #
    for p in fit_face:
        if p[0] not in id_data:
            id_data.append(p[0])
            pt_data.append(p)

# Fit all remaining triangles
for i in row_face_coll[strtrow+1::]:
    indi=0 # whether the previous fit is a three-point fit or an awkward fit
    #istrainfit=0 # whether the previous fit is a strain fit
    out_face=[]
    for j in i:
        share_pt=[]
        unshared_pt=[]
        for k in j:
            if k[0] in id_data:
                share_pt.append(k)
            else:
                unshared_pt.append(k)
        #
        if len(share_pt)==2 and indi==0:
            fit_face,indi=rh.rigid_fit(datapool=pt_data,shared_pt=share_pt,unshared_pt=unshared_pt)
            face_bef_fit=share_pt+unshared_pt
            #istrainfit=0
            fit_face_coll.append(fit_face)
        elif len(share_pt)==2 and indi==1:
            fit_face=rh.strain_fit(datapool=pt_data,shared_pt=share_pt,unshared_pt=unshared_pt,
                                   previous_fit_tri=fit_face_coll[-1],previous_unfit_tri=face_bef_fit)
            indi=0
            #istrainfit=1
            fit_face_coll.append(fit_face)
        elif len(share_pt)==3:
            fit_face=rh.three_pt_fit(datapool=pt_data,shared_pt=share_pt)
            face_bef_fit=share_pt
            indi=1
            #istrainfit=0
            fit_face_coll.append(fit_face)
        else:
            out_face.append(j)
        # update datapool
        id_data,pt_data=rh.update_datapool(id_data=id_data,pt_data=pt_data,fit_face=fit_face)
    # Fit out_face
    if len(out_face)>0:
        for w in out_face[::-1]:
            fit_face=rh.fit_out_face(one_face=w,id_data=id_data,pt_data=pt_data)
            fit_face_coll.append(fit_face)
            id_data,pt_data=rh.update_datapool(id_data=id_data,pt_data=pt_data,fit_face=fit_face)
            
for i in row_face_coll[0:strtrow][::-1]:
    indi=0 # whether the previous fit is a three-point fit or an awkward fit
    #istrainfit=0 # whether the previous fit is a strain fit
    out_face=[]
    for j in i[::-1]:
        share_pt=[]
        unshared_pt=[]
        for k in j:
            if k[0] in id_data:
                share_pt.append(k)
            else:
                unshared_pt.append(k)
        #
        if len(share_pt)==2 and indi==0:
            fit_face,indi=rh.rigid_fit(datapool=pt_data,shared_pt=share_pt,unshared_pt=unshared_pt)
            face_bef_fit=share_pt+unshared_pt
            #istrainfit=0
            fit_face_coll.append(fit_face)
        elif len(share_pt)==2 and indi==1:
            fit_face=rh.strain_fit(datapool=pt_data,shared_pt=share_pt,unshared_pt=unshared_pt,
                                   previous_fit_tri=fit_face_coll[-1],previous_unfit_tri=face_bef_fit)
            indi=0
            #istrainfit=1
            fit_face_coll.append(fit_face)
        elif len(share_pt)==3:
            fit_face=rh.three_pt_fit(datapool=pt_data,shared_pt=share_pt)
            face_bef_fit=share_pt
            indi=1
            #istrainfit=0
            fit_face_coll.append(fit_face)
        else:
            out_face.append(j)
        # update datapool
        id_data,pt_data=rh.update_datapool(id_data=id_data,pt_data=pt_data,fit_face=fit_face)
    # Fit out_face
    if len(out_face)>0:
        for w in out_face[::-1]:
            fit_face=rh.fit_out_face(one_face=w,id_data=id_data,pt_data=pt_data)
            fit_face_coll.append(fit_face)
            id_data,pt_data=rh.update_datapool(id_data=id_data,pt_data=pt_data,fit_face=fit_face)
            

# Output
np.savetxt('Output/restored_top_point_data.txt',np.array(pt_data),
           header='New vertex id, X, Y, Z_meter')
