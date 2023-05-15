# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 14:45:01 2020

@author: syang34
"""

import numpy as np
import map_back_geo_lib as mbg
from joblib import Parallel, delayed


#
flat_geo=np.loadtxt('Output/Interpolation_nondipping_results.txt',skiprows=1)
flat_geo=flat_geo[np.where(flat_geo[:,3] != 300.0)[0], :]
flat_geo=flat_geo[np.where(flat_geo[:,3] != 10.0)[0], :]

top_surface=np.loadtxt('Input/top_surface_point.txt',skiprows=1)
bot_surface=np.loadtxt('Input/bot_surface_point.txt',skiprows=1)

flat_top_surface=np.loadtxt('Output/restored_top_point_data.txt',skiprows=1)
flat_top_surface=flat_top_surface[flat_top_surface[:,0].argsort()]
flat_top_surface[:,-1]=300. #Set flat top surface elevation_meter
flat_bot_surface=np.loadtxt('Output/restored_bot_point_data.txt',skiprows=1)
flat_bot_surface=flat_bot_surface[flat_bot_surface[:,0].argsort()]
flat_bot_surface[:,-1]=10. #Set flat bottom surface elevation_meter

faces=np.loadtxt('Input/triangle_faces.txt',skiprows=1).astype('int32')

# Organiz geology in flat unit column by column
columnID=np.unique(flat_geo[:,0])
def organize_geology(cid):
    ind=np.where(flat_geo[:,0]==cid)[0]
    column=flat_geo[ind]
    return column[:,1::].tolist()

flat_col_geo=Parallel(n_jobs=-1,verbose=1)(delayed(organize_geology)(i) for i in columnID)
problem=[]
for column in flat_col_geo:
    seq=[p[2] for p in column]
    if all(m > n for m, n in zip(seq, seq[1::])) != 1:
        problem.append(column)

if len(problem)>0:
    print('problem in organize geology in flat unit')

# Search possible triangles for each column
def search_triangles(col,surface_pt):
    col_xy=[col[0][0],col[0][1]]
    dist2=(surface_pt[:,1]-col_xy[0])**2 + (surface_pt[:,2]-col_xy[1])**2
    ind=np.where(dist2 == min(dist2))[0]
    npt=int(surface_pt[ind[0]][0]) # nearest point id
    p_tri=[]
    for i in [faces[:,1], faces[:,2], faces[:,3]]:
        inde=np.where(i==npt)[0]
        p_tri += faces[inde].tolist()
    #
    return np.array(p_tri)#,np.sqrt(min(dist2))

#test=search_triangles(grid2d[5000:5001,1::], flat_top_surface)
top_tri_coll=Parallel(n_jobs=-1,verbose=1)(delayed(search_triangles)(col,flat_top_surface) for col in flat_col_geo)
bot_tri_coll=Parallel(n_jobs=-1,verbose=1)(delayed(search_triangles)(col,flat_bot_surface) for col in flat_col_geo)

# Map one column in flat unit to dipped unit
def ptintri(tri2d,pt):
    x1,y1=tri2d[0][0],tri2d[0][1]
    x2,y2=tri2d[1][0],tri2d[1][1]
    x3,y3=tri2d[2][0],tri2d[2][1]
    xp,yp=pt[0],pt[1]
    c1 = (x2-x1)*(yp-y1)-(y2-y1)*(xp-x1)
    c2 = (x3-x2)*(yp-y2)-(y3-y2)*(xp-x2)
    c3 = (x1-x3)*(yp-y3)-(y1-y3)*(xp-x3)
    if (c1<0 and c2<0 and c3<0) or (c1>0 and c2>0 and c3>0):
        return 1
    else:
        return 0

#
def map_one_column(col, ptri_top, ptri_bot):
    col_xy=[col[0][0],col[0][1]]
    fall_top_face_vertex=[]
    fall_bot_face_vertex=[]
    #point = Point(col_xy[0], col_xy[1])
    for i in ptri_top: # Search top triangle where col_xy falls
        face_pt=flat_top_surface[i[1::]].tolist()
        tri=[p[1:3] for p in face_pt]
        judge=ptintri(tri2d=tri, pt=col_xy)
        if judge==1:
            fall_top_face_vertex=i
    #
    for i in ptri_bot: # Search bot triangle where col_xy falls
        face_pt=flat_bot_surface[i[1::]].tolist()
        tri=[p[1:3] for p in face_pt]
        judge=ptintri(tri2d=tri, pt=col_xy)
        if judge==1:
            fall_bot_face_vertex=i
    #
    if len(fall_top_face_vertex)>0 and len(fall_bot_face_vertex)>0:
        fall_face_flat_top=flat_top_surface[fall_top_face_vertex[1::]].tolist() # Top of flat face where col_xy falls
        fall_face_flat_bot=flat_bot_surface[fall_bot_face_vertex[1::]].tolist() # Bottom of flat face where col_xy falls
        pt_top_flat=col[0][0:2]+[300.0]#col[0][0:3] # intersect point of col and the top flat face
        pt_bot_flat=col[-1][0:2]+[10.0]#col[-1][0:3] # intersect point of col and the bot flat face
        #
        fall_face_top=top_surface[fall_top_face_vertex[1::]] # Top of dipped face where col_xy falls
        fall_face_bot=bot_surface[fall_bot_face_vertex[1::]] # Bottom of dipped face where col_xy falls
        #
        pt_top=mbg.calculate_dip_surface_pt(dipped_surface=fall_face_top,
                                            flat_surface=fall_face_flat_top,pt=pt_top_flat)[0] # location of col top on top dip face
        pt_bot=mbg.calculate_dip_surface_pt(dipped_surface=fall_face_bot,
                                            flat_surface=fall_face_flat_bot,pt=pt_bot_flat)[0] # location of col bot on bottom dip face
        #
        mapped_geo=mbg.map_col_geo(col_geo=col,dip_top=pt_top.tolist(),dip_bot=pt_bot.tolist(),
                                   flat_top=pt_top_flat,flat_bot=pt_bot_flat)
        #
        return mapped_geo
    else:
        return []

#
# Main parallel program
mapped_geo_coll_para=Parallel(n_jobs=-1,verbose=1)(delayed(map_one_column)(col, ptri_top, ptri_bot) for col, ptri_top, ptri_bot in zip(flat_col_geo,top_tri_coll,bot_tri_coll))


mapped_geo_coll_clean=[]
for i in mapped_geo_coll_para:
    if len(i) > 0:
        mapped_geo_coll_clean.append(i)

# Export
#  Add layer id of each column of geology before export
mapped_geo_coll_clean_addid=[]
for i in mapped_geo_coll_clean:
    col_plus_id=[]
    p=0
    for j in i:
        col_plus_id.append([p]+j)
        p+=1
    #
    mapped_geo_coll_clean_addid.append(col_plus_id)
#  End of adding layer id
#
np.savetxt('Output/Interpolation_results_map_back.txt',np.array([q for p in mapped_geo_coll_clean_addid for q in p]),
           header='layerID X Y Z Indicator',fmt='%d %.2f %.2f %.2f %d')
