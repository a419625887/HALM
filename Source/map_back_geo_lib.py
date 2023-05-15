# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 17:39:49 2020

@author: syang34
"""

# This module is used to map elogs
import math
import numpy as np
from scipy.interpolate import griddata
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

#
def vertex_interpo(sample_info,query_info):
    x=np.array(sample_info)[:,1]
    y=np.array(sample_info)[:,2]
    z=np.array(sample_info)[:,3]
    #
    xq=query_info[0]
    yq=query_info[1]
    zq=griddata((x,y),z,(xq,yq),method='linear')
    #
    result=[xq, yq, float(zq)]
    return result

def tri_area(tri_list):
    verte=np.array(tri_list)
    side1=verte[0]-verte[1]
    side2=verte[0]-verte[2]
    cos_angle=np.dot(side1,side2)/(np.sqrt(sum(side1**2))*np.sqrt(sum(side2**2)))
    area=0.5*np.sqrt(sum(side1**2))*np.sqrt(sum(side2**2))*(1-cos_angle**2)**0.5
    return area

def calculate_area(vertex_list,pt):
    tri1=[pt, vertex_list[0][1::], vertex_list[1][1::]]
    tri2=[pt, vertex_list[0][1::], vertex_list[2][1::]]
    tri3=[pt, vertex_list[1][1::], vertex_list[2][1::]]
    s1=tri_area(tri_list=tri1)
    s2=tri_area(tri_list=tri2)
    s3=tri_area(tri_list=tri3)
    return s1,s2,s3


def calculate_distance(dipped_face,flat_face,pt):
    fs1,fs2,fs3=calculate_area(vertex_list=flat_face,pt=pt)
    fall_face_dipped_area=tri_area(tri_list=dipped_face[:,1::].tolist()) # area of dipped face
    #
    ls1=np.sqrt(sum((dipped_face[0][1::]-dipped_face[1][1::])**2)) # lenth of the first side of dipped triangle
    ls2=np.sqrt(sum((dipped_face[0][1::]-dipped_face[2][1::])**2)) # second side
    ls3=np.sqrt(sum((dipped_face[1][1::]-dipped_face[2][1::])**2)) # third side
    #
    d1=2*fall_face_dipped_area*fs1/(fs1+fs2+fs3)/ls1
    d2=2*fall_face_dipped_area*fs2/(fs1+fs2+fs3)/ls2
    #d3=2*fall_face_dipped_area*fs3/(fs1+fs2+fs3)/ls3
    #
    return d1,d2

def normalize_v3(vec):
    ''' Normalize a numpy array of 2 component vectors '''
    lens = np.sqrt( vec[0]**2 + vec[1]**2 +vec[2]**2)
    vec[0] /= lens
    vec[1] /= lens  
    vec[2] /= lens             
    return np.array(vec.tolist())


def find_move_dire(target,line_start,line_end):
    t=-np.dot(line_start-target,line_end-line_start)/sum((line_end-line_start)**2)
    perpend_pt=line_start+(line_end-line_start)*t
    dire=target-perpend_pt
    #
    return  normalize_v3(dire)
    
def calculate_intersect(moved_pt1,moved_pt2,v1,v2):
    x1,y1,z1=moved_pt1
    x2,y2,z2=moved_pt2
    a1,b1,c1=v1
    a2,b2,c2=v2
    #
    u=((b1/a1)*(x1-x2)+y2-y1)/(a2*b1/a1-b2)
    lam=(y2+b2*u-y1)/b1
    check=z1+c1*lam-z2-c2*u
    #
    intersect=moved_pt1+lam*v1
    #
    return intersect,check
    

def calculate_dip_surface_pt(dipped_surface,flat_surface,pt):
    d1,d2=calculate_distance(dipped_surface,flat_surface,pt=pt)
    #
    dire1=find_move_dire(target=dipped_surface[2][1::], line_start=dipped_surface[0][1::], line_end=dipped_surface[1][1::])
    dire2=find_move_dire(target=dipped_surface[1][1::], line_start=dipped_surface[0][1::], line_end=dipped_surface[2][1::])
    #
    moved_pt1=dipped_surface[0][1::]+dire1*d1
    moved_pt2=dipped_surface[0][1::]+dire2*d2
    v1=dipped_surface[1][1::]-dipped_surface[0][1::]
    v2=dipped_surface[2][1::]-dipped_surface[0][1::]
    for ii in range(len(v1)):
        if v1[ii] == 0:
            v1[ii] = 0.001
    for ii in range(len(v2)):
        if v2[ii] == 0:
            v2[ii] = 0.001
    #
    intersect,check=calculate_intersect(moved_pt1, moved_pt2, v1, v2)
    return intersect,check


# Linearly map column geology in flat unit to dipped unit
def map_col_geo(col_geo,dip_top,dip_bot,flat_top,flat_bot):
    flat_dist=np.sqrt(sum((np.array(flat_top)-np.array(flat_bot))**2))
    v=np.array(dip_bot)-np.array(dip_top)
    #
    map_geo_seq=[]
    for i in col_geo:
        f=i[0:3]
        dist=np.sqrt(sum((np.array(flat_top)-np.array(f))**2))
        t=dist/flat_dist
        map_pt=np.array(dip_top)+v*t
        map_geo_seq.append(map_pt.tolist()+[i[3]])
    #
    return map_geo_seq

#mapped_log=map_log_seq(cut_log=cut_log,dip_top=pt_top,dip_bot=pt_bot,
#                       flat_top=pt_top_flat,flat_bot=pt_bot_flat)
    
        