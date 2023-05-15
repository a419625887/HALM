# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 17:39:49 2020

@author: syang34
"""

# This module is used to map elogs
import math
import numpy as np
from numpy.linalg import lstsq
import csv
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

def get_line_function(points):
    x_coords, y_coords = zip(*points)
    A = np.vstack([x_coords,np.ones(len(x_coords))]).T
    m, c = lstsq(A, y_coords)[0] #Line Solution is y = {m}x + {c}
    #print("Line Solution is y = {m}x + {c}".format(m=m,c=c))
    return m,c

def calculate_intersection(slope1,slope2,intercept1,intercept2):
    x=(intercept2-intercept1)/(slope1-slope2)
    y=slope1*x+intercept1
    return [x,y]

def calculate_distance(dipped_face,flat_face,pt):
    s1,s2,s3=calculate_area(vertex_list=dipped_face,pt=pt)
    fall_face_flat_area=tri_area(tri_list=flat_face[:,1::].tolist())
    #
    ls1=np.sqrt(sum((flat_face[0][1::]-flat_face[1][1::])**2)) # lenth of the first side of flat triangle
    ls2=np.sqrt(sum((flat_face[0][1::]-flat_face[2][1::])**2)) # second side
    ls3=np.sqrt(sum((flat_face[1][1::]-flat_face[2][1::])**2)) # third side
    #
    d1=2*fall_face_flat_area*s1/(s1+s2+s3)/ls1
    d2=2*fall_face_flat_area*s2/(s1+s2+s3)/ls2
    #d3=2*fall_face_flat_area*s3/(s1+s2+s3)/ls3
    #
    return d1,d2


def calculate_flat_surface_pt(dipped_surface,flat_surface,pt):
    d1,d2=calculate_distance(dipped_surface,flat_surface,pt=pt)
    #
    m1,c1=get_line_function(points=[tuple(flat_surface[0][1:3]),tuple(flat_surface[1][1:3])])
    m2,c2=get_line_function(points=[tuple(flat_surface[0][1:3]),tuple(flat_surface[2][1:3])])
    angle1=math.atan(m1)
    angle2=math.atan(m2)
    diff1=abs(d1/math.cos(angle1))
    diff2=abs(d2/math.cos(angle2))
    intersecion=[]
    intersecion.append(calculate_intersection(slope1=m1,slope2=m2,intercept1=c1+diff1,intercept2=c2+diff2))
    intersecion.append(calculate_intersection(slope1=m1,slope2=m2,intercept1=c1+diff1,intercept2=c2-diff2))
    intersecion.append(calculate_intersection(slope1=m1,slope2=m2,intercept1=c1-diff1,intercept2=c2+diff2))
    intersecion.append(calculate_intersection(slope1=m1,slope2=m2,intercept1=c1-diff1,intercept2=c2-diff2))
    #
    return intersecion

def select_intersect(intersect,flat_surface):
    pt_on_flat=[]
    for i in intersect:
        point = Point(i[0], i[1])
        polygon = Polygon([tuple(p[1:3]) for p in flat_surface.tolist()])
        judge=polygon.contains(point)
        if judge==1:
            pt_on_flat.append(i+[flat_surface[0][-1]])
    #
    if len(pt_on_flat)==1:
        return pt_on_flat[0]
    else:
        return ['Problem']

# Cut a Elog sequence given upper bound and lower bound
def cut_elog(log,up_bound,lw_bound): # input of up_bound and lw_bound should be in foot
    if '' in log[0]:
        log_clean=[log[0][0:log[0].index('')],log[1][0:log[0].index('')]]
    else:
        log_clean=log
    #
    datum=float(log_clean[1][3])
    dep_seq=np.array([float(p) for p in log_clean[0][4::]])
    soil_seq=[p for p in log_clean[1][4::]]
    #
    if datum-dep_seq[-1]<up_bound and datum>lw_bound:
        if datum-lw_bound>dep_seq[-1]:
            lw_bound=round(datum-dep_seq[-1],1)
    #
        inde=np.where(np.around(datum-dep_seq,1)>=up_bound)[0]
        indx=np.where(np.around(datum-dep_seq,1)<=lw_bound)[0]
    #
        if len(inde)==0 and up_bound<=datum:
            cut_dep_seq=[datum-up_bound]+dep_seq[0:indx[0]].tolist()+[datum-lw_bound]
            cut_soil_seq=soil_seq[0:indx[0]+1]
        elif len(inde)==0 and up_bound>datum:
            cut_dep_seq=[0]+dep_seq[0:indx[0]].tolist()+[datum-lw_bound]
            cut_soil_seq=soil_seq[0:indx[0]+1]
        else:
            cut_dep_seq=[datum-up_bound]+dep_seq[inde[-1]+1:indx[0]].tolist()+[datum-lw_bound]
            cut_soil_seq=soil_seq[inde[-1]+1:indx[0]+1]
    #
        return [log_clean[0][0:3]+cut_dep_seq,log_clean[1][0:4]+cut_soil_seq]
    else:
        return ['no data']

#cut_log=cut_elog(log=log,up_bound=-237,lw_bound=-260)

# Descritilize cut log sequence into 1ft by 1ft data
def descritilize_log(cut_log,resl):
    dep_seq=np.array(cut_log[0][4::])
    soil_seq=cut_log[1][4::]
    #
    des_dep_seq=[]
    des_soil_seq=[]
    stop=0
    dep=cut_log[0][3]
    while stop==0:
        ind=np.where(dep_seq>=dep+resl)[0]
        if len(ind)>0:
            des_dep_seq.append(dep+resl)
            des_soil_seq.append(soil_seq[ind[0]])
        else:
            stop=1
        #
        dep+=resl
    #
    if des_dep_seq[-1]<dep_seq[-1]:
        des_dep_seq.append(dep_seq[-1])
        des_soil_seq.append(soil_seq[-1])
    #
    return [cut_log[0][0:4]+des_dep_seq,cut_log[1][0:4]+des_soil_seq]


#des_log=descritilize_log(cut_log=cut_log,resl=1.5)

# Linearly map log sequence data in dipped unit to flat unit
def map_log_seq(cut_log,dip_top,dip_bot,flat_top,flat_bot):
    datum=float(cut_log[1][3])
    elev_seq=(datum-np.array(cut_log[0][3::]))*0.3048 # elevation in meter
    soil_seq=cut_log[1][4::]
    dip_dist=np.sqrt(sum((np.array(dip_top)-np.array(dip_bot))**2))
    #flat_dist=np.sqrt(sum((np.array(flat_top)-np.array(flat_bot))**2))
    v=np.array(flat_bot)-np.array(flat_top)
    #
    map_pt_seq=[]
    for i in elev_seq:
        f=dip_top[0:2]+[i]
        dist=np.sqrt(sum((np.array(dip_top)-np.array(f))**2))
        t=dist/dip_dist
        map_pt=np.array(flat_top)+v*t
        map_pt_seq.append(map_pt.tolist())
    #
    map_soil_seq=['']+soil_seq
    return [map_pt_seq,map_soil_seq]

#mapped_log=map_log_seq(cut_log=cut_log,dip_top=pt_top,dip_bot=pt_bot,
#                       flat_top=pt_top_flat,flat_bot=pt_bot_flat)
    
        