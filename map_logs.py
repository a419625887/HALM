# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 16:45:10 2020

@author: syang34
"""

import numpy as np
import csv
import matplotlib.pyplot as plt
import map_elogs_lib as mel

# make up elog sequence


with open('Input/Logs.csv', newline='') as f:
    reader = csv.reader(f,delimiter=',')
    elogs = list(reader)
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

#
tri_coll=[]
for i in faces:
    face_pt=top_surface[i[1::]].tolist()
    tri_coll.append([p[1:3] for p in face_pt])

#
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
cutlog_len=[]
mapped_log_coll=[]
intersection = []
c=1
for count in range(0,len(elogs),2):
    if c%100 == 0:
        print(c)
    c+=1
    log=[elogs[count],elogs[count+1]]
    log_xy=[float(log[1][1]),float(log[1][2])]
    fall_face_vertex=[]
    #point = Point(log_xy[0], log_xy[1])
    for i,tri in zip(faces,tri_coll):
        judge=ptintri(tri2d=tri, pt=log_xy)
        if judge==1:
            fall_face_vertex=i
    #
    #
    if len(fall_face_vertex) > 0:
        fall_face_top=top_surface[fall_face_vertex[1::]].tolist() # Top of dipped face where elog falls
        pt_top=mel.vertex_interpo(sample_info=fall_face_top,query_info=log_xy) # intersect point of elog and the top face
        fall_face_bot=bot_surface[fall_face_vertex[1::]].tolist() # Bottom of dipped face where elog falls
        pt_bot=mel.vertex_interpo(sample_info=fall_face_bot,query_info=log_xy) # intersect point of elog and the bottom face
        #print([pt_top[2], pt_bot[2]])
        intersection.append(pt_top + pt_bot)
        #
        fall_face_flat_top=flat_top_surface[fall_face_vertex[1::]] # Top of flat face where elog falls
        fall_face_flat_bot=flat_bot_surface[fall_face_vertex[1::]] # Bottom of flat face where elog falls
        #
        intersect_top=mel.calculate_flat_surface_pt(dipped_surface=fall_face_top,
                                                    flat_surface=fall_face_flat_top,pt=pt_top)
        intersect_bot=mel.calculate_flat_surface_pt(dipped_surface=fall_face_bot,
                                                    flat_surface=fall_face_flat_bot,pt=pt_bot)
        #
        pt_top_flat=mel.select_intersect(intersect=intersect_top,flat_surface=fall_face_flat_top) # Elog location in top flat surface
        pt_bot_flat=mel.select_intersect(intersect=intersect_bot,flat_surface=fall_face_flat_bot) # Elog location in bottom flat surface
        #
        cut_log=mel.cut_elog(log=log,up_bound=pt_top[-1]/0.3048,lw_bound=pt_bot[-1]/0.3048) # Upper and lower bound are in foot
        if cut_log != ['no data']:
            cutlog_len.append(cut_log[0][-1]-cut_log[0][3])
            mapped_log=mel.map_log_seq(cut_log=cut_log,dip_top=pt_top,dip_bot=pt_bot,
                                       flat_top=pt_top_flat,flat_bot=pt_bot_flat)
            mapped_log_coll.append([[log[0][0]]+m+[n] for m,n in zip(mapped_log[0],mapped_log[1])])
    else:
        print(log_xy)

# Export
with open('Output/Mapped_logs.csv', "w",newline='') as output:
    writer = csv.writer(output)
    writer.writerows([['The maximum lengh of cut log is {} ft'.format(max(cutlog_len))]])
    writer.writerows([q for p in mapped_log_coll for q in p])

with open('Output/log_intersection.txt', "w",newline='') as output:
    writer = csv.writer(output)
    writer.writerow(['X', 'Y', 'Z_top', 'X', 'Y', 'Z_bot'])
    writer.writerows(intersection)
