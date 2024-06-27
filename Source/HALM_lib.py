# -*- coding: utf-8 -*-
"""
Created on Sun Aug  6 22:17:52 2023

@author: syang34
"""

# Package for HALM
import pyvista as pv
import numpy as np
import math
from scipy.interpolate import griddata
from scipy.spatial import Delaunay
from metpy.interpolate import geometry
from metpy.interpolate.points import natural_neighbor_point
import vtk

# module for horizon restoration method assuming flexural slip deformation
class Flexural_slip():
    def __init__(self, vertices, faces):
        surf = pv.PolyData(np.around(vertices, 2), faces)
        surf.compute_normals(cell_normals=True, point_normals=False, inplace=True)
        self.faces = surf.faces.reshape(-1, 4)
        self.vertices=surf.points
        self.norm=surf['Normals']
        self.centers=surf.cell_centers().points
        
        # resolution of horizon surface
        self.sur_resol = np.sum((self.vertices[0, 0:2] - self.vertices[1, 0:2]) ** 2) ** 0.5
        
        #
        self.pt_data = np.zeros(self.vertices.shape) # vertex of fitted triangles
        self.pt_data = np.hstack([np.arange(len(self.pt_data)).reshape(-1, 1), self.pt_data])
        self.fit_flag = np.zeros(len(self.vertices)).astype('int64') # a flag whether a vertex has been fitted
        self.strtrow = None # starting row to fit triangles
        
        #
        self.rotated_face = self.rotate_triangles()
        self.row_faceid_coll, self.n_row_faces, self.i_row_faces_gap = self.get_row_face()
        
    # flatten triangle faces by rotation
    def rotate_triangles(self):
        rotated_face = np.zeros((self.faces.shape[0] * 3, 4))
        for i in range(self.faces.shape[0]):
            n = self.norm[i].copy()
            n[2] = np.abs(n[2])
            centroid = self.centers[i].copy()
            face_pt = self.vertices[self.faces[i][1::], :].copy()
            #
            x = (n[2]**2 / (n[1]**2 + n[2]**2)) ** 0.5
            d = np.zeros(3)
            d[0] = x*n[0]
            d[1] = x*n[1]
            d[2] = -(n[0]**2 + n[1]**2) ** 0.5
            s = np.cross(n,d)
            #
            if face_pt[0, 2] == face_pt[1, 2] == face_pt[2, 2]:
                # No need to rotate if there is no dip for a face
                rotated_pt = np.array([p[:] for p in face_pt])
            else:
                rotated_pt = np.zeros((3, 3))
                for j in range(face_pt.shape[0]):
                    p = face_pt[j]
                    t = -np.dot(centroid - p, s) / sum(s**2)
                    v = centroid + s * t
                    dist = sum((p-v)**2)**0.5
                    dire = normalize_v2((p-v)[0:2])
                    #
                    if np.isnan(dire[0]) == 1 and round(dist, 3) == 0:
                        # Handle the speical case that the strike passes through a node of a triangle
                        node = p.copy()
                    else:
                        node = v + dire*dist
                    #
                    rotated_pt[j] = node.copy()
            #
            rotated_pt = np.hstack([self.faces[i][1::].reshape(-1, 1), rotated_pt])
            rotated_face[3*i:3*i+3, :] = rotated_pt
        #
        return rotated_face # vertex ID, x, y, z
    #
    # search triangles for each row with the ascending order of x
    def get_row_face(self):
        xs = np.unique(self.vertices[:, 0])
        ys = np.unique(self.vertices[:, 1])
        row_faceid_coll = np.ones((len(ys) - 1, len(xs) * 2)) * (-1)
        row_faceid_coll = row_faceid_coll.astype('int64') # face id for each row
        n_row_faces = np.zeros(len(ys) - 1).astype('int64') # number of faces for each row
        i_row_faces_gap = np.zeros(len(ys) - 1).astype('int64') # where there are gaps of face
                                                                # in a row, like a hole on surface.
        #
        for i in range(len(ys) - 1):
            y1 = ys[i]
            y2 = ys[i+1]
            ind = np.where((self.centers[:, 1] > y1) & (self.centers[:, 1] < y2))[0]
            center_x = self.centers[ind, 0]
            ind_sort = ind[np.argsort(center_x)] # face id with increasing x order in a row
            #
            row_faceid_coll[i][0:len(ind_sort)] = ind_sort.copy()
            n_row_faces[i] = len(ind_sort)
            #
            center_x_sort = center_x[np.argsort(center_x)]
            x_dist = center_x_sort[1::] - center_x_sort[0:-1]
            assert all(m > 0 for m in x_dist) == 1
            if np.sum(x_dist > self.sur_resol) > 0:
                i_row_faces_gap[i] = 1
        #
        return row_faceid_coll, n_row_faces, i_row_faces_gap
    #
    # Find starting row
    def find_start_row(self):
        ind = np.where(self.n_row_faces == np.max(self.n_row_faces))[0]
        strtrow = ind[0]
        #
        return strtrow
    #
    #
    def find_start_row_v2(self):
        ind_sort = np.argsort(self.n_row_faces)[::-1] # decreasing order
        for i in ind_sort:
            strtrow = i
            i_gap = self.i_row_faces_gap[i]
            if i_gap == 0:
                break
        #
        return strtrow
    #
    # fit triangles in the starting row
    def fit_first_row(self, startrow = None):
        if startrow is None:
            self.strtrow = self.find_start_row_v2()
        else:
            self.strtrow = startrow
        #
        onerow = self.row_faceid_coll[self.strtrow].copy()
        onerow = onerow[onerow >= 0]
        # fix the first face
        i = onerow[0]
        tri_pt = self.rotated_face[3*i:3*i+3]
        tri_pt_id = tri_pt[:, 0].astype('int64')
        self.pt_data[tri_pt_id] = tri_pt
        self.fit_flag[tri_pt_id] = 1
        #
        for i in onerow:
            tri_pt = self.rotated_face[3*i:3*i+3]
            tri_pt_id = tri_pt[:, 0].astype('int64')
            i_fit = self.fit_flag[tri_pt_id]
            #
            share_pt = tri_pt[i_fit == 1]
            unshared_pt = tri_pt[i_fit == 0]
            #
            if len(share_pt) == 3:
                fit_face = share_pt
            elif len(share_pt) == 2:
                fit_face, _ = self.rigid_fit(shared_pt=share_pt,unshared_pt=unshared_pt)
            #
            _ = self.update_datapool(fit_face = fit_face)
        #
        return 1
    #
    # fit remaining rows
    def fit_remaining_rows(self):
        pending_faceid_coll1 = []
        for row in self.row_faceid_coll[self.strtrow+1::]:
            _, pending_faceid = self.fit_other_row(row = row)
            if len(pending_faceid) > 0:
                pending_faceid_coll1.append(pending_faceid)
        #
        pending_faceid_coll2 = []
        for row in self.row_faceid_coll[0:self.strtrow][::-1]:
            _, pending_faceid = self.fit_other_row(row = row[::-1])
            if len(pending_faceid) > 0:
                pending_faceid_coll2.append(pending_faceid)
        #
        pending_faceid_coll = pending_faceid_coll1[::-1] + pending_faceid_coll2[::-1]
        #
        t = 1
        while len(pending_faceid_coll) > 0:
            temp = []
            for row in pending_faceid_coll:
                _, pending_faceid = self.fit_pending_row(row = row)
                if len(pending_faceid) > 0:
                    temp.append(pending_faceid)
            #
            pending_faceid_coll = [p[:] for p in temp]
            if t % 2 == 1:
                pending_faceid_coll = pending_faceid_coll[::-1]
            #
            print('Iteration {} to fit pending faces.'.format(t))
            t += 1
        #
        return 1
    #
    # fit a row other than the starting row
    def fit_other_row(self, row):
        onerow = row.copy()
        onerow = onerow[onerow >= 0]
        #
        indi=0 # whether the previous fit is a three-point fit or an awkward fit
        #istrainfit=0 # whether the previous fit is a strain fit
        out_faceid = np.zeros(len(onerow)).astype('int64')
        out_faceid[:] = -1
        cc = 0
        #
        for i in onerow:
            tri_pt = self.rotated_face[3*i:3*i+3]
            tri_pt_id = tri_pt[:, 0].astype('int64')
            i_fit = self.fit_flag[tri_pt_id]
            #
            share_pt = tri_pt[i_fit == 1]
            unshared_pt = tri_pt[i_fit == 0]
            #
            i_outface = 0 # whether the current triangle is a outface
            if len(share_pt)==2 and indi==0:
                fit_face, indi = self.rigid_fit(shared_pt=share_pt,unshared_pt=unshared_pt)
                face_bef_fit = np.vstack([share_pt, unshared_pt])
                previous_fit_tri = np.array([p[:] for p in fit_face])
                #istrainfit=0
            elif len(share_pt)==2 and indi==1:
                fit_face = self.strain_fit(shared_pt=share_pt,unshared_pt=unshared_pt,
                                           previous_fit_tri=previous_fit_tri,previous_unfit_tri=face_bef_fit)
                indi=0
                #istrainfit=1
            #elif len(share_pt)==2 and indi==0 and istrainfit==1:
            #    fit_face=rh.awk_fit(datapool=pt_data,shared_pt=share_pt,unshared_pt=unshared_pt)
            #    face_bef_fit=share_pt+unshared_pt
            #    indi=1
            #    istrainfit=0
            #    fit_face_coll.append(fit_face)
            elif len(share_pt)==3:
                fit_face = self.pt_data[tri_pt_id]
                face_bef_fit= np.array([p[:] for p in share_pt])
                previous_fit_tri = np.array([p[:] for p in fit_face])
                indi=1
                #istrainfit=0
            else:
                out_faceid[cc] = i
                i_outface = 1
            #
            cc += 1
            # update datapool
            if i_outface == 0:
                _ = self.update_datapool(fit_face = fit_face)
        # Fit out_face
        out_faceid = out_faceid[out_faceid >= 0]
        pending_faceid = np.zeros(len(out_faceid)).astype('int64')
        pending_faceid[:] = -1
        tt = 0
        if len(out_faceid) > 0:
            for w in out_faceid[::-1]:
                fit_face, pending = self.fit_out_face(one_faceid = w)
                if pending == 0:
                    _ = self.update_datapool(fit_face = fit_face)
                else:
                    pending_faceid[tt] = w
                    tt += 1
        #
        pending_faceid = pending_faceid[pending_faceid >= 0]
        #
        return 1, pending_faceid
    #
    # Rigid translation to fit triangles with two shared nodes
    def rigid_fit(self, shared_pt, unshared_pt):
        p1 = shared_pt[0].copy()
        p2 = shared_pt[1].copy()
        p3 = unshared_pt[0].copy()
        dp1 = self.pt_data[int(p1[0])].copy()
        dp2 = self.pt_data[int(p2[0])].copy()
        #
        trans=dp1[1::]-p1[1::]
        p1_trans=p1[1::]+trans
        p2_trans=p2[1::]+trans
        #
        v1=dp2[1::]-dp1[1::]
        v2=p2_trans-p1_trans
        cos=np.dot(v1,v2)/(np.sqrt(sum(v1**2))*np.sqrt(sum(v2**2)))
        cos = round(cos, 3) # precision issue
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
            fit_results = np.zeros((3, 4))
            fit_results[0, 0] = p1[0]
            fit_results[0, 1::] = p1_trans
            fit_results[1, 0] = p2[0]
            fit_results[1, 1::] = np.array(p2_fit)
            fit_results[2, 0] = p3[0]
            fit_results[2, 1::] = np.array(p3_fit)
            #fit_results=[[p1[0]]+p1_trans.tolist(),[p2[0]]+p2_fit,[p3[0]]+p3_fit]
            distort=0
        else:
            #trans2=dp2[1::]-np.array(p2_fit)
            #p2_fit=np.array(p2_fit)+trans2
            #p3_fit=np.array(p3_fit)+trans2
            #fit_results=[[p1[0]]+p1_trans.tolist(),[p2[0]]+p2_fit.tolist(),[p3[0]]+p3_fit.tolist()]
            fit_results=self.awk_fit(shared_pt,unshared_pt)
            distort=1
        #
        return fit_results, distort
    #
    # Fit face without three shared points following a strain-fitted face
    def awk_fit(self, shared_pt, unshared_pt):
        p1 = shared_pt[0].copy()
        p2 = shared_pt[1].copy()
        p3 = unshared_pt[0].copy()
        dp1 = self.pt_data[int(p1[0])].copy()
        dp2 = self.pt_data[int(p2[0])].copy()
        #
        trans1=dp1[1::]-p1[1::]
        p1_trans=p1[1::]+trans1
        p2_med=p2[1::]+trans1
        trans2=dp2[1::]-p2_med
        p2_trans=p2_med+trans2
        p3_trans=p3[1::]+trans1+trans2
        #
        fit_results = np.zeros((3, 4))
        fit_results[0, 0] = p1[0]
        fit_results[0, 1::] = p1_trans
        fit_results[1, 0] = p2[0]
        fit_results[1, 1::] = p2_trans
        fit_results[2, 0] = p3[0]
        fit_results[2, 1::] = p3_trans
        #fit_results=[[p1[0]]+p1_trans.tolist(),[p2[0]]+p2_trans.tolist(),[p3[0]]+p3_trans.tolist()]
        #
        return fit_results
    #
    # Strain fit for the triangle following a triangle fixed by 3 nodes
    def strain_fit(self, shared_pt, unshared_pt, previous_fit_tri, previous_unfit_tri):
        p1 = shared_pt[0].copy()
        p2 = shared_pt[1].copy()
        p3 = unshared_pt[0].copy()
        dp1 = self.pt_data[int(p1[0])].copy()
        dp2 = self.pt_data[int(p2[0])].copy()
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
        tri2_unfit_area=tri_area(tri_list=np.vstack([shared_pt, unshared_pt])) # Area of the focused triangle before fitting
        #
        t=-np.dot(p1_trans-p3_trans,p2_trans-p1_trans)/sum((p2_trans-p1_trans)**2)
        v=p1_trans+(p2_trans-p1_trans)*t # Projection of p3_trans on the line passing p1_trans_and p2_trans
        d=2*(tri1_unfit_area+tri2_unfit_area-tri1_fit_area)/np.sqrt(sum((p2_trans-p1_trans)**2)) # Distance v and p3_fit
        D=np.sqrt(sum((p3_trans-v)**2)) # Distance between v and p3_trans
        p3_fit=v+(p3_trans-v)*(d/D)
        #
        fit_results = np.zeros((3, 4))
        fit_results[0, 0] = p1[0]
        fit_results[0, 1::] = p1_trans
        fit_results[1, 0] = p2[0]
        fit_results[1, 1::] = p2_trans
        fit_results[2, 0] = p3[0]
        fit_results[2, 1::] = p3_fit
        #fit_results=[[p1[0]]+p1_trans.tolist(),[p2[0]]+p2_trans.tolist(),[p3[0]]+p3_fit.tolist()]
        #
        return fit_results
    #
    # Fit out_face using rigid fit method
    def fit_out_face(self, one_faceid):
        tri_pt = self.rotated_face[3*one_faceid:3*one_faceid+3]
        tri_pt_id = tri_pt[:, 0].astype('int64')
        i_fit = self.fit_flag[tri_pt_id]
        #
        share_pt = tri_pt[i_fit == 1]
        unshared_pt = tri_pt[i_fit == 0]
        #
        if len(share_pt) == 2:
            fit_face, _ = self.rigid_fit(shared_pt=share_pt,unshared_pt=unshared_pt)
            pending = 0
        elif len(share_pt) == 3: # to handle special grid configuration, should rarely occur.
            fit_face = self.pt_data[tri_pt_id]
            pending = 0
        else:
            fit_face = None
            pending = 1
        #
        return fit_face, pending
    #
    # Fit pending faces
    def fit_pending_row(self, row):
        onerow = row.copy()
        assert all(p >= 0 for p in onerow) == 1
        #
        indi=0 # whether the previous fit is a three-point fit or an awkward fit
        out_faceid = np.zeros(len(onerow)).astype('int64')
        out_faceid[:] = -1
        cc = 0
        #
        for i in onerow:
            tri_pt = self.rotated_face[3*i:3*i+3]
            tri_pt_id = tri_pt[:, 0].astype('int64')
            i_fit = self.fit_flag[tri_pt_id]
            #
            share_pt = tri_pt[i_fit == 1]
            unshared_pt = tri_pt[i_fit == 0]
            #
            i_outface = 0 # whether the current triangle is a outface
            if len(share_pt)==2 and indi==0:
                fit_face, indi = self.rigid_fit(shared_pt=share_pt,unshared_pt=unshared_pt)
                face_bef_fit = np.vstack([share_pt, unshared_pt])
                previous_fit_tri = np.array([p[:] for p in fit_face])
            elif len(share_pt)==2 and indi==1:
                fit_face = self.strain_fit(shared_pt=share_pt,unshared_pt=unshared_pt,
                                           previous_fit_tri=previous_fit_tri,previous_unfit_tri=face_bef_fit)
                indi=0
            elif len(share_pt)==3:
                fit_face = self.pt_data[tri_pt_id]
                face_bef_fit= np.array([p[:] for p in share_pt])
                previous_fit_tri = np.array([p[:] for p in fit_face])
                indi=1
            else:
                out_faceid[cc] = i
                i_outface = 1
            #
            cc += 1
            # update datapool
            if i_outface == 0:
                _ = self.update_datapool(fit_face = fit_face)
        # Fit out_face
        out_faceid = out_faceid[out_faceid >= 0]
        pending_faceid = np.zeros(len(out_faceid)).astype('int64')
        pending_faceid[:] = -1
        tt = 0
        if len(out_faceid) > 0:
            for w in out_faceid[::-1]:
                fit_face, pending = self.fit_out_face(one_faceid = w)
                if pending == 0:
                    _ = self.update_datapool(fit_face = fit_face)
                else:
                    pending_faceid[tt] = w
                    tt += 1
                    #print('Problem occurs while fitting pending faces.')
        #
        pending_faceid = pending_faceid[pending_faceid >= 0]
        #
        return 1, pending_faceid
    #
    # update pt_data after fiting a face
    def update_datapool(self, fit_face):
        fit_face_pt_id = fit_face[:, 0].astype('int64')
        self.pt_data[fit_face_pt_id] = fit_face
        self.fit_flag[fit_face_pt_id] = 1
        #
        return 1
    #
        

#
def normalize_v2(vec):
    ''' Normalize a numpy array of 2 component vectors '''
    lens = np.sqrt( vec[0]**2 + vec[1]**2 )
    vec[0] /= lens
    vec[1] /= lens               
    return np.array(vec.tolist()+[0])
#
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
#
# Compute triangle area given three vertexs
def tri_area(tri_list):
    verte=np.array(tri_list)[:,1::]
    side1=verte[0]-verte[1]
    side2=verte[0]-verte[2]
    cos_angle=np.dot(side1,side2)/(np.sqrt(sum(side1**2))*np.sqrt(sum(side2**2)))
    area=0.5*np.sqrt(sum(side1**2))*np.sqrt(sum(side2**2))*(1-cos_angle**2)**0.5
    #
    return area


# Module for mapping log data after horizon restoration
class Map_logs():
    def __init__(self, faces, top_surface_path, bot_surface_path,
                  flat_top_surface_path, flat_bot_surface_path):
        self.faces = np.array([p[:] for p in faces]).astype('int64')
        
        self.top_surface = np.loadtxt(top_surface_path)
        self.top_surface = np.hstack((np.arange(len(self.top_surface)).reshape(-1, 1), self.top_surface))
        self.bot_surface = np.loadtxt(bot_surface_path)
        self.bot_surface = np.hstack((np.arange(len(self.bot_surface)).reshape(-1, 1), self.bot_surface))

        self.flat_top_surface = np.loadtxt(flat_top_surface_path)
        self.flat_top_surface = self.flat_top_surface[self.flat_top_surface[:,0].argsort()]
        self.flat_top_surface[:,-1] = 300. #Set flat top surface elevation_meter
        self.flat_bot_surface = np.loadtxt(flat_bot_surface_path)
        self.flat_bot_surface = self.flat_bot_surface[self.flat_bot_surface[:,0].argsort()]
        self.flat_bot_surface[:,-1] = 10. #Set flat bottom surface elevation_meter
        
        #
    #
    # Map logs into non-dipping domain 
    def map_logs(self, logdata):
        cutlog_len = np.zeros(int(len(logdata) / 2))
        mapped_log_coll = np.empty((0, 5))
        count = 0
        #
        for i in range(0, len(logdata), 2):
            if int(i/2) % 1000 == 0:
                print(int(i/2), 'out of', int(len(logdata)/2), 'logs.')
            #
            onelog = [logdata[i], logdata[i+1]]
            cut_len, i_map, mapped_log = self.map_a_log(onelog = onelog)
            #
            if i_map == 1:
                cutlog_len[int(i/2)] = cut_len
                count += 1
                mapped_log_coll = np.vstack([mapped_log_coll, mapped_log])
        #
        print('Maximum length of cut log is', round(np.max(cutlog_len), 2), 'm.')
        print(count, 'out of ', int(len(logdata) / 2), 'logs are mapped.')
        return mapped_log_coll
    #
    # Map one log
    def map_a_log(self, onelog):
        log_xy, log = self.format_logdata(onelog = onelog)
        i_find, fall_face_vertex = self.locate_face_v2(log_xy = log_xy)
        #
        if i_find == 1:
            fall_face_top = self.top_surface[fall_face_vertex, 1::] # Top of dipped face where a log falls
            pt_top = self.locate_dip_surface_pt(tri_pt = fall_face_top, log_xy = log_xy) # intersect point of log and the top face
            fall_face_bot = self.bot_surface[fall_face_vertex, 1::] # Bottom of dipped face where a log falls
            pt_bot = self.locate_dip_surface_pt(tri_pt = fall_face_bot, log_xy = log_xy) # intersect point of log and the bottom face
            #
            fall_face_flat_top = self.flat_top_surface[fall_face_vertex, 1::] # Top of flat face where a log falls
            fall_face_flat_bot = self.flat_bot_surface[fall_face_vertex, 1::] # Bottom of flat face where a log falls
            #
            intersect_top=self.find_intersections_v2(dip_face=fall_face_top,
                                                  flat_face=fall_face_flat_top, pt=pt_top)
            intersect_bot=self.find_intersections_v2(dip_face=fall_face_bot,
                                                  flat_face=fall_face_flat_bot, pt=pt_bot)
            #
            status, pt_top_flat = self.select_intersect(intersect=intersect_top,
                                                        flat_face=fall_face_flat_top) # log location in top flat surface
            status, pt_bot_flat = self.select_intersect(intersect=intersect_bot,
                                                        flat_face=fall_face_flat_bot) # log location in bottom flat surface
            #
            i_cut, cut_log = self.cut_a_log(log=log, up_bound=pt_top[-1], lw_bound=pt_bot[-1])
            #
            if i_cut == 1:
                cut_len = float(cut_log[0, 3]) - float(cut_log[-1, 4])
                mapped_log = self.linear_log_map(cut_log=cut_log,dip_top=pt_top,dip_bot=pt_bot,
                                                 flat_top=pt_top_flat,flat_bot=pt_bot_flat)
                i_map = 1 # whether a log is mapped
            else:
                cut_len = 0.0
                mapped_log = None
                i_map = 0
        else:
            cut_len = 0.0
            mapped_log = None
            i_map = 0
        #
        #print(cut_log, intersect_bot, fall_face_flat_bot, pt_bot_flat, status)
        return cut_len, i_map, mapped_log
    #
    # organize a log seqence
    def format_logdata(self, onelog):
        logname = onelog[0][0]
        logx = onelog[1][1]
        logy = onelog[1][2]
        log_xy = np.array([logx, logy]).astype('float64')
        log_xy = np.around(log_xy, 1)
        datum = float(onelog[1][3])
        #
        dep_seq = onelog[0][3::]
        dep_seq = np.array(dep_seq).astype('float64')
        elev_seq = np.around((datum - dep_seq) * 0.3048, 1) # to meter
        data_seq = onelog[1][4::]
        data_seq = np.array(data_seq)
        #
        seq_len = len(data_seq)
        format_log = np.empty((seq_len, 6), dtype = 'U64')
        #
        format_log[:, 0] = logname
        format_log[:, 1] = logx
        format_log[:, 2] = logy
        format_log[:, 3] = elev_seq[0:-1]
        format_log[:, 4] = elev_seq[1::]
        format_log[:, 5] = data_seq
        #
        return log_xy, format_log
    #
    # find a face where a point falls in
    def locate_face(self, log_xy):
        i_find = 0 # whether a face is located that enclose the point
        fall_face_vertex = np.ones(3).astype('int64') # vertex ID of the located face
        fall_face_vertex[:] = -1
        for i in range(len(self.faces)):
            tri_pt = self.top_surface[self.faces[i][1::]]
            tri_pt = tri_pt[:, 1::]
            i_find = ptintri(tri2d = tri_pt, pt = log_xy)
            if i_find == 1:
                fall_face_vertex = self.faces[i][1::].copy()
                break
        #
        return i_find, fall_face_vertex
    #
    # find a face where a point falls in, more efficient version
    def locate_face_v2(self, log_xy):
        i_find = 0 # whether a face is located that enclose the point
        fall_face_vertex = np.ones(3).astype('int64') # vertex ID of the located face
        fall_face_vertex[:] = -1
        #
        x = log_xy[0]
        y = log_xy[1]
        dist2 = (self.top_surface[:, 1] - x)**2 + (self.top_surface[:, 2] - y)**2
        ind = np.where(dist2 == np.min(dist2))[0]
        near_ptid = ind[0] # nearest node id
        #
        inde = np.where(self.faces[:, 1::] == near_ptid)[0]
        near_faces = self.faces[inde]
        #
        for i in range(len(near_faces)):
            tri_pt = self.top_surface[near_faces[i][1::]]
            tri_pt = tri_pt[:, 1::]
            i_find = ptintri(tri2d = tri_pt, pt = log_xy)
            if i_find == 1:
                fall_face_vertex = near_faces[i][1::].copy()
                break
        #
        return i_find, fall_face_vertex
    #
    # locate intersection point on dipping face
    def locate_dip_surface_pt(self, tri_pt, log_xy):
        x = tri_pt[:, 0]
        y = tri_pt[:, 1]
        z = tri_pt[:, 2]
        #
        xq = log_xy[0]
        yq = log_xy[1]
        zq = griddata((x,y),z,(xq,yq), method='linear')
        #
        result = np.array([xq, yq, zq])
        result = np.around(result, 1)
        #
        return result
    #
    # find candidate intersection points on flat face
    def find_intersections(self, dip_face, flat_face, pt):
        d1,d2 = self.calculate_distance(dip_face, flat_face, pt)
        #
        m1,c1 = self.get_line_function(points = np.vstack([flat_face[0][0:2], flat_face[1][0:2]]))
        m2,c2 = self.get_line_function(points = np.vstack([flat_face[0][0:2], flat_face[2][0:2]]))
        angle1=math.atan(m1)
        angle2=math.atan(m2)
        diff1=abs(d1/math.cos(angle1))
        diff2=abs(d2/math.cos(angle2))
        intersecion = np.zeros((4, 2))
        intersecion[:] = np.nan
        intersecion[0] = self.calculate_intersection(slope1=m1,slope2=m2,intercept1=c1+diff1,intercept2=c2+diff2)
        intersecion[1] = self.calculate_intersection(slope1=m1,slope2=m2,intercept1=c1+diff1,intercept2=c2-diff2)
        intersecion[2] = self.calculate_intersection(slope1=m1,slope2=m2,intercept1=c1-diff1,intercept2=c2+diff2)
        intersecion[3] = self.calculate_intersection(slope1=m1,slope2=m2,intercept1=c1-diff1,intercept2=c2-diff2)
        #
        return intersecion
    #
    # more stable verison
    def find_intersections_v2(self, dip_face, flat_face, pt):
        d1,d2 = self.calculate_distance(dip_face, flat_face, pt)
        #
        m1,c1 = self.get_line_function_v2(points = np.vstack([flat_face[0][0:2], flat_face[1][0:2]]))
        m2,c2 = self.get_line_function_v2(points = np.vstack([flat_face[0][0:2], flat_face[2][0:2]]))
        #
        if m1 == 'inf':
            diff1 = d1
        else:
            angle1=math.atan(m1)
            diff1=abs(d1/math.cos(angle1))
        if m2 == 'inf':
            diff2 = d2
        else:            
            angle2=math.atan(m2)
            diff2=abs(d2/math.cos(angle2))
        #
        intersecion = np.zeros((4, 2))
        intersecion[:] = np.nan
        intersecion[0] = self.calculate_intersection_v2(slope1=m1,slope2=m2,intercept1=c1+diff1,intercept2=c2+diff2)
        intersecion[1] = self.calculate_intersection_v2(slope1=m1,slope2=m2,intercept1=c1+diff1,intercept2=c2-diff2)
        intersecion[2] = self.calculate_intersection_v2(slope1=m1,slope2=m2,intercept1=c1-diff1,intercept2=c2+diff2)
        intersecion[3] = self.calculate_intersection_v2(slope1=m1,slope2=m2,intercept1=c1-diff1,intercept2=c2-diff2)
        #
        return intersecion
    #
    #
    def calculate_distance(self, dip_face, flat_face, pt):
        s1,s2,s3 = self.calculate_area(vertex_list=dip_face, pt=pt)
        fall_face_flat_area = self.tri3d_area(tri_list = flat_face)
        #
        ls1=np.sqrt(sum((flat_face[0] - flat_face[1])**2)) # lenth of the first side of flat triangle
        ls2=np.sqrt(sum((flat_face[0] - flat_face[2])**2)) # second side
        ls3=np.sqrt(sum((flat_face[1] - flat_face[2])**2)) # third side
        #
        d1=2*fall_face_flat_area*s1/(s1+s2+s3)/ls1
        d2=2*fall_face_flat_area*s2/(s1+s2+s3)/ls2
        #d3=2*fall_face_flat_area*s3/(s1+s2+s3)/ls3
        #
        return d1,d2
    #
    #
    def calculate_area(self, vertex_list, pt):
        tri1=[pt, vertex_list[0], vertex_list[1]]
        tri2=[pt, vertex_list[0], vertex_list[2]]
        tri3=[pt, vertex_list[1], vertex_list[2]]
        s1=self.tri3d_area(tri_list=tri1)
        s2=self.tri3d_area(tri_list=tri2)
        s3=self.tri3d_area(tri_list=tri3)
        #
        return s1,s2,s3
    #
    #
    def tri3d_area(self, tri_list):
        verte=np.array(tri_list)
        side1=verte[0]-verte[1]
        side2=verte[0]-verte[2]
        cos_angle=np.dot(side1,side2)/(np.sqrt(sum(side1**2))*np.sqrt(sum(side2**2)))
        area=0.5*np.sqrt(sum(side1**2))*np.sqrt(sum(side2**2))*(1-cos_angle**2)**0.5
        #
        return area
    #
    #
    def get_line_function(self, points):
        pt_arr = np.array(points)
        x_coord = pt_arr[:, 0]
        y_coord = pt_arr[:, 1]
        if np.abs(x_coord[1] - x_coord[0]) < 1e-20:
            m = (y_coord[1] - y_coord[0]) / (1e-20) # slope of line
        else:
            m = (y_coord[1] - y_coord[0]) / (x_coord[1] - x_coord[0]) # slope of line
        c = y_coord[0] - m * x_coord[0] # intercept with y axis
        # Line function is y = {m}x + {c}
        #
        return m, c
    #
    # more stable version
    def get_line_function_v2(self, points):
        pt_arr = np.array(points)
        x_coord = pt_arr[:, 0]
        y_coord = pt_arr[:, 1]
        if np.abs(x_coord[1] - x_coord[0]) < 1e-10:
            m = 'inf' # line is vertical
            c = (x_coord[0] + x_coord[1]) / 2
            # Line function is x = {a}, a is denoted as c
        else:
            m = (y_coord[1] - y_coord[0]) / (x_coord[1] - x_coord[0]) # slope of line
            c = y_coord[0] - m * x_coord[0] # intercept with y axis
            # Line function is y = {m}x + {c}
        #
        return m, c
    #
    #
    def calculate_intersection(self, slope1,slope2,intercept1,intercept2):
        x=(intercept2-intercept1)/(slope1-slope2)
        y=slope1*x+intercept1
        return np.array([x,y])
    #
    # more stable version
    def calculate_intersection_v2(self, slope1,slope2,intercept1,intercept2):
        if slope1 != 'inf' and slope2 != 'inf':
            x=(intercept2-intercept1)/(slope1-slope2)
            y=slope1*x+intercept1
        elif slope1 == 'inf':
            x = intercept1
            y = slope2 * x + intercept2
        elif slope2 == 'inf':
            x = intercept2
            y = slope1 * x + intercept1
        #
        return np.array([x,y])
    #
    # select a intersection point with flat face among 4 candidate points
    def select_intersect(self, intersect, flat_face):
        pt_on_flat = np.zeros(3)
        pt_on_flat[:] = np.nan
        count = 0
        for i in range(len(intersect)):
            judge = ptintri(tri2d = flat_face[:, 0:2], pt = intersect[i])
            if judge == 1:
                pt_on_flat[0:2] = intersect[i].copy()
                pt_on_flat[2] = flat_face[0, 2]
                pt_on_flat = np.around(pt_on_flat, 1)
                count += 1
        #
        if count == 1:
            status = 'Normal'
        else:
            status = 'Problem'
        #
        return status, pt_on_flat
    #
    # Cut a log sequence given upper bound and lower bound
    def cut_a_log(self, log, up_bound, lw_bound):
        top_arr = log[:, 3].astype('float64')
        bot_arr = log[:, 4].astype('float64')
        #
        i_cut = 1 # whether a log is cut by the bound values
        if top_arr[0] > lw_bound and bot_arr[-1] < up_bound:
            ind1 = np.where(bot_arr < up_bound)[0]
            ind2 = np.where(top_arr > lw_bound)[0]
            #
            cut_log = [p[:] for p in log[ind1[0]:ind2[-1]+1]]
            cut_log = np.array(cut_log)
            #
            if up_bound < top_arr[0]:
                cut_log[0, 3] = up_bound
            if lw_bound > bot_arr[-1]:
                cut_log[-1, 4] = lw_bound
        else:
            i_cut = 0
            cut_log = np.zeros((1, 6))
            cut_log[:] = np.nan
        #
        return i_cut, cut_log
    #
    # Linearly map log sequence data from dipping domain to flat domain
    def linear_log_map(self, cut_log, dip_top, dip_bot, flat_top, flat_bot):
        elev_arr = np.zeros(cut_log.shape[0] + 1)
        elev_arr[0:-1] = cut_log[:, 3].astype('float64')
        elev_arr[-1] = float(cut_log[-1, 4])
        data_arr = cut_log[:, 5]
        #
        dip_dist = dip_top[2] - dip_bot[2]
        v = flat_bot - flat_top
        t = (dip_top[2] - elev_arr) / dip_dist
        mapped_coords = flat_top + t.reshape(-1, 1) @ v.reshape(1, -1)
        #
        mapped_log = np.empty((cut_log.shape[0]+1, 5), dtype = 'U64')
        mapped_log[:, 0] = cut_log[0, 0]
        mapped_log[:, 1:3] = np.around(mapped_coords[:, 0:2], 2)
        mapped_log[:, 3] = np.around(mapped_coords[:, 2], 1)
        mapped_log[0, 4] = ''
        mapped_log[1::, 4] = data_arr
        #
        return mapped_log
    #

    

#
# judge whether a 2D triangle enclose a point
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
#
#
# A module for lithologic modeling in non-dipping domain
class Litho_modeling():
    def __init__(self, log_nondip, faces, 
                 flat_top_surface_path, flat_bot_surface_path):
        
        #
        self.flat_top_surface = np.loadtxt(flat_top_surface_path, skiprows=1)
        self.flat_top_surface = self.flat_top_surface[self.flat_top_surface[:,0].argsort()]
        self.flat_top_surface[:,-1] = 300. #Set flat top surface elevation_meter
        self.flat_bot_surface = np.loadtxt(flat_bot_surface_path, skiprows=1)
        self.flat_bot_surface = self.flat_bot_surface[self.flat_bot_surface[:,0].argsort()]
        self.flat_bot_surface[:,-1] = 10. #Set flat bottom surface elevation_meter
        
        #
        pt_id = np.unique(faces[:, 1::])
        self.flat_top_surface = self.flat_top_surface[pt_id]
        self.flat_bot_surface = self.flat_bot_surface[pt_id]
        
        #
        self.log_nondip = np.array(log_nondip)
        self.lognames = np.unique(self.log_nondip[:, 0])
        
        #
        self.grid2d_nondip = None # 2D grid in non-dipping domain
        self.layer_elev = None # layer elevations for quasi-3D modeling
        
        #
        
    #
    # Construct 2D grid in non-dipping domain
    def generate_2dgrid(self, model_res, n_layer, h_res = None):
        '''
        model_res: Horizontal resolution of the final model
        h_res: Horizontal resolution of the non-dipping model
        n_layer: number of layer for quasi-3D lithologic modeling in non-dipping domain
        
        '''
        
        top_xmax,top_xmin=max(self.flat_top_surface[:,1]),min(self.flat_top_surface[:,1])
        top_ymax,top_ymin=max(self.flat_top_surface[:,2]),min(self.flat_top_surface[:,2])
        bot_xmax,bot_xmin=max(self.flat_bot_surface[:,1]),min(self.flat_bot_surface[:,1])
        bot_ymax,bot_ymin=max(self.flat_bot_surface[:,2]),min(self.flat_bot_surface[:,2])
        #
        xmax=max([top_xmax,bot_xmax])
        xmin=min([top_xmin,bot_xmin])
        ymax=max([top_ymax,bot_ymax])
        ymin=min([top_ymin,bot_ymin])
        zmax = self.flat_top_surface[0, -1]
        zmin = self.flat_bot_surface[0, -1]
        #
        if h_res is None:
            h_res = model_res * 0.8
        #
        x_vec = np.linspace(xmin - h_res, xmax + h_res, int((xmax-xmin + 2*h_res)/h_res))
        y_vec = np.linspace(ymin - h_res, ymax + h_res, int((ymax-ymin + 2*h_res)/h_res))
        z_vec = np.linspace(zmax, zmin, n_layer)
        self.layer_elev = np.around(z_vec, 1)
        #
        grid_2d=np.array([[i,j] for i in x_vec for j in y_vec])
        #
        unstru_grid_2d = np.zeros(grid_2d.shape)
        i_domain = np.zeros(grid_2d.shape[0]).astype('int64')
        all_pt = np.vstack([self.flat_top_surface, self.flat_bot_surface])
        for i in range(len(grid_2d)):
            x = grid_2d[i][0]
            y = grid_2d[i][1]
            dis = (all_pt[:, 1] - x)**2 + (all_pt[:, 2] - y)**2
            if min(dis) < (model_res)**2:
                unstru_grid_2d[i] = grid_2d[i].copy()
                i_domain[i] = 1
        #
        unstru_grid_2d = unstru_grid_2d[i_domain == 1]
        unstru_grid_2d = np.around(unstru_grid_2d, 2)
        self.grid2d_nondip = np.hstack([np.arange(0, len(unstru_grid_2d)).reshape(-1, 1), unstru_grid_2d]) # 2Dgrid ID, X, Y
        #
        return 1
    #
    # quasi-3D lithologic modeling
    def model_3D(self, method):
        litho_layer_all = np.zeros((len(self.grid2d_nondip), 5, len(self.layer_elev)))
        litho_layer_all[:] = np.nan
        #
        for i in range(len(self.layer_elev)):
            if (i+1) % 10 == 0:
                print(i+1, 'out of ', len(self.layer_elev), 'layers.')
            #
            elev = self.layer_elev[i]
            #
            if method == 'Natural_neighbor':
                litho_layer = self.nn_interpo_layer(m = elev, grid2d_nondip = self.grid2d_nondip)
            elif method == 'Nearest_neighbor':
                litho_layer = self.nearest_interpo_layer(m = elev, grid2d_nondip = self.grid2d_nondip)
            #
            litho_layer_all[:, :, i] = litho_layer
        #
        litho_upscale = []
        for i in range(len(self.grid2d_nondip)):
            geo_col = litho_layer_all[i, :, :].T
            #
            assert all(p == i for p in geo_col[:, 0]) == 1
            assert len(np.unique(geo_col[:, 1])) == 1
            assert len(np.unique(geo_col[:, 2])) == 1
            assert all(m > n for m,n in zip(geo_col[0:-1, 3], geo_col[1::, 3])) == 1
            #
            col_upscale = self.upscale_a_column(geo_col = geo_col)
            litho_upscale.append(col_upscale)
        #
        litho_upscale = np.vstack(litho_upscale)
        #    
        return litho_upscale
    #
    # Natural Neighbor interpolation for one layer
    def nn_interpo_layer(self, m, grid2d_nondip): # m is elevation of the layer
        sample_pt = self.get_layer_data(elev=m)
        #
        tri = Delaunay(np.array([p[1:3] for p in sample_pt]))
        members, circumcenters = geometry.find_natural_neighbors(tri, [tuple(p[1::]) for p in grid2d_nondip.tolist()])
        #
        xp = sample_pt[:, 1].astype('float64')
        yp = sample_pt[:, 2].astype('float64')
        zp = sample_pt[:, 4].astype('float64')
        #
        members_list=[[key,value] for key, value in members.items()]
        #
        nn_layer_results = np.zeros((len(grid2d_nondip), 5))
        nn_layer_results[:, 0:3] = grid2d_nondip
        nn_layer_results[:, 3] = m
        nn_layer_results[:, 4] = -1
        #
        i_nn = np.zeros(len(grid2d_nondip)).astype('int64') # whether natural neighbor can be used at a point
        #
        for i in range(len(grid2d_nondip)):
            if len(members_list[i][1])>0:
                val = natural_neighbor_point(xp, yp, zp, (grid2d_nondip[i][1], grid2d_nondip[i][2]),
                                             tri, members_list[i][1],circumcenters)
                #
                if val>=0.5:
                    interpo_class=1
                else:
                    interpo_class=0
                #
                nn_layer_results[i, 4] = interpo_class # Grid2D X Y Z indicator
                #
                i_nn[i] = 1
            else:
                i_nn[i] = 0
        #
        # Nearest neighbor interpolation to handle extrapolation
        nn_extra = self.nearest_interpo(sample_info = nn_layer_results[i_nn == 1], 
                                        query_info = nn_layer_results[i_nn == 0])
        nn_layer_results[i_nn == 0] = nn_extra
        #
        return nn_layer_results
    #
    # Nearest neighbor interpolation for one layer
    def nearest_interpo_layer(self, m, grid2d_nondip): # m is elevation of the layer
        sample_pt = self.get_layer_data(elev=m)
        #
        sample_info = np.array([p[:] for p in sample_pt])
        sample_info[:, 0] = np.nan
        sample_info = sample_info.astype('float64')
        #
        query_info = np.zeros((len(grid2d_nondip), 5))
        query_info[:, 0:3] = grid2d_nondip
        query_info[:, 3] = m
        query_info[:, 4] = -1
        #
        nn_layer_results = self.nearest_interpo(sample_info = sample_info, 
                                                query_info = query_info)
        #
        return nn_layer_results
    #
    # Get point data at each log given elevation of a layer
    def get_layer_data(self, elev):
        layer_data = np.empty((len(self.lognames), 5), dtype = 'U64')
        i_ava = np.zeros(len(self.lognames)).astype('int64') # whether log data is available on a layer
        #
        for i in range(len(self.lognames)):
            logname = self.lognames[i]
            indx = np.where(self.log_nondip[:, 0] == logname)[0]
            onelog = self.log_nondip[indx]
            #
            arr = onelog[:, 1:4].astype('float64')
            elev_seq = onelog[:, 3].astype('float64')
            soil_seq = onelog[:, 4]
            #
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
                #
                layer_data[i][0] = logname
                layer_data[i][1:4] = log_point.copy()
                layer_data[i][4] = indicator
                i_ava[i] = 1
        #
        layer_data = layer_data[i_ava == 1]
        #
        return layer_data
    #
    # Nearest neighbor interpolation
    def nearest_interpo(self, sample_info, query_info):
        x = sample_info[:,1].copy()
        y = sample_info[:,2].copy()
        z = sample_info[:,4].copy()
        #
        xq = query_info[:,1].copy()
        yq = query_info[:,2].copy()
        zq = griddata((x,y),z,(xq,yq),method='nearest')
        #
        zq = np.where(zq >= 0.5, 1, 0)
        #
        result = np.array([p[:] for p in query_info])
        result[:, 4] = zq
        #
        return result
    #
    # Upscale a column of lithology
    def upscale_a_column(self, geo_col):
        col_upscale = np.zeros(geo_col.shape)
        col_upscale[:] = np.nan
        #
        t = 1
        col_upscale[0] = geo_col[0].copy()
        #
        for i in range(1, len(geo_col)):
            if geo_col[i, 4] != col_upscale[t-1, 4]:
                col_upscale[t] = geo_col[i].copy()
                t += 1
        #
        if col_upscale[t-1, 3] != 10.0:
            col_upscale[t, 0:3] = geo_col[0, 0:3].copy()
            col_upscale[t, 3] = 10.0
            col_upscale[t, 4] = -999
            t += 1
        else:
            col_upscale[t-1, 4] = -999
        #
        col_upscale = col_upscale[0:t]
        #
        return col_upscale
    #

#
#
# Module for transforming non-dipping model to dipping domain
class Back_transform():
    def __init__(self, model_nondip, grid2d_nondip, grid2d_dip,
                 faces, top_surface_path, bot_surface_path,
                 flat_top_surface_path, flat_bot_surface_path,
                 NHG = False, vtu2d_path = None):
        
        self.model_nondip = np.array([p[:] for p in model_nondip])
        self.grid2d_nondip = np.array([p[:] for p in grid2d_nondip])
        
        self.faces = np.array([p[:] for p in faces]).astype('int64')
        
        self.top_surface = np.loadtxt(top_surface_path)
        self.top_surface = np.hstack((np.arange(len(self.top_surface)).reshape(-1, 1), self.top_surface))
        self.bot_surface = np.loadtxt(bot_surface_path)
        self.bot_surface = np.hstack((np.arange(len(self.bot_surface)).reshape(-1, 1), self.bot_surface))

        self.flat_top_surface = np.loadtxt(flat_top_surface_path)
        self.flat_top_surface = self.flat_top_surface[self.flat_top_surface[:,0].argsort()]
        self.flat_top_surface[:,-1] = 300. #Set flat top surface elevation_meter
        self.flat_bot_surface = np.loadtxt(flat_bot_surface_path)
        self.flat_bot_surface = self.flat_bot_surface[self.flat_bot_surface[:,0].argsort()]
        self.flat_bot_surface[:,-1] = 10. #Set flat bottom surface elevation_meter
        
        #
        if NHG == True:
            # For model built on NHG grid:
            # "grid2d_dip" is filled with gridID following by two vertexID of surface
            # the two vertexID locate upperleft vertex and lowerright vertex of a NHG cell
            self.grid2d_dip = self.get_grid2d_dip_NHG(vtu2d_path = vtu2d_path)
        else:
            self.grid2d_dip = np.array([p[:] for p in grid2d_dip])
        
        #
        self.organized_model_nondip = []
        cell2d = np.unique(self.model_nondip[:, 0])
        for i in cell2d:
            ind = np.where(self.model_nondip[:, 0] == i)[0]
            self.organized_model_nondip.append(self.model_nondip[ind])
        
        #
        
    #
    # Transform all columns
    def transform_columns(self, vertsol):
        model_dip = []
        trans_flag = np.zeros(len(self.grid2d_dip)).astype('int64') # whether a column is transformed
        #
        for i in range(len(self.grid2d_dip)):
            if i % 1000 == 0:
                print(i)
            #
            agrid = self.grid2d_dip[i]
            i_transform, acol = self.transform_a_column(agrid = agrid, vertsol = vertsol)
            model_dip.append(acol)
            trans_flag[i] = i_transform
        #
        print(sum(trans_flag == 0), 'out of ', 
              len(self.grid2d_dip), 'column not transformed.')
        #
        query_grid = self.grid2d_dip[trans_flag == 0]
        sample_grid = self.grid2d_dip[trans_flag == 1]
        for i in query_grid:
            cid = int(i[0]) - 1
            x = i[1]
            y = i[2]
            dist2 = (x - sample_grid[:, 1])**2 + (y - sample_grid[:, 2])**2
            ind = np.where(dist2 == min(dist2))[0]
            nearest_cid = int(sample_grid[ind[0]][0]) - 1
            #
            nearest_col = model_dip[nearest_cid]
            cp_col = np.array([p[:] for p in nearest_col])
            cp_col[:, 0] = i[0]
            cp_col[:, 1] = x
            cp_col[:, 2] = y
            #
            model_dip[cid] = np.array([p[:] for p in cp_col])
        #
        return trans_flag, model_dip
    #
    #
    # Transform for one column
    def transform_a_column(self, agrid, vertsol):
        col_id = agrid[0]
        col_xy = agrid[1::]
        i_find, fall_face_vertex, fall_edge_vertex, fc = self.locate_face_edge(col_xy = col_xy)
        #
        if i_find == 1:
            fall_face_top = self.top_surface[fall_face_vertex, 1::] # Top of dipped face where a log falls
            pt_top = self.locate_dip_surface_pt(tri_pt = fall_face_top, col_xy = col_xy) # intersect point of log and the top face
            fall_face_bot = self.bot_surface[fall_face_vertex, 1::] # Bottom of dipped face where a log falls
            pt_bot = self.locate_dip_surface_pt(tri_pt = fall_face_bot, col_xy = col_xy) # intersect point of log and the bottom face
            #
            assert pt_top[-1] > pt_bot[-1]
            #
            fall_face_flat_top = self.flat_top_surface[fall_face_vertex, 1::] # Top of flat face where a log falls
            fall_face_flat_bot = self.flat_bot_surface[fall_face_vertex, 1::] # Bottom of flat face where a log falls
            #
            intersect_top=self.find_intersections_v2(dip_face=fall_face_top,
                                                  flat_face=fall_face_flat_top, pt=pt_top)
            intersect_bot=self.find_intersections_v2(dip_face=fall_face_bot,
                                                  flat_face=fall_face_flat_bot, pt=pt_bot)
            #
            status, pt_top_flat = self.select_intersect(intersect=intersect_top,
                                                        flat_face=fall_face_flat_top) # log location in top flat surface
            status, pt_bot_flat = self.select_intersect(intersect=intersect_bot,
                                                        flat_face=fall_face_flat_bot) # log location in bottom flat surface
            #
            coords, tar_grid_nondip, col = self.linear_transform(dip_top = pt_top, dip_bot = pt_bot,
                                        nondip_top = pt_top_flat, nondip_bot = pt_bot_flat,
                                        vertsol = vertsol)
            i_transform = 1
            #
            col_upscale = self.upscale_a_column_dip(col = col)
            #
            col_upscale_addid = np.zeros((len(col_upscale), 5))
            col_upscale_addid[:, 1::] = col_upscale
            col_upscale_addid[:, 0] = col_id
        #
        elif i_find == 2:
            fall_edge_top = self.top_surface[fall_edge_vertex, 1::] 
            pt_top = fall_edge_top[0] + fc * (fall_edge_top[1] - fall_edge_top[0])
            pt_top = np.around(pt_top, 2)
            fall_edge_bot = self.bot_surface[fall_edge_vertex, 1::]
            pt_bot = fall_edge_bot[0] + fc * (fall_edge_bot[1] - fall_edge_bot[0])
            pt_bot = np.around(pt_bot, 2)
            #
            assert pt_top[-1] > pt_bot[-1]
            #
            fall_edge_flat_top = self.flat_top_surface[fall_edge_vertex, 1::]
            fall_edge_flat_bot = self.flat_bot_surface[fall_edge_vertex, 1::]
            #
            pt_top_flat = fall_edge_flat_top[0] + fc * (fall_edge_flat_top[1] - fall_edge_flat_top[0])
            pt_bot_flat = fall_edge_flat_bot[0] + fc * (fall_edge_flat_bot[1] - fall_edge_flat_bot[0])
            pt_top_flat = np.around(pt_top_flat, 2)
            pt_bot_flat = np.around(pt_bot_flat, 2)
            #
            coords, tar_grid_nondip, col = self.linear_transform(dip_top = pt_top, dip_bot = pt_bot,
                                        nondip_top = pt_top_flat, nondip_bot = pt_bot_flat,
                                        vertsol = vertsol)
            i_transform = 1
            #
            col_upscale = self.upscale_a_column_dip(col = col)
            #
            col_upscale_addid = np.zeros((len(col_upscale), 5))
            col_upscale_addid[:, 1::] = col_upscale
            col_upscale_addid[:, 0] = col_id
        else:
            coords = tar_grid_nondip = col = col_upscale_addid = None
            i_transform = 0
        #
        return i_transform, col_upscale_addid
        #return i_transform, coords, tar_grid_nondip, col, col_upscale_addid
    #
    # find a face where a point falls in
    def locate_face(self, col_xy):
        i_find = 0 # whether a face is located that enclose the point
        fall_face_vertex = np.ones(3).astype('int64') # vertex ID of the located face
        fall_face_vertex[:] = -1
        for i in range(len(self.faces)):
            tri_pt = self.top_surface[self.faces[i][1::]]
            tri_pt = tri_pt[:, 1::]
            i_find = ptintri(tri2d = tri_pt, pt = col_xy)
            if i_find == 1:
                fall_face_vertex = self.faces[i][1::].copy()
                break
        #
        return i_find, fall_face_vertex
    #
    # find a face or a edge where a point falls
    def locate_face_edge(self, col_xy):
        i_find = 0 # whether a face or edge is located that enclose the point
                   # 1 means the point falls in a face
                   # 2 means the point falls on a edge
        fall_face_vertex = np.ones(3).astype('int64') # vertex ID of the located face
        fall_edge_vertex = np.ones(2).astype('int64') # vertex ID of the located edge
        fall_face_vertex[:] = -1
        fall_edge_vertex[:] = -1
        fc = None # a factor to locate the point on edge
        #
        x = col_xy[0]
        y = col_xy[1]
        dist2 = (self.top_surface[:, 1] - x)**2 + (self.top_surface[:, 2] - y)**2
        ind = np.where(dist2 == np.min(dist2))[0]
        near_ptid = ind[0] # nearest node id
        #
        inde = np.where(self.faces[:, 1::] == near_ptid)[0]
        near_faces = self.faces[inde]
        # try to search faces
        for i in range(len(near_faces)):
            tri_pt = self.top_surface[near_faces[i][1::]]
            tri_pt = tri_pt[:, 1::]
            i_find = ptintri(tri2d = tri_pt, pt = col_xy)
            if i_find == 1:
                fall_face_vertex = near_faces[i][1::].copy()
                break
        #
        # try to search edges
        if i_find == 0:
            med1 = near_faces[:, [1, 2]]
            med2 = near_faces[:, [1, 3]]
            med3 = near_faces[:, [2, 3]]
            pt_pairs = np.vstack([med1, med2, med3])
            #
            for pt_pair in pt_pairs:
                line_pt = self.top_surface[pt_pair]
                line_pt = line_pt[:, 1::]
                if self.ptinline(line_pt = line_pt, pt = col_xy) == 1:
                    i_find = 2
                    fall_edge_vertex[:] = pt_pair.copy()
                    #
                    line2d_len = np.sum((line_pt[0][0:2] - line_pt[1][0:2]) ** 2) ** 0.5
                    dd = np.sum((line_pt[0][0:2] - col_xy) ** 2) ** 0.5
                    fc = dd / line2d_len
                    break
        #
        return i_find, fall_face_vertex, fall_edge_vertex, fc
    #
    # judge whether a point is on a line
    def ptinline(self, line_pt, pt):
        xx = line_pt[:, 0]
        yy = line_pt[:, 1]
        pt_x = pt[0]
        pt_y = pt[1]
        judge = 0
        toler = 0.01 # tolerance 
        #
        if np.abs(xx[1] - xx[0]) < 1e-10:
            # the line is vertical
            if pt_y >= np.min(yy) and pt_y <= np.max(yy) and np.abs(pt_x - xx.mean()) <= toler:
                judge = 1
        else:
            # Line function is y = {m}x + {c}
            m = (yy[1] - yy[0]) / (xx[1] - xx[0]) # slope of line
            c = yy[0] - m * xx[0] # intercept with y axis
            y_bar = m * pt_x + c
            if pt_y >= np.min(yy) and pt_y <= np.max(yy) and pt_x >= np.min(xx) and pt_x <= np.max(xx) and np.abs(y_bar - pt_y) <= toler:
                judge = 1
        #
        return judge
    #
    # locate intersection point on dipping face
    def locate_dip_surface_pt(self, tri_pt, col_xy):
        x = tri_pt[:, 0]
        y = tri_pt[:, 1]
        z = tri_pt[:, 2]
        #
        xq = col_xy[0]
        yq = col_xy[1]
        zq = griddata((x,y),z,(xq,yq), method='linear')
        #
        result = np.array([xq, yq, zq])
        result = np.around(result, 2)
        #
        return result
    #
    # find candidate intersection points on flat face
    def find_intersections(self, dip_face, flat_face, pt):
        d1,d2 = self.calculate_distance(dip_face, flat_face, pt)
        #
        m1,c1 = self.get_line_function(points = np.vstack([flat_face[0][0:2], flat_face[1][0:2]]))
        m2,c2 = self.get_line_function(points = np.vstack([flat_face[0][0:2], flat_face[2][0:2]]))
        angle1=math.atan(m1)
        angle2=math.atan(m2)
        diff1=abs(d1/math.cos(angle1))
        diff2=abs(d2/math.cos(angle2))
        intersecion = np.zeros((4, 2))
        intersecion[:] = np.nan
        intersecion[0] = self.calculate_intersection(slope1=m1,slope2=m2,intercept1=c1+diff1,intercept2=c2+diff2)
        intersecion[1] = self.calculate_intersection(slope1=m1,slope2=m2,intercept1=c1+diff1,intercept2=c2-diff2)
        intersecion[2] = self.calculate_intersection(slope1=m1,slope2=m2,intercept1=c1-diff1,intercept2=c2+diff2)
        intersecion[3] = self.calculate_intersection(slope1=m1,slope2=m2,intercept1=c1-diff1,intercept2=c2-diff2)
        #
        return intersecion
    #
    # more stable verison
    def find_intersections_v2(self, dip_face, flat_face, pt):
        d1,d2 = self.calculate_distance(dip_face, flat_face, pt)
        #
        m1,c1 = self.get_line_function_v2(points = np.vstack([flat_face[0][0:2], flat_face[1][0:2]]))
        m2,c2 = self.get_line_function_v2(points = np.vstack([flat_face[0][0:2], flat_face[2][0:2]]))
        #
        if m1 == 'inf':
            diff1 = d1
        else:
            angle1=math.atan(m1)
            diff1=abs(d1/math.cos(angle1))
        if m2 == 'inf':
            diff2 = d2
        else:            
            angle2=math.atan(m2)
            diff2=abs(d2/math.cos(angle2))
        #
        intersecion = np.zeros((4, 2))
        intersecion[:] = np.nan
        intersecion[0] = self.calculate_intersection_v2(slope1=m1,slope2=m2,intercept1=c1+diff1,intercept2=c2+diff2)
        intersecion[1] = self.calculate_intersection_v2(slope1=m1,slope2=m2,intercept1=c1+diff1,intercept2=c2-diff2)
        intersecion[2] = self.calculate_intersection_v2(slope1=m1,slope2=m2,intercept1=c1-diff1,intercept2=c2+diff2)
        intersecion[3] = self.calculate_intersection_v2(slope1=m1,slope2=m2,intercept1=c1-diff1,intercept2=c2-diff2)
        #
        return intersecion
    #
    #
    def calculate_distance(self, dip_face, flat_face, pt):
        s1,s2,s3 = self.calculate_area(vertex_list=dip_face, pt=pt)
        fall_face_flat_area = self.tri3d_area(tri_list = flat_face)
        #
        ls1=np.sqrt(sum((flat_face[0] - flat_face[1])**2)) # lenth of the first side of flat triangle
        ls2=np.sqrt(sum((flat_face[0] - flat_face[2])**2)) # second side
        ls3=np.sqrt(sum((flat_face[1] - flat_face[2])**2)) # third side
        #
        d1=2*fall_face_flat_area*s1/(s1+s2+s3)/ls1
        d2=2*fall_face_flat_area*s2/(s1+s2+s3)/ls2
        #d3=2*fall_face_flat_area*s3/(s1+s2+s3)/ls3
        #
        return d1,d2
    #
    #
    def calculate_area(self, vertex_list, pt):
        tri1=[pt, vertex_list[0], vertex_list[1]]
        tri2=[pt, vertex_list[0], vertex_list[2]]
        tri3=[pt, vertex_list[1], vertex_list[2]]
        s1=self.tri3d_area(tri_list=tri1)
        s2=self.tri3d_area(tri_list=tri2)
        s3=self.tri3d_area(tri_list=tri3)
        #
        return s1,s2,s3
    #
    #
    def tri3d_area(self, tri_list):
        verte=np.array(tri_list)
        side1=verte[0]-verte[1]
        side2=verte[0]-verte[2]
        cos_angle=np.dot(side1,side2)/(np.sqrt(sum(side1**2))*np.sqrt(sum(side2**2)))
        area=0.5*np.sqrt(sum(side1**2))*np.sqrt(sum(side2**2))*(1-cos_angle**2)**0.5
        #
        return area
    #
    #
    def get_line_function(self, points):
        pt_arr = np.array(points)
        x_coord = pt_arr[:, 0]
        y_coord = pt_arr[:, 1]
        if np.abs(x_coord[1] - x_coord[0]) < 1e-20:
            m = (y_coord[1] - y_coord[0]) / (1e-20) # slope of line
        else:
            m = (y_coord[1] - y_coord[0]) / (x_coord[1] - x_coord[0]) # slope of line
        c = y_coord[0] - m * x_coord[0] # intercept with y axis
        # Line function is y = {m}x + {c}
        #
        return m, c
    #
    # more stable version
    def get_line_function_v2(self, points):
        pt_arr = np.array(points)
        x_coord = pt_arr[:, 0]
        y_coord = pt_arr[:, 1]
        if np.abs(x_coord[1] - x_coord[0]) < 1e-10:
            m = 'inf' # line is vertical
            c = (x_coord[0] + x_coord[1]) / 2
            # Line function is x = {a}, a is denoted as c
        else:
            m = (y_coord[1] - y_coord[0]) / (x_coord[1] - x_coord[0]) # slope of line
            c = y_coord[0] - m * x_coord[0] # intercept with y axis
            # Line function is y = {m}x + {c}
        #
        return m, c
    #
    #
    def calculate_intersection(self, slope1,slope2,intercept1,intercept2):
        x=(intercept2-intercept1)/(slope1-slope2)
        y=slope1*x+intercept1
        return np.array([x,y])
    #
    # more stable version
    def calculate_intersection_v2(self, slope1,slope2,intercept1,intercept2):
        if slope1 != 'inf' and slope2 != 'inf':
            x=(intercept2-intercept1)/(slope1-slope2)
            y=slope1*x+intercept1
        elif slope1 == 'inf':
            x = intercept1
            y = slope2 * x + intercept2
        elif slope2 == 'inf':
            x = intercept2
            y = slope1 * x + intercept1
        #
        return np.array([x,y])
    #
    # select a intersection point with flat face among 4 candidate points
    def select_intersect(self, intersect, flat_face):
        pt_on_flat = np.zeros(3)
        pt_on_flat[:] = np.nan
        count = 0
        for i in range(len(intersect)):
            judge = ptintri(tri2d = flat_face[:, 0:2], pt = intersect[i])
            if judge == 1:
                pt_on_flat[0:2] = intersect[i].copy()
                pt_on_flat[2] = flat_face[0, 2]
                pt_on_flat = np.around(pt_on_flat, 2)
                count += 1
        #
        if count == 1:
            status = 'Normal'
        else:
            status = 'Problem'
        #
        return status, pt_on_flat
    #
    # Linearly transform a column of geology
    def linear_transform(self, dip_top, dip_bot, nondip_top, nondip_bot, vertsol):
        zz = np.linspace(dip_top[-1], dip_bot[-1], int((dip_top[-1] - dip_bot[-1]) / vertsol))
        zz = np.around(zz, 2)
        col = np.zeros((len(zz), 4)) # X, Y, Z, Litho
        col[:] = np.nan
        col[:, 0] = dip_top[0]
        col[:, 1] = dip_top[1]
        col[:, 2] = zz.copy()
        #
        t = (zz[0] - zz) / (zz[0] - zz[-1])
        v = nondip_bot - nondip_top
        coords = nondip_top + t.reshape(-1, 1) @ v.reshape(1, -1) # correspongding X Y Z in non-dipping domain
        #
        tar_grid_nondip = self.search_grid2d(coords = coords) # locate grid2D ID for each point in non-dipping domain
        #
        for i in range(len(coords)):
            elev = coords[i, 2]
            cid = int(tar_grid_nondip[i])
            col_nondip = self.organized_model_nondip[cid] # column in nondipping domain to find lithology
            litho = self.get_litho(elev, col_nondip)
            #
            col[i, 3] = litho
        #
        return coords, tar_grid_nondip, col
    #
    #
    def search_grid2d(self, coords):
        grid_arr = np.zeros(len(coords))
        grid_arr[:] = np.nan
        rad_factor = 2.5
        #print(coords)
        # locate top point
        x = coords[0, 0]
        y = coords[0, 1]
        dist2 = (x - self.grid2d_nondip[:, 1]) ** 2 + (y - self.grid2d_nondip[:, 2]) ** 2
        ind = np.where(dist2 == min(dist2))[0]
        grid_arr[0] = self.grid2d_nondip[ind[0]][0]
        #
        # locate bottom point
        x = coords[-1, 0]
        y = coords[-1, 1]
        dist2 = (x - self.grid2d_nondip[:, 1]) ** 2 + (y - self.grid2d_nondip[:, 2]) ** 2
        ind = np.where(dist2 == min(dist2))[0]
        grid_arr[-1] = self.grid2d_nondip[ind[0]][0]
        #
        if grid_arr[0] == grid_arr[-1]:
            grid_arr[:] = grid_arr[0]
        #
        else:
            top_id = int(grid_arr[0])
            bot_id = int(grid_arr[-1])
            rad = np.sqrt((self.grid2d_nondip[bot_id, 1] - self.grid2d_nondip[top_id, 1]) ** 2 + 
                         (self.grid2d_nondip[bot_id, 2] - self.grid2d_nondip[top_id, 2]) ** 2)
            #
            x = self.grid2d_nondip[top_id, 1]
            y = self.grid2d_nondip[top_id, 2]
            dd_1 = np.sqrt((x - self.grid2d_nondip[:, 1]) ** 2 + (y - self.grid2d_nondip[:, 2]) ** 2)
            #
            x = self.grid2d_nondip[bot_id, 1]
            y = self.grid2d_nondip[bot_id, 2]
            dd_2 = np.sqrt((x - self.grid2d_nondip[:, 1]) ** 2 + (y - self.grid2d_nondip[:, 2]) ** 2)
            #
            inde = np.where((dd_1 <= rad_factor * rad) & (dd_2 <= rad_factor * rad))[0]
            sub = self.grid2d_nondip[inde]
            #
            for i in range(1, len(coords) - 1):
                x = coords[i, 0]
                y = coords[i, 1]
                #dist2 = (x - self.grid2d_nondip[:, 1]) ** 2 + (y - self.grid2d_nondip[:, 2]) ** 2
                dist2 = (x - sub[:, 1]) ** 2 + (y - sub[:, 2]) ** 2
                ind = np.where(dist2 == min(dist2))[0]
                grid_arr[i] = sub[ind[0]][0]
        #
        return grid_arr
    #
    #
    def get_litho(self, elev, col_nondip):
        elev_seq = col_nondip[:, 1].copy()
        ind = np.where(elev <= elev_seq)[0]
        litho = col_nondip[ind[-1]][2]
        if litho == -999:
            litho = col_nondip[ind[-1] - 1][2]
        #
        return litho
    #
    # Upscale a column in dipping domain
    def upscale_a_column_dip(self, col):
        bot = col[-1][2]
        col_upscale = np.zeros(col.shape)
        col_upscale[:] = np.nan
        #
        t = 1
        col_upscale[0] = col[0].copy()
        #
        for i in range(1, len(col)):
            if col[i, 3] != col_upscale[t-1, 3]:
                col_upscale[t] = col[i].copy()
                t += 1
        #
        if col_upscale[t-1, 2] != bot:
            col_upscale[t, 0:2] = col[0, 0:2].copy()
            col_upscale[t, 2] = bot
            col_upscale[t, 3] = -999
            t += 1
        else:
            col_upscale[t-1, 3] = -999
        #
        col_upscale = col_upscale[0:t]
        #
        return col_upscale
    #
    # Prepare grid2d information in dipping domain on NHG
    def get_grid2d_dip_NHG(self, vtu2d_path):
        vertex_id = self.get_NHG_vertex(vtu2d_path = vtu2d_path)
        cell_center = self.get_cell_center_xy(vertex_id = vertex_id)
        grid2d_dip = np.zeros((len(vertex_id), 3)).astype('int64')
        #
        for i in range(len(vertex_id)):
            center = cell_center[i]
            pt_id = vertex_id[i]
            pt_coords = self.top_surface[pt_id]
            #
            center_x = center[0]
            center_y = center[1]
            ind = np.where((pt_coords[:, 1] < center_x) & (pt_coords[:, 2] < center_y))[0]
            lowleft_id = pt_id[ind[0]]
            ind = np.where((pt_coords[:, 1] > center_x) & (pt_coords[:, 2] > center_y))[0]
            upright_id = pt_id[ind[0]]
            #
            remain_id = np.array([p for p in pt_id if p not in [lowleft_id, upright_id]])
            #
            grid2d_dip[i, 0] = i+1
            grid2d_dip[i, 1::] = remain_id
        #
        return grid2d_dip
    #
    # Get vertex id that defines each NHG cell
    def get_NHG_vertex(self, vtu2d_path):
        # Read 2D VTU
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtu2d_path)
        reader.Update()  # Needed because of GetScalarRange
        output = reader.GetOutput()
        #
        vertex_id = self.get_vertex_id(vtu_obj=output) # Vertex ids defining each 2D cell
        
        assert all(len(p) == 4 for p in vertex_id) == 1
        vertex_id = np.array(vertex_id)
        
        #
        return vertex_id
    #
    # Get vertex id that defines each cell
    def get_vertex_id(self, vtu_obj):
        id_ls=[]
        numberOfCells=vtu_obj.GetNumberOfCells()
        #
        for cellIndex in range(numberOfCells): # for every cell
            cellIds = vtk.vtkIdList() # cell ids store to
            vtu_obj.GetCellPoints(cellIndex, cellIds) # get ids of points of the given cell
            id_ls.append([cellIds.GetId(i) for i in range(0, cellIds.GetNumberOfIds())])
        #
        return id_ls
    #
    # calculate 2D cell center given vertex XY
    def get_cell_center_xy(self, vertex_id):
        cellxy = np.zeros((len(vertex_id), 2))
        for i in range(len(vertex_id)):
            ind = vertex_id[i]
            pts = self.top_surface[ind]
            x = np.mean(pts[:, 1])
            y = np.mean(pts[:, 2])
            cellxy[i] = x, y
        #
        return cellxy
    #
    # Transform all NHG columns
    def transform_NHG_columns(self, vertsol):
        model_dip = []
        #
        for i in range(len(self.grid2d_dip)):
            if i % 1000 == 0:
                print(i)
            #
            agrid = self.grid2d_dip[i]
            acol = self.transform_a_NHG_column(agrid = agrid, vertsol = vertsol)
            model_dip.append(acol)
        #
        return model_dip
    #
    # Transform for a column defined on NHG
    def transform_a_NHG_column(self, agrid, vertsol):
        col_id = agrid[0]
        v_id = agrid[1::]
        #
        pt_top = self.top_surface[v_id, 1::].mean(axis = 0) # NHG center in the top face
        pt_bot = self.bot_surface[v_id, 1::].mean(axis = 0) # NHG center in the bottom face
        #
        pt_top_flat = self.flat_top_surface[v_id, 1::].mean(axis = 0) # corresponding point in top flat surface
        pt_bot_flat = self.flat_bot_surface[v_id, 1::].mean(axis = 0) # corresponding point in bottom flat surface
        #
        coords, tar_grid_nondip, col = self.linear_transform(dip_top = pt_top, dip_bot = pt_bot,
                                                             nondip_top = pt_top_flat, nondip_bot = pt_bot_flat,
                                                             vertsol = vertsol)
        #
        col_upscale = self.upscale_a_column_dip(col = col)
        #
        col_upscale_addid = np.zeros((len(col_upscale), 5))
        col_upscale_addid[:, 1::] = col_upscale
        col_upscale_addid[:, 0] = col_id
        #
        return col_upscale_addid
    #
#
#

#
# A module for constructing top and bottom surfaces of a hydrogeologic unit
# defined on NHG.
class Mesh_generator():
    def __init__(self, vtu2d_path, vtu3d_path):
        
        self.vertex_xy, self.vertex_id = self.get_vertex_data(vtu2d_path = vtu2d_path)
        self.cell_center = self.get_cell_center() # X and Y of 2D cell
        
        self.cell_top, self.cell_bot = self.get_cell_top_bot(vtu3d_path = vtu3d_path)
        assert len(self.cell_center) == len(self.cell_top) == len(self.cell_bot)
        
        self.top_surface_nodes = self.interpo_vertex_z(sample_z = self.cell_top)
        self.bot_surface_nodes = self.interpo_vertex_z(sample_z = self.cell_bot)
        
        self.faces = self.define_faces()
        
        #
    #
    #
    def get_vertex_data(self, vtu2d_path):
        # Read 2D VTU
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtu2d_path)
        reader.Update()  # Needed because of GetScalarRange
        output = reader.GetOutput()
        #
        vertex_xy = np.array(output.GetPoints().GetData()) # Vertex coordinate from 2D VTU
        vertex_id = self.get_cell_vertex_id(vtu_obj=output) # Vertex ids defining each 2D cell
        
        assert all(len(p) == 4 for p in vertex_id) == 1
        vertex_id = np.array(vertex_id)
        
        #
        return vertex_xy, vertex_id
    #
    # Get vertex id that defines each cell
    def get_cell_vertex_id(self, vtu_obj):
        id_ls=[]
        numberOfCells=vtu_obj.GetNumberOfCells()
        #
        for cellIndex in range(numberOfCells): # for every cell
            cellIds = vtk.vtkIdList() # cell ids store to
            vtu_obj.GetCellPoints(cellIndex, cellIds) # get ids of points of the given cell
            id_ls.append([cellIds.GetId(i) for i in range(0, cellIds.GetNumberOfIds())])
        #
        return id_ls
    #
    # calculate 2D cell center given vertex XY
    def get_cell_center(self):
        cellxy = np.zeros((len(self.vertex_id), 2))
        for i in range(len(self.vertex_id)):
            ind = self.vertex_id[i]
            pts = self.vertex_xy[ind]
            x = np.mean(pts[:, 0])
            y = np.mean(pts[:, 1])
            cellxy[i] = x, y
        #
        return cellxy
    #
    # Get cell top and bottom from 3D VTU
    def get_cell_top_bot(self, vtu3d_path):
        # Read 3D VTU
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtu3d_path)
        reader.Update()  # Needed because of GetScalarRange
        output = reader.GetOutput()
        #
        cell_center_top = np.array(output.GetCellData().GetArray("XMS:Cell Top"))
        cell_center_bot = np.array(output.GetCellData().GetArray("XMS:Cell Bottom"))
        #
        return cell_center_top, cell_center_bot
    #
    # Inpolate vertex elevation using cell center data
    def interpo_vertex_z(self, sample_z):
        x = self.cell_center[:, 0].copy()
        y = self.cell_center[:, 1].copy()
        z = sample_z.copy()
        #
        xq = self.vertex_xy[:, 0]
        yq = self.vertex_xy[:, 1]
        #
        zq_linear=griddata((x,y),z,(xq,yq),method='linear')
        zq_near=griddata((x,y),z,(xq,yq),method='nearest')
        zq = np.where(np.isnan(zq_linear) == 1, zq_near, zq_linear)
        #
        results = np.zeros(self.vertex_xy.shape)
        results[:, 0:2] = self.vertex_xy[:, 0:2]
        results[:, 2] = zq
        
        results = np.around(results, 2)
        #
        return results
    #
    # Define faces of triangle meshes
    def define_faces(self):
        faces = np.zeros((len(self.vertex_id) * 2, 4)).astype('int64')
        faces[:, 0] = 3
        #
        for i in range(len(self.vertex_id)):
            center = self.cell_center[i]
            pt_id = self.vertex_id[i]
            pt_coords = self.vertex_xy[pt_id]
            #
            center_x = center[0]
            center_y = center[1]
            ind = np.where((pt_coords[:, 0] < center_x) & (pt_coords[:, 1] < center_y))[0]
            lowleft_id = pt_id[ind[0]]
            ind = np.where((pt_coords[:, 0] > center_x) & (pt_coords[:, 1] > center_y))[0]
            upright_id = pt_id[ind[0]]
            #
            remain_id = np.array([p for p in pt_id if p not in [lowleft_id, upright_id]])
            #
            faces[2*i, 1:3] = remain_id
            faces[2*i, 3] = lowleft_id
            faces[2*i+1, 1:3] = remain_id
            faces[2*i+1, 3] = upright_id
        #
        return faces
    #

#
            
            
            
    