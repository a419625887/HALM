# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 14:21:38 2023

@author: syang34
"""

import numpy as np
import csv
from HALM_lib import Flexural_slip as fs
from HALM_lib import Map_logs as mlgs
from HALM_lib import Litho_modeling as lmg
from HALM_lib import Back_transform as btf

#
# flat top and bottom surfaces
top_surface_path = 'Input/top_surface_point.txt'
bot_surface_path = 'Input/bot_surface_point.txt'
flat_top_surface_path = 'Output/restored_top_point_data.txt'
flat_bot_surface_path = 'Output/restored_bot_point_data.txt'

dip_sur_paths = [top_surface_path, bot_surface_path]
nondip_sur_paths = [flat_top_surface_path, flat_bot_surface_path]

faces = np.loadtxt('Input/triangle_faces.txt').astype('int64')

for i in range(len(dip_sur_paths)):
    
    dip_sur_path = dip_sur_paths[i]
    vertices = np.loadtxt(dip_sur_path)
    
    #
    fs_obj = fs(vertices = vertices, faces = faces)
    
    fs_faces = fs_obj.faces
    fs_vertices = fs_obj.vertices
    fs_norm = fs_obj.norm
    fs_centers = fs_obj.centers
    
    fs_rotated_face = fs_obj.rotated_face
    fs_row_faceid_coll = fs_obj.row_faceid_coll
    
    fs_obj.fit_first_row(startrow = 0)
    fs_obj.fit_remaining_rows()
    
    fs_pt_data = fs_obj.pt_data
    fs_fit_flag = fs_obj.fit_flag
    
    #
    nondip_sur_path = nondip_sur_paths[i]
    np.savetxt(nondip_sur_path, fs_pt_data, fmt = '%d %.2f %.2f %.2f',
               header='Vertex id, X, Y, Z_meter')

#
# map log data
with open('Input/Logs.csv', newline='') as f:
    reader = csv.reader(f,delimiter=',')
    logdata = list(reader)

for count in range(len(logdata)):
    if '' in logdata[count]:
        logdata[count]=logdata[count][0:logdata[count].index('')]

#
mlgs_obj = mlgs(faces = faces, top_surface_path = top_surface_path, bot_surface_path = bot_surface_path,
                flat_top_surface_path = flat_top_surface_path, flat_bot_surface_path = flat_bot_surface_path)

mapped_log_coll = mlgs_obj.map_logs(logdata = logdata)

#
with open('Output/Logs_nondip.csv', "w",newline='') as output:
    writer = csv.writer(output)
    writer.writerows(mapped_log_coll.tolist())


#
# lithological interpolation in non-dipping domain

with open('Output/Logs_nondip.csv', newline='') as f:
    reader = csv.reader(f, delimiter=',')
    log_nondip = list(reader)

lmg_obj = lmg(log_nondip = log_nondip, faces = faces,
              flat_top_surface_path = flat_top_surface_path, flat_bot_surface_path = flat_bot_surface_path)

lmg_obj.generate_2dgrid(model_res = 15.0, n_layer = 300)

lmg_unstru_grid_2d = lmg_obj.grid2d_nondip
lmg_layer_elev = lmg_obj.layer_elev

np.savetxt('Output/Grid2D_nondip.txt', lmg_unstru_grid_2d, fmt = '%d %.2f %.2f')

litho_upscale = lmg_obj.model_3D(method = 'Nearest_neighbor')

geomodel_nondip = np.zeros((len(litho_upscale), 3))
geomodel_nondip[:, 0] = litho_upscale[:, 0] + 1 # make CellID starts from 1
geomodel_nondip[:, 1::] = litho_upscale[:, 3::]

np.savetxt('Output/Geo-model_nondip.dat', litho_upscale, fmt = '%d %.2f %.2f %.2f %d',
           header = 'Cell2D, Z_m, indicator')

#
# Back transform lithofacies to dipping domain
model_nondip = np.array([p[:] for p in geomodel_nondip])

grid2d_dip = np.loadtxt('Input/domain_discretization.txt')
grid2d_dip = grid2d_dip[:, 0:3]

grid2d_nondip = np.array([p[:] for p in lmg_unstru_grid_2d])

btf_obj = btf(model_nondip = model_nondip, grid2d_nondip = grid2d_nondip, grid2d_dip = grid2d_dip,
              faces = faces, top_surface_path = top_surface_path, bot_surface_path = bot_surface_path,
              flat_top_surface_path = flat_top_surface_path, flat_bot_surface_path = flat_bot_surface_path)




organized_model_nondip = btf_obj.organized_model_nondip


trans_flag, model_dip = btf_obj.transform_columns(vertsol = 0.5)

model_dip_w_xy = np.vstack(model_dip)

model_dip_clean = np.array([[p[0], p[3], p[4]] for p in model_dip_w_xy])

np.savetxt('Output/Geo-model.dat', model_dip_w_xy, fmt = '%d %.2f %.2f %.2f %d',
           header = 'Cell2D, Z_m, indicator')