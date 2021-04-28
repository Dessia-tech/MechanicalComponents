
"""
Created on Tue Apr  6 12:14:22 2021

@author: dasilva
"""

import mechanical_components.optimization.meshes as meshes_opt
import mechanical_components.meshes as me
import numpy as npy

connections = [(0, 1)]

rigid_links = []



gear_speeds = {0: (1878.1453579221634, 1974.460504482274),
               1: (449.8807309958231, 472.95153771355757)}
               
center_dists = [(0.11134984458664793, 0.1293457790652981)]

torques = {0: 'output', 1: 40}

rack = meshes_opt.RackOpti(transverse_pressure_angle_0=[20/180.*npy.pi,20/180.*npy.pi], module=[2*1e-3,2*1e-3],
             coeff_gear_addendum=[1,1],coeff_gear_dedendum=[1.25,1.25],coeff_root_radius=[0.38,0.38],
             coeff_circular_tooth_thickness=[0.5,0.5],helix_angle=[21,60])
        
meshoptis = []
for i, speed_input in enumerate(gear_speeds.values()):
    meshoptis.append(meshes_opt.MeshOpti(rack = rack, torque_input = torques[i], speed_input = speed_input))

center_distances = []
for i , center_distance in enumerate(center_dists):

    center_distances.append(meshes_opt.CenterDistanceOpti(center_distance = center_distance, meshes = [meshoptis[connections[i][0]], meshoptis[connections[i][1]]]))


cycles = {0: 1272321481513.054}

list_rack = {0:{'name':'Catalogue_A','module':[2*1e-3,2*1e-3],
              'transverse_pressure_angle_rack':[20/180.*npy.pi,20/180.*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
rack_choices = {0:[0], 1:[0]}


GA = meshes_opt.MeshAssemblyOptimizer(center_distances,cycles)

#Optimization for gear set with center-distance closed to the minimum boundary
GA.Optimize(nb_sol = 50, verbose=True)
print('Number of solutions:',len(GA.solutions))

solutions = GA.solutions

# solution.pos_axis({0:(0,0,0)})


    
for mesh_assembly in GA.solutions:
    m=mesh_assembly.mesh_combinations
#solution=GA.solutions[-1]
#solution.SVGExport('name.txt',{0 : [0,0], 1 : [0.5,0]})
#solution.FreeCADExport('meshes3')
