#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy

#3 gears meshes test
list_cd = [[0.08,0.12],[0.07,0.12]]
connections = [[(0,1)],[(1,2)]]
list_speed = {0:[1000*npy.pi/30,1200*npy.pi/30],1:[2000*npy.pi/30,
               2100*npy.pi/30],2:[1000*npy.pi/30,1200*npy.pi/30]}
list_rack = {0:{'name':'Racks_A','module':[2*1e-3,2*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi, 20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
list_rack_choice={0 : [0],1 : [0]}
list_helix_angle={0 : [0,0]}
list_material={0 : meshes_opt.hardened_alloy_steel}
list_torque={1 : 186, 0 : 'output', 2 : 20}
list_cycle={1 : 1e12}

GA = meshes_opt.MeshAssemblyOptimizer(Z={},
                                  connections = connections, 
                                  gear_speeds = list_speed,
                                  center_distances = list_cd,
                                  rack_list = list_rack,
                                  rack_choice = list_rack_choice,
                                  helix_angle = list_helix_angle,
                                  material = list_material,torques = list_torque,
                                  cycles = list_cycle, verbose = True)

#Optimization for gear set with center-distance closed to the minimum boundary
GA.OptimizeCD(nb_sol = 2, verbose = True)
print('Number of solutions:',len(GA.solutions))
solution=GA.solutions[-1]
#solution.SVGExport('meshes2.txt',{0 : [0,0], 2 : [0.15,0]})
#solution.FreeCADExport('meshes2')
 
#d = solution.Dict()
#json.dumps(d)
d = solution.Dict()
json_file = json.dumps(d)
d2 = json.loads(json_file)
obj = meshes_opt.MeshAssembly.DictToObject(d2)


obj = meshes_opt.MeshAssembly.DictToObject(d)
#obj.SVGExport('meshes2.txt',{0 : [0,0], 2 : [0.15,0]})
print(obj.mesh_combinations[0].meshes[0] == solution.mesh_combinations[0].meshes[0])
print(obj == solution)
print(hash(obj))
print(hash(solution))
