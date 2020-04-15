#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 12:43:48 2020

@author: launay
"""

import matplotlib.pyplot as plt
import mechanical_components.planetary_gears_2 as pg
import networkx as nx
import numpy as npy
import dectree 
import math as m
import copy
import volmdlr as vm
import volmdlr.primitives3D as p3d
import volmdlr.primitives2D as p2d
import mechanical_components.meshes as meshes

volumic_mass=7800
data_coeff_YB_Iso={'data':[[0.0,1.0029325508401201],
                           [4.701492563229561,0.9310850480431024],
                           [23.955224059269884,0.7609970656504502],
                           [40.0,0.7492668574805859]
                          ], 'x':'Linear','y':'Linear'}
data_wholer_curve={'data':[[4.307791955971963,1.6419147590563592],
                       [6.518240063668731,1.431665495290182],
                       [7.989456220850952,1.4353220033111185]
                      ], 'x':'Log','y':'Log'}
data_gear_material={'data':[[1.313871566195314,0.7858874572688317],
                      [1.4294457009773085,0.8802021097895326],
                      [1.4551288380965028,0.9097910273994609]
                     ], 'x':'Log','y':'Log'}
material1=meshes.Material(volumic_mass, data_coeff_YB_Iso,
                           data_wholer_curve, data_gear_material)
rack=meshes.Rack(0.34,)

meshes_1=meshes.Mesh(20,0.06,0.01,rack)
meshes_1.Contour()

center_distances=[0.09713462117912072]
connections = [(0,1)]
torques = {0: -16.380372067375156, 1: 'output'}
cycles={0:1e8}
# mesh_assembly=meshes.MeshCombination(center_distances,connections,{0:meshes_1,1:copy.copy(meshes_1)},torques,cycles)
# volumemodel=mesh_assembly.VolumeModel({0:(0,0,0),1:(0,0.2,0.3)})   
# volumemodel.babylonjs()                               
# pos = vm.Point3D((0, 0, 1))
# axis = vm.Vector3D((0,0,1))
# radius=0.2
# length=0.5
# cylinder = p3d.Cylinder(pos, axis, radius, length)

# volumemodel = vm.Contour2D(meshes_1.Contour(1) )
# volumemodel.MPLPlot() 
sun=pg.Planetary('sun',36,'Sun')
sun_2=pg.Planetary('sun_2',60,'Sun')
ring= pg.Planetary('ring',84,'Ring')
planet_carrier= pg.PlanetCarrier('planet_carrier')
planet_1=pg.Planet('planet_1','Simple',12)



planet_2=pg.Planet('planet_2','Simple',12)
planet_3=pg.Planet('planet_3','Simple',12)
planet_4=pg.Planet('planet_4','Double',5)
planet_5=pg.Planet('planet_5','Double',5)
planetary_gears_1= pg.PlanetaryGears('pl_1', [sun,ring,sun_2], [planet_1,planet_2], planet_carrier,[[sun,planet_1,'GE'],[planet_1,planet_2,'GE'],[planet_2,ring,'GE'],[planet_2,sun_2,'GI']])
torque_solution=planetary_gears_1.torque_solve({sun:0,planet_carrier:500})
speed_solution=planetary_gears_1.speed_solve({sun_2:0,planet_carrier:500})
print(torque_solution)
print(speed_solution)

# volume=vm.VolumeModel(planetary_gears_1.volume_plot([0,0],0,[0.02,0.01],0.1,0.5,0.05))
# volume.babylonjs()




# planetary_gears_2= pg.PlanetaryGears('pl_2', [sun,ring], [planet_1], planet_carrier)
# planetary_gears_3= pg.PlanetaryGears('pl_3', sun, ring, [planet_1], planet_carrier)
# planetary_gears_4= pg.PlanetaryGears('pl_4', sun, ring, [planet_1], planet_carrier)
# assembly_planetary_gear=pg.AssemblyPlanetaryGears('assembly_planetary_gear', 
#                                                   [planetary_gears_1,planetary_gears_2,planetary_gears_3,planetary_gears_4], 
#                                                   [[[sun,planetary_gears_2],[sun,planetary_gears_3]],
#                                                     [[planet_carrier,planetary_gears_3],[ring,planetary_gears_4]],
#                                                     [[planet_carrier,planetary_gears_1],[planet_carrier,planetary_gears_4]],
#                                                     [[sun,planetary_gears_1],[ring,planetary_gears_2]],
#                                                     [[ring,planetary_gears_1],[sun,planetary_gears_4]]])

# print(planetary_gears_1)

# planet_carrier=pg.PlanetCarrier('PlanetCarrier')
# planetary_1=pg.Planetary('Planetary_1',7,'Ring')
# planetary_2=pg.Planetary('Planetary_2',7,'Ring')
# list_element={'Planet_Carrier': planet_carrier, 'Planetary_1' : planetary_1,'Planetary_2':planetary_2}
# planets=[]
# for i in range(2):
#             planets.append(pg.Planet('Planet'+str(i),'Double',7)) 
                                                   
# pg.test_ratio_max_ratio_min([2,69,60,42,9],planetary_1,planetary_2,planet_carrier,planets,{'Planetary_1':540,'Planet_Carrier':0},'Planetary_2',200,[7 , 80],3)




# planetary_gears_1.solve(500, ring,planet_3)
# assembly_planetary_gear.plot()
# print(assembly_planetary_gear.system_equations()[0])print(node)
                        
#solution,system_matrix=assembly_planetary_gear.solve({planet_carrier:[500,planetary_gears_2] ,ring:[0,planetary_gears_3],sun:[0,planetary_gears_2]},)

#print(assembly_planetary_gear.solve(500,planet_carrier,planetary_gears_2,[ring,sun],[planetary_gears_3,planetary_gears_2]))

# optimizer=pg.OptimizerPlanetaryGears([[200,250],[100,150], [300,350],[200,300],[200,350],[11,12],[13,15]],0.1,2,45)

# optimizer.decission_tree_()



# solutions=[]
# solutions = pg.decision_tree_planetary_gears({'Planetary_1':500,'Planet_Carrier':0},'Planetary_2',230,2,0.1,3,[7 , 80],[0.1, 1],0.1)

# x = vm.Vector3D((1,0,0))
# y = vm.Vector3D((0,1,0))
# z = vm.Vector3D(npy.cross(x.vector, y.vector))



# Gears3D={0:meshes_1.Contour(3)}
# export=[]
# center_2=(0,0)
# center = vm.Point2D(center_2)
# model_trans =Gears3D[0][0].Translation(center)
# model_trans_rot = model_trans.Rotation(center, 0.1)
# Gears3D_Rotate=[model_trans_rot]
# list_gear=[Gears3D[0]]
# list_center=[center_2]
# list_rot=[-1]
# export=[]
# for (i,center,k) in zip(list_gear,list_center,list_rot):
#             model_export=[]
            
#             for m in i:
                
#                 center = vm.Point2D(center)
#                 model_trans = m.Translation(center)
#                 model_trans_rot = model_trans.Rotation(center, k)
#                 model_export.append(model_trans_rot)
#             export.append(model_export)
# Gears3D_Rotate=export

# vect_x = -0.5*0.01*x 
# extrusion_vector1 = 0.01*x
# C1=vm.Contour2D(Gears3D_Rotate[0])
# i=vm.Vector3D(vect_x)
# t1=p3d.ExtrudedProfile(vm.Vector3D(vect_x), y, z, C1, [], vm.Vector3D(extrusion_vector1))
# modle=vm.VolumeModel([t1],'d')
# modle.babylonjs()

