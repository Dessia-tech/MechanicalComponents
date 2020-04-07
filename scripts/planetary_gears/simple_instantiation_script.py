#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 12:43:48 2020

@author: launay
"""

import matplotlib.pyplot as plt
import mechanical_components.planetary_gears as pg
import networkx as nx
import numpy as npy
import dectree 
import math as m
import copy
import volmdlr as vm
import volmdlr.primitives3D as p3d
import volmdlr.primitives2D as p2d
import mechanical_components.meshes as meshes
rack=meshes.Rack(0.34,2)

meshes_1=meshes.Mesh(15,0.06,0.01,rack)
meshes_1.Contour()
speeds = {0: [243.620149773094, 253.56382935566927], 1: [92.0103009101568, 95.76582339628565]}
center_distances=[0.09713462117912072, 0.12713462117912072]
connections = {(0,1):(0,1),(0,1):(0,1)}
torques = {0: -16.380372067375156, 1: 'output'}
mesh_assembly=meshes.MeshCombination(center_distances,connections,{0:meshes_1,1:copy.copy(meshes_1)})
                                  
pos = vm.Point3D((0, 0, 1))
axis = vm.Vector3D((0,0,1))
radius=0.2
length=0.5
cylinder = p3d.Cylinder(pos, axis, radius, length)

volumemodel = vm.Contour2D(meshes_1.Contour(1) )
volumemodel.MPLPlot() 
sun=pg.Planetary('sun',18,'Sun')
ring= pg.Planetary('ring',25,'Sun')
planet_carrier= pg.PlanetCarrier('planet_carrier')
planet_1=pg.Planet('planet_1','Double',10)
planet_2=pg.Planet('planet_2','Double',20)
planet_3=pg.Planet('planet_3','Simple',5)
planet_4=pg.Planet('planet_4','Double',5)
planet_5=pg.Planet('planet_5','Double',5)
planetary_gears_1= pg.PlanetaryGears('pl_1', sun, ring, [planet_1 ,planet_2], planet_carrier)




planetary_gears_2= pg.PlanetaryGears('pl_2', sun, ring, [planet_1], planet_carrier)
planetary_gears_3= pg.PlanetaryGears('pl_3', sun, ring, [planet_1], planet_carrier)
planetary_gears_4= pg.PlanetaryGears('pl_4', sun, ring, [planet_1], planet_carrier)
assembly_planetary_gear=pg.AssemblyPlanetaryGears('assembly_planetary_gear', 
                                                  [planetary_gears_1,planetary_gears_2,planetary_gears_3,planetary_gears_4], 
                                                  [[[sun,planetary_gears_2],[sun,planetary_gears_3]],
                                                    [[planet_carrier,planetary_gears_3],[ring,planetary_gears_4]],
                                                    [[planet_carrier,planetary_gears_1],[planet_carrier,planetary_gears_4]],
                                                    [[sun,planetary_gears_1],[ring,planetary_gears_2]],
                                                    [[ring,planetary_gears_1],[sun,planetary_gears_4]]])

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
# solutions = pg.decision_tree_planetary_gears({'Planetary_1':500,'Planet_Carrier':0},'Planetary_2',-230,2,1,3,[7 , 80],[0.1, 1],0.1)



