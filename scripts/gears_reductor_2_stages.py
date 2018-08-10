#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 18:49:12 2018

@author: masfaraud
"""

import mechanical_components.meshes as meshes
import numpy as npy
import itertools as it
import matplotlib.pyplot as plt

"""
Test case to see how results og gears are distributed

"""

nsols = 200

w_error=0.02
list_cd1=[(0.08800844625355042, 0.11800844625355042)]
list_cd2=[(0.08800844625355042, 0.11800844625355042)]

list_gear_set1=[(0,1)]
list_gear_set2=[(0,1)]

list_speed1 = {0: [251.86607746064126, 262.14632552025927], 1: [117.75600124605177, 122.5623686438498]}
list_speed2 = {0: [117.75600124605177, 122.5623686438498], 1: [43.05958236055149, 44.817116334451555]}

list_rack={0:{'name':'Catalogue_A','module':[1*1e-3,3*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}

list_choice={0:[0],1:[0]}

list_torque={0:106,1:'output'}
list_torque={0:106,1:'output'}

MAO1=meshes.MeshAssemblyOptimizer(connections = list_gear_set1, 
                                  gear_speed = list_speed1,
                                  center_distance = list_cd1,
                                  rack_list = list_rack,
                                  torque = list_torque,
                               rack_choice=list_choice)

MAO2=meshes.MeshAssemblyOptimizer(connections = list_gear_set2, 
                                  gear_speed = list_speed2,
                                  center_distance = list_cd2,
                                  rack_list = list_rack,
                                  torque = list_torque,
                               rack_choice=list_choice)


#Recherche triée des nb_sol architecture ayant un entraxe mini (nb_sol=-1 pour analyser l'ensemble des solutions)
MAO1.Optimize(nb_sol=nsols, verbose = True)
MAO2.Optimize(nb_sol=nsols, verbose = True)

ratios1 = []
cd1 = []
w1 = []
for sol in MAO1.solutions:
    Z1 = sol.meshes[0][0].Z
    Z2 = sol.meshes[0][1].Z
    ratios1.append(Z1/Z2)
    cd1.append(0.5*(sol.DF[0][0]+sol.DF[0][1]))
ratio1_min = list_speed1[1][0]/list_speed1[0][1]
ratio1_max = list_speed1[1][1]/list_speed1[0][0]
cd1_min = list_cd1[0][0]
cd1_max = list_cd1[0][1]

ratios2 = []
cd2 = []
for sol in MAO2.solutions:
    Z1 = sol.meshes[0][0].Z
    Z2 = sol.meshes[0][1].Z
    ratios2.append(Z1/Z2)
    cd2.append(0.5*(sol.DF[0][0]+sol.DF[0][1]))
    
ratio2_min = list_speed2[1][0]/list_speed2[0][1]
ratio2_max = list_speed2[1][1]/list_speed2[0][0]
cd2_min = list_cd2[0][0]
cd2_max = list_cd2[0][1]
    
ratios = []
for r1,r2 in it.product(ratios1, ratios2):
    ratios.append(r1 * r2)    

ratio_min = ratio1_min * ratio2_min
ratio_max = ratio1_max * ratio2_max

    
f, subplots = plt.subplots(3, 1, sharex = True)
subplots[0].hist(ratios1, bins = int(nsols/20))
subplots[0].plot([ratio1_min, ratio1_max], [0,0], 'o-r')
subplots[0].set_xlabel('Ratio 1')
subplots[1].hist(ratios2, bins = int(nsols/20))
subplots[1].plot([ratio2_min, ratio2_max], [0,0], 'o-r')
subplots[1].set_xlabel('Ratio 2')
subplots[2].hist(ratios, bins = 200)
subplots[2].plot([ratio_min, ratio_max], [0,0], 'o-r')
subplots[2].set_xlabel('ratio1 * ratio2')


f, subplots = plt.subplots(2, 1, sharex = True)
subplots[0].hist(cd1, bins = int(nsols/10))
subplots[0].plot([cd1_min, cd1_max], [0,0], 'o-r')
subplots[0].set_xlabel('Center distance 1')
subplots[1].hist(cd2, bins = int(nsols/10))
subplots[1].plot([cd2_min, cd2_max], [0,0], 'o-r')
subplots[1].set_xlabel('Center distance 2')
