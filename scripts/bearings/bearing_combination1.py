#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""


import mechanical_components.bearings as bearings
import mechanical_components.optimization.bearings as bearings_opt
#from mechanical_components.load import *

b0 = bearings.AngularBallBearing(d = 0.02, D = 0.04, B = 0.01, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0.1, Cr = 1e3)
b1 = bearings.RadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, Cr = 1e3)
b2 = bearings.RadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, Cr = 1e3)
b3 = bearings.AngularBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0.1, Cr = 1e3)
b4 = bearings.AngularBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0.1, Cr = 1e3)
b5 = bearings.RadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, Cr = 1e3)
list_bearing = [b1, b2, b2, b2]
BA = bearings.BearingCombination(list_bearing, directions = [1, 1, 1, 1], radial_load_linkage = [True]*4, internal_pre_load = 0, 
                 connection_bi = ['left', 'right'], connection_be = ['left', 'right'], behavior_link = 'both')
print(hash(BA))
print(BA == BA)


li_bg_results = []
for bearing in BA.bearings:
    li_bg_results.append(bearings.BearingSimulationResult())
bcs = bearings.BearingCombinationSimulationResult(li_bg_results,
                    axial_loads = [1000, 2000], radial_loads = [2500, 3000],
                    speeds = [100, 200], operating_times = [1e6, 1e7])
BA.BaseLifeTime(bcs)

BCO = bearings_opt.BearingCombinationOptimizer(radial_loads = [100, 2000], 
                                           axial_loads = [0, 0], 

                                           speeds = [100, 150], 
                                           operating_times = [1e6, 1e8],
                                           inner_diameter = 0.04,
                                           outer_diameter = 0.1,
                                           length = 0.1,
                                           linkage_types = ['ball_joint', 'cylindric_joint'],
                                           mounting_types = ['free', 'right', 'both'],
                                           number_bearings = [3],
                                           bearing_classes = [bearings.RadialBallBearing, 
                                               bearings.AngularBallBearing,
                                               bearings.TaperedRollerBearing,
                                               bearings.NUP
                                               ],)

print(hash(BCO))
print(BCO == BCO)
    
BCO.Optimize(10)

for num_sol, bc_simulation in enumerate(BCO.bearing_combination_simulations):
    print(num_sol, bc_simulation.bearing_combination.mass, bc_simulation.bearing_combination_simulation_result.L10)
    bc_simulation.bearing_combination.Plot()
    print(hash(bc_simulation))
    print(bc_simulation == bc_simulation)
