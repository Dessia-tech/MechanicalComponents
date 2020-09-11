#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""


import mechanical_components.bearings as bearings
import mechanical_components.optimization.bearings as bearings_opt
#from mechanical_components.load import *


BCO = bearings_opt.BearingCombinationOptimizer(radial_loads = [643.2158843252842], 
                                           axial_loads = [6.643929650868828e-14], 

                                           speeds = [1924.046566484553], 
                                           operating_times = [1274657003.7725844],
                                           inner_diameter = 0.03,
                                           outer_diameter = 0.15,
                                           length = 0.05,
                                           linkage_types = ['ball_joint', 'cylindric_joint'],
                                           mounting_types = ['both', 'right', 'left', 'free'],
                                           number_bearings = [1, 2],
                                           bearing_classes = [bearings.RadialBallBearing, 
                                                              bearings.AngularBallBearing,
                                                              bearings.TaperedRollerBearing,
                                                              bearings.NUP
                                                              ],
                                           )

    
BCO.Optimize(10)

for num_sol, bc_simulation in enumerate(BCO.bearing_combination_simulations):
    print(num_sol, bc_simulation.bearing_combination.mass, bc_simulation.bearing_combination_simulation_result.L10)
    bc_simulation.bearing_combination.Plot()
