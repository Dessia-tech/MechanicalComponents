#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""
#import sys
#import mechanical_components.bearings as bearings
import mechanical_components.optimization.bearings as bearings_opt
#import numpy as npy

bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer(
                    loads = [[[(0.1, 0, 0), (0, 10000, 0), (0, 0, 0)]]], 
                    speeds = [300],
                    operating_times = [10000*3600],
                    inner_diameters = [0.02, 0.02],
                    axial_positions = [0, 0.4], 
                    outer_diameters = [0.07, 0.07], 
                    lengths = [0.05, 0.05],
#                    linkage_types = [['cylindric_joint'], ['cylindric_joint']],
#                    mounting_types = [['free', 'both'], ['right', 'left']],
                    number_bearings = [[1, 2, 3], [1, 2, 3]],
                    bearing_classes = [bearings_opt.RadialBallBearing, 
                                       bearings_opt.AngularBallBearing,
                                       bearings_opt.TaperedRollerBearing,
                                       bearings_opt.NUP
                                       ])

d = bearing_assembly_opt.Dict()
del bearing_assembly_opt
bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer.DictToObject(d)

bearing_assembly_opt.Optimize(2)

for num_sol, ba_simulation in enumerate(bearing_assembly_opt.bearing_assembly_simulations):
    print(num_sol, ba_simulation.bearing_assembly.mass, ba_simulation.bearing_assembly_simulation_result.L10)
    ba_simulation.bearing_assembly.Plot()
    
#d = bearing_assembly_opt.Dict()
#del bearing_assembly_opt
#bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer.DictToObject(d)
#bearing_assembly_opt.Optimize(3)