#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""
#import sys
import mechanical_components.bearings as bearings
import mechanical_components.optimization.bearings as bearings_opt

#import numpy as npy

#bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer(
#                    loads = [[[(0.1, 0, 0), (2000, 10000, 0), (0, 0, 0)], 
#                              [(0.3, 0, 0), (1000, 3000, 0), (0, 0, 0)]]], 
#                    speeds = [1000],
#                    operating_times = [10000*3600],
#                    inner_diameters = [0.02, 0.02],
#                    axial_positions = [0, 0.2], 
#                    outer_diameters = [0.1, 0.1], 
#                    lengths = [0.08, 0.08],
#                    linkage_types = [['cylindric_joint'], ['cylindric_joint']],
#                    mounting_types = [['free', 'both'], ['right', 'left']],
#                    number_bearings = [[2], [2]],
#                    bearing_classes = [bearings.RadialBallBearing, 
#                                       bearings.AngularBallBearing,
#                                       bearings.TaperedRollerBearing,
#                                       bearings.NUP, bearings.N, bearings.NU,
##                                       bearings_opt.NF
#                                       ])


bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer(
                    loads = [[[[0.2, 0, 0], [1000, 2000, 0], [0, 0, 0]]]], 
                    speeds = [157],
                    operating_times = [3600000],
                    inner_diameters = [0.03, 0.03],
                    axial_positions = [0, 0.25], 
                    outer_diameters = [0.08, 0.09], 
                    lengths = [0.04, 0.1],
                    linkage_types = [['ball_joint'], ['cylindric_joint']],
                    mounting_types = [['right', 'left']],
                    number_bearings = [[1, 2], [1, 2]],

                    bearing_classes = [bearings.RadialBallBearing, 
                                       bearings.AngularBallBearing,
                                       bearings.TaperedRollerBearing,
                                       bearings.NUP, bearings.N, bearings.NU,
#                                       bearings_opt.NF
                                       ])


#'axial_positions': ,
# 'bearing_assembly_simulations': [],
# 'inner_diameters': ,
# 'lengths': ,
# 'linkage_types': ,
# 'loads': ,
# 'mounting_types': ,
# 'number_bearings': ,
# 'operating_times': ,
# 'outer_diameters': ,
# 'speeds': 

print(hash(bearing_assembly_opt))
print(bearing_assembly_opt == bearing_assembly_opt)

#d = bearing_assembly_opt.Dict()
#del bearing_assembly_opt
#bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer.DictToObject(d)

bearing_assembly_opt.Optimize(12)

for num_sol, ba_simulation in enumerate(bearing_assembly_opt.bearing_assembly_simulations):
    print(num_sol, ba_simulation.bearing_assembly.mass, ba_simulation.bearing_assembly_simulation_result.L10)
    ba_simulation.bearing_assembly.Plot()    
    print(hash(ba_simulation))
    print(ba_simulation == ba_simulation)
    
#d = bearing_assembly_opt.Dict()
#del bearing_assembly_opt
#bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer.DictToObject(d)
#bearing_assembly_opt.Optimize(3)