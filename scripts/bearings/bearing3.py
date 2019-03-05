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
                    loads = [[[(0.1, 0, 0), (2000, 10000, 0), (0, 0, 0)], 
                              [(0.3, 0, 0), (1000, 3000, 0), (0, 0, 0)]]], 
                    speeds = [1000],
                    operating_times = [10000*3600],
                    inner_diameters = [0.02, 0.02],
                    axial_positions = [0, 0.2], 
                    outer_diameters = [0.1, 0.1], 
                    lengths = [0.08, 0.08],
                    linkage_types = [['cylindric_joint'], ['cylindric_joint']],
                    mounting_types = [['free', 'both'], ['right', 'left']],
                    number_bearings = [[1, 2], [1, 2]],
                    bearing_classes = [bearings_opt.RadialBallBearing, 
                                       bearings_opt.AngularBallBearing,
                                       bearings_opt.TaperedRollerBearing,
                                       bearings_opt.NUP, bearings_opt.N, bearings_opt.NU,
                                       bearings_opt.NF
                                       ])

d = bearing_assembly_opt.Dict()
del bearing_assembly_opt
bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer.DictToObject(d)

bearing_assembly_opt.Optimize(10)

for num_sol, ba_simulation in enumerate(bearing_assembly_opt.bearing_assembly_simulations):
    print(num_sol, ba_simulation.bearing_assembly.mass, ba_simulation.bearing_assembly_simulation_result.L10)
    ba_simulation.bearing_assembly.Plot()
    
#d = bearing_assembly_opt.Dict()
#del bearing_assembly_opt
#bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer.DictToObject(d)
#bearing_assembly_opt.Optimize(3)