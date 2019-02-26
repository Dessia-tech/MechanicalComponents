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
                    loads = [[[(-0.001, 0.005, 0), (5000, -250, 100), (0, 100, 0)], 
                               [(0, 0.002, 0), (2000, -5000, 1000), (0, 100, 0)]],
                              [[(-0.001, 0.005, 0), (3000, -2500, 100), (0, 100, 0)], 
                               [(0, 0.002, 0), (0, -5000, 1000), (0, 100, 0)]]], 
                    speeds = [1000],
                    operating_times = [1e7],
                    inner_diameters = [0.02, 0.025],
                    axial_positions = [0, 0.2], 
                    outer_diameters = [0.1, 0.1], 
                    lengths = [0.05, 0.1],
                    linkage_types = [['ball_joint'], ['cylindric_joint']],
                    mounting_types = [['free', 'both'], ['right', 'left']],
                    number_bearings = [[1], [1, 2]],
                    bearing_classes = bearings_opt.bearing_classes)

d = bearing_assembly_opt.Dict()
del bearing_assembly_opt
bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer.DictToObject(d)

bearing_assembly_opt.Optimize(10)

for num_sol, ba_simulation in enumerate(bearing_assembly_opt.bearing_assembly_simulations):
    print(num_sol, ba_simulation.bearing_assembly.mass, ba_simulation.bearing_assembly_simulation_result.L10)
    ba_simulation.bearing_assembly.Plot()
    
d = bearing_assembly_opt.Dict()
del bearing_assembly_opt
bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer.DictToObject(d)
#bearing_assembly_opt.Optimize(3)