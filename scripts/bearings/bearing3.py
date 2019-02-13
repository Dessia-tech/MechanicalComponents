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

bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer(loads = [[[(-0.001, 0.005, 0), (2000, -2500, 100), (0, 100, 0)], 
                                                      [(0, 0.002, 0), (2000, -5000, 1000), (0, 100, 0)]]], 
                    speeds = [1000], operating_times = [1e6],
                    inner_diameter = [0.02, 0.025], axial_positions = [0, 0.1], 
                    outer_diameter = [0.1, 0.1], 
                    length = [0.09, 0.08],
                    linkage_types = [['both'], ['cylindric_joint']],
                    mounting_types = [['left', 'right'], ['right', 'right'], ['free', 'both']] ,
                    number_bearings = [[1, 2], [1, 2]],
                    number_solutions = [-1, 5, 3])

d = bearing_assembly_opt.Dict()
del bearing_assembly_opt
bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer.DictToObject(d)
#
bearing_assembly_opt.Optimize()
#
#results = bearing_assembly_opt.bearing_assembly_simulation_results
#        
for num_sol, result in enumerate(bearing_assembly_opt.bearing_assembly_simulation_results):
    print(num_sol, result.bearing_assembly.mass, result.bearing_assembly_simulation_result.L10)
    result.bearing_assembly.Plot()