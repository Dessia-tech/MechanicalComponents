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

bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer(loads = [[[(-0.001, 0.005, 0), (2000, -2500, 100), (0, 100, 0)], 
                                                      [(0, 0.002, 0), (2000, -5000, 1000), (0, 100, 0)]]], 
                    speeds = [100], operating_times = [1e6],
                    inner_diameter = [0.02, 0.025], axial_positions = [0, 0.1], outer_diameter = [0.1, 0.1], 
                    length = [0.09, 0.08],
                    linkage_types = [['all'], ['cylindric_joint']],
                    mounting_types = [['pn', 0]],
                    number_bearings=[[1, 2, 3], [2]],
                    number_solutions = [-1, 3, 3],
                    sort_optim = {'typ': 'L10', 'min':300, 'max':1e10})
bearing_assembly_opt.Optimize()

results = bearing_assembly_opt.bearing_assembly_results

d = bearing_assembly_opt.Dict()
obj = bearings_opt.BearingAssemblyOptimizer.DictToObject(d)
d = obj.Dict()

        
for num_sol, result in enumerate(obj.bearing_assembly_results):
    print(num_sol, result.bearing_assembly_simulation_result.L10)
    result.bearing_assembly.Plot()
    
