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

S1 = bearings_opt.BearingAssemblyOptimizer(loads = [[[(-0.001, 0.005, 0), (2000, -2500, 100), (0, 100, 0)], 
                                                      [(0, 0.002, 0), (2000, -5000, 1000), (0, 100, 0)]]], 
                    speeds = [100], operating_times = [1e6],
                    inner_diameter = [0.02, 0.025], axial_positions = [0, 0.1], outer_diameter = [0.1, 0.1], 
                    length = [0.09, 0.08],
                    linkage_types = [['all'], ['cylindric_joint']],
                    mounting_types = [['pn', 0]],
                    number_bearings=[[1, 2, 3], [2]],
                    number_solutions = [-1, 3, 3],
                    sort_optim = {'typ': 'L10', 'min':300, 'max':1e10})
S1.Optimize()

results = S1.bearing_assembly_results

d = S1.Dict()
obj = bearings_opt.BearingAssemblyOptimizer.DictToObject(d)
d = obj.Dict()

#for num_sol, result in enumerate(results):
#    bcs = result.bearing_assembly_simulation_result.bearing_combination_simulation_results
#    for bc in bcs:
#        for bg in bc.bearing_simulation_results:
#            print(num_sol, bg.L10)
#            
for num_sol, result in enumerate(obj.bearing_assembly_results):
    print(num_sol, result.bearing_assembly_simulation_result.L10)
    result.bearing_assembly.Plot()
    
#d = results[0].Dict()
#obj = bearings.BearingAssembly.DictToObject(d)
#obj.Dict()
#obj.Plot(typ = 'Load', ind_load_case = 0)
#
#a=results[0].bearing_combinations[0].Dict()
#obj_a = bearings.BearingCombination.DictToObject(a)
#obj_a.Plot(typ = 'Load', box=False)
#results[0].bearing_combinations[0].Plot(typ = 'Load')
#
#ba = obj.bearing_assemblies[0]
#bc_result = obj.results[ba][0]['bearing_combinations'][1]
#bc = ba.bearing_combinations[1]
#
#bc.Plot(typ='Load', bearing_combination_result = bc_result)
#    sol.Graph()
#    sol.list_bearing_assembly[0].list_bearing[0].FreeCADExport('extrusion2',python_path = '/Applications/FreeCAD.app/Contents/MacOS/FreeCADCmd',
#            path_lib_freecad = '/Applications/FreeCAD.app/Contents/lib', export_types=['step'])
#for num_sol, sol in enumerate(S1.solutions):
#    S1.Plot(sol)
#    S1.Export(sol, num_sol)


#d = results.Dict()
#import json
#print(json.dumps(d))
#obj = bearings.BearingAssemblyOptimizationResults.DictToObject(d)
#obj.bearing_assemblies[0].Plot(typ='Load', box=True)
#
##optim = obj.DefOptimizer()
#
#print(json.dumps(sol.PlotData()))