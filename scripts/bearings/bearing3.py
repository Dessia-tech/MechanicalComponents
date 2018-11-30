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

S1 = bearings_opt.BearingAssemblyOptimizer(boundaries = [[(-0.001, 0.005, 0), (2000, -2500, 100), (0, 100, 0)]], 
                    speeds = [200], times = [1e6],
                    d_shaft_min = [0.02, 0.025], axial_pos = [0, 0.1], d_ext = [0.05, 0.07], 
                    length = [0.07, 0.04],
                    typ_linkage = [['cylindric_joint'], ['cylindric_joint']],
                    typ_mounting = None,
                    number_bearing=[[1], [2]],
                    nb_sol = [2, 1, 1])

S1.Optimize(number_solutions = 1, verbose = True)
results = S1.results

for ba in results.bearing_assemblies:
    bc = ba.bearing_combinations[0]
    bc_analyze = results.results[ba][0]['bearing_combinations'][0]
    bc.Plot(typ='Load', bearing_combination_result = bc_analyze)
    
d = results.Dict()
obj = bearings.BearingAssemblyOptimizationResults.DictToObject(d)
obj.Dict()

ba = obj.bearing_assemblies[0]
bc_result = obj.results[ba][0]['bearing_combinations'][1]
bc = ba.bearing_combinations[1]
bc.PlotData(typ='Load', bearing_combination_result = bc_result)
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