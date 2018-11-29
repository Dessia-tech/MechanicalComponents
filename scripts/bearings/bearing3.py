#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""
import sys as sys
import mechanical_components.bearings as bearings
import mechanical_components.optimization.bearings as bearings_opt
import numpy as npy

S1 = bearings_opt.BearingAssemblyOptimizer(boundaries = [[(-0.001, 0.005, 0), (2000, -2500, 100), (0, 100, 0)]], 
                    speeds = [200], times = [1e6],
                    d_shaft_min = [0.02, 0.025], axial_pos = [0, 0.1], d_ext = [0.05, 0.07], 
                    length = [0.07, 0.04],
                    typ_linkage = [['cylindric_joint'], ['cylindric_joint']],
                    typ_mounting = None,
                    number_bearing=[[1], [2]],
                    nb_sol = [2, 1, 1])


S1.Optimize(nb_sol = 1, verbose = True)
for sol in S1.architectures:
    sol.Plot(typ='Load')
#    sol.Graph()
#    sol.list_bearing_assembly[0].list_bearing[0].FreeCADExport('extrusion2',python_path = '/Applications/FreeCAD.app/Contents/MacOS/FreeCADCmd',
#            path_lib_freecad = '/Applications/FreeCAD.app/Contents/lib', export_types=['step'])
#for num_sol, sol in enumerate(S1.solutions):
#    S1.Plot(sol)
#    S1.Export(sol, num_sol)

results = S1.results
d = results.Dict()
import json
print(json.dumps(d))
obj = bearings.BearingAssemblyResults.DictToObject(d)
obj.architectures[0].Plot(typ='Load', box=True)

#optim = obj.DefOptimizer()

print(json.dumps(sol.PlotData()))