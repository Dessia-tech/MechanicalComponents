#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""
import sys as sys
del sys.modules['mechanical_components.optimization']
import mechanical_components.optimization.bearings as bearings
import numpy as npy

S1 = bearings.CompositeBearingAssemblyOptimizer(list_pos_unknown = [[-0.001,0.005,0]], 
                    list_load = [[2000, -2500, 0]], list_torque = [[0, 100, 0]],
                    list_speed = [200], list_time = [1e6],
                    d_shaft_min = [0.02, 0.025], axial_pos = [0, 0.1], d_ext = [0.05, 0.07], 
                    length = [0.07, 0.04],
                    typ_linkage = [['all'], ['all']],
                    typ_mounting = None,
                    number_bearing=[[1, 2], [1, 2]],
                    nb_sol = [None, None, 3])

print(len(S1.solutions))
#for sol in S1.solutions:
#    print(sol.mass)
    
S1.Optimize(nb_sol = 20, verbose = True)
for sol in S1.solutions:
    sol.Plot()
    print(sol.pos_x)
#    sol.list_bearing_assembly[0].list_bearing[0].FreeCADExport('extrusion2',python_path = '/Applications/FreeCAD.app/Contents/MacOS/FreeCADCmd',
#            path_lib_freecad = '/Applications/FreeCAD.app/Contents/lib', export_types=['step'])
#for num_sol, sol in enumerate(S1.solutions):
#    S1.Plot(sol)
#    S1.Export(sol, num_sol)