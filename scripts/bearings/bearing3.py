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

S1 = bearings.CompositeBearingAssemblyOptimizer(list_pos_unknown = [[0.09,0.004,0]], 
                    list_load = [[200, -2500, 0]], list_torque = [[0, 100, 0]],
                    list_speed = [200], list_time = [1e6],
                    d_shaft_min = [0.02, 0.025], axial_pos = [0, 0.1], d_ext = [0.05, 0.07], 
                    length = [0.04, 0.04],
                    typ_linkage = [['ball_joint'], ['cylindric_joint']],
                    typ_mounting = [(0, 'pn')],
                    number_bearing=[[1, 2], [1, 2]],
                    nb_sol = [20, 10, 1])

for sol in S1.solutions:
    print(sol.mass)
    
S1.Optimize(index_sol = [1,2,3,4], nb_sol = 5, verbose = True)
for sol in S1.solutions:
    sol.Plot()
#for num_sol, sol in enumerate(S1.solutions):
#    S1.Plot(sol)
#    S1.Export(sol, num_sol)