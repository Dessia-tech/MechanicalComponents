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
                    typ_linkage = [['cylindric_joint'], ['cylindric_joint']],
                    typ_mounting = [(0, 'pn')],
                    number_bearing=[[2], [1]])

#S1.Optimization(10, Lnm_min = 1e8)
#for num_sol, sol in enumerate(S1.solutions):
#    S1.Plot(sol)
#    S1.Export(sol, num_sol)