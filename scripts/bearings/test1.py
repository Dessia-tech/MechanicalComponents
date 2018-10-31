#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""
import sys as sys
#del sys.modules['mechanical_components.optimization']
import mechanical_components.optimization.bearings as bearings
import numpy as npy

B1 = bearings.BearingAssemblyOptimizer(linkage = 'cylindric_joint', behavior_link = 'pn',
                                       nb_rlts = 3, d = [0.02, 0.03], D = [0.04, 0.06], 
                                       length = [0, 0.08], nb_sol = [10, 2])

for bg in B1.solutions:
    bg.Plot(box = False)
    bg.PlotGraph()