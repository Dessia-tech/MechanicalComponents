#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""
import sys as sys
#del sys.modules['mechanical_components.optimization']
import mechanical_components.bearings as bearings
import numpy as npy
import volmdlr as vm
import copy
from mechanical_components.load import *

b0 = bearings.AngularBallBearing(d = 0.02, D = 0.04, B = 0.01, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0.1, Cr = 1e3)
b1 = bearings.RadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, Cr = 1e3)
b2 = bearings.RadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, Cr = 1e3)
b3 = bearings.AngularBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0.1, Cr = 1e3)
b4 = bearings.AngularBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0.1, Cr = 1e3)
b5 = bearings.TaperedRollerBearing(d = 0.02, D = 0.04, B = 0.015, alpha = 0.2, i = 1, 
                                       Z = 20, Dw = 0.005, Cr = 1e3)
list_bearing = [b1, b2, b2, b2]
BA = bearings.BearingCombination(list_bearing, directions = [1, 1, 1, 1], radial_load_linkage = [True]*4, internal_pre_load = 0, 
                 connection_bi = ['left', 'right'], connection_be = ['left', 'right'], behavior_link = 'both')

#print(hash(b1))
#b1_copy = b1.Copy()
#print(hash(b0))
#print(b1 == b1_copy)
#print(hash(b1), hash(b1_copy))
#BA.SolveAxialLoad()
#BA.PlotGraph()
#fa = BA.SearchBestGraph()
#BA.BearingCombinationLoad(fr=1, fa=0)
#axial_load, axial_pre_load = BA.CheckViabilityAxialPath(list_bearing, 0)
#d = BA.Dict()
#obj = bearings.BearingCombination.DictToObject(d)
#d = obj.Dict()
#print(d['bearings'])
#obj.Plot(typ=None, box=False)
#print(obj.PlotData(typ='Load'))
#BA.Plot(box = False, typ = 'Load')

#BA.PlotGraph()

#d = BA.bearings[3].Dict()
#obj = bearings.RadialBearing.DictToObject(d)
#obj.Plot()
#import json
##print(json.dumps(d))
#
#sol = bearings.BearingCombination.DictToObject(d)
#sol.Plot(typ='Load', box=False)
#
#bg = BA.bearings_solution[0].Plot(typ=None)
#export = BA.bearings[2].PlotData()
#
#sol = BA.PlotData()
##print(BA.PlotD3())
##print(json.dumps(export))
#print(json.dumps(BA.PlotData()))
#
##ax = bg.MPLPlot(style='-ob')
##li = []
##for item in bg.basis_primitives:
##    if 'Arc2D' in str(item.__class__):
##        li.extend(item.Discret())
##c=vm.Contour2D(li)
##c.MPLPlot(style='ob')