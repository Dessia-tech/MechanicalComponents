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
import volmdlr as vm

b0 = bearings.ConceptAngularBallBearing(d = 0.02, D = 0.04, B = 0.01, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0, direction = -1)
b1 = bearings.ConceptRadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0)
b2 = bearings.ConceptRadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0)
b3 = bearings.ConceptAngularBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0, direction = 1)
b4 = bearings.ConceptAngularBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0, direction = -1)
b5 = bearings.ConceptRadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0)
list_bearing = [b1, b2, b3, b4, b5]
BA = bearings.BearingCombination(list_bearing, radial_load_linkage = [True]*5, internal_pre_load = 0, 
                 connection_bi = ['p'], connection_be = ['n', 'p'], behavior_link = 'pn')

#BA.PlotGraph()
fa = BA.BearingCombinationLoad(fa = 0, fr = 200)
#BA.Plot(box = False, typ = 'Load')

#BA.PlotGraph()

d = BA.Dict()
import json
#print(json.dumps(d))

sol = bearings.BearingCombination.Dict2Obj(d)
sol.Plot(typ='Load', box=False)

bg = BA.list_bearing[2].Plot()
export = BA.list_bearing[2].PlotData()

sol = BA.PlotData()
#print(BA.PlotD3())
#print(json.dumps(export))
print(json.dumps(BA.list_bearing[0].PlotData()))

#ax = bg.MPLPlot(style='-ob')
#li = []
#for item in bg.basis_primitives:
#    if 'Arc2D' in str(item.__class__):
#        li.extend(item.Discret())
#c=vm.Contour2D(li)
#c.MPLPlot(style='ob')