#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""
import sys as sys
#del sys.modules['mechanical_components.optimization']
import mechanical_components.bearings as bearings
import mechanical_components.optimization.bearings as bearings_opt
import numpy as npy
import volmdlr as vm
import copy
from mechanical_components.load import *

cbc = bearings.ConceptualBearingCombination([bearings.AngularBallBearing, bearings.AngularBallBearing, 
                                             bearings.NUP],
                                            [-1, 1, -1], 'left')
print(cbc.CheckKinematic())

bco = bearings_opt.BearingCombinationOptimizer(linkage = 'ball_joint', mounting = 'left', d = [0.02, 0.025],
                                     D = [0.02, 0.025], length = [0.02, 0.025], bearing_classes=bearings.bearing_classes)
sol = bco.ConceptualBearingCombinations(2)

bc = bearings.BearingCatalogue()

b0 = bearings.AngularBallBearing(d = 0.02, D = 0.04, B = 0.01, i = 1, 
                                 Z = 20, Dw = 0.005, alpha = 0, Cr = 1e3)
d = b0.Dict()
obj = bearings.AngularBallBearing.DictToObject(d)
d = obj.Dict()
b0.Plot(direction = 1)
#
b1 = bearings.RadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0, Cr = 1e3)
d = b1.Dict()
obj = bearings.RadialBallBearing.DictToObject(d)
d = obj.Dict()
b1.Plot()

b2 = bearings.TaperedRollerBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0.1, Cr = 1e3)
d = b2.Dict()
obj = bearings.TaperedRollerBearing.DictToObject(d)
d = obj.Dict()
b2.Plot(direction = 1)
#
#b2 = bearings.RadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
#                                Z = 20, Dw = 0.005, alpha = 0, Cr = 1e3)
#
#b3 = bearings.AngularBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
#                                       Z = 20, Dw = 0.005, alpha = 0.15, Cr = 1e3)
#b4 = bearings.AngularBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
#                                       Z = 20, Dw = 0.005, alpha = 0.15, Cr = 1e3)
#b5 = bearings.RadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
#                                       Z = 20, Dw = 0.005, alpha = 0, Cr = 1e3)
#list_bearing = [b1, b2, b3, b4]
#BC1 = bearings.BearingCombination(list_bearing, radial_load_linkage = [True]*4, internal_pre_load = 0, 
#                 connection_bi = ['n'], connection_be = ['p'], behavior_link = 'pn')
#d = BC1.Dict()
#obj = bearings.BearingCombination.DictToObject(d)
#d = obj.Dict()
#
#list_bearing = [b1, b2]
#BC2 = bearings.BearingCombination(list_bearing, radial_load_linkage = [True]*2, internal_pre_load = 0, 
#                 connection_bi = ['n'], connection_be = ['p'], behavior_link = 'pn')
#
#BA = bearings.BearingAssembly([BC1, BC2], axial_positions = [0, 0.5])
#d = BA.Dict()
#obj = bearings.BearingAssembly.DictToObject(d)
#d = obj.Dict()
#                
#loads = [[[(-0.001, 0.005, 0), (2000, -2500, 100), (0, 100, 0)], 
#          [(0, 0.002, 0), (2000, -50, 1000), (0, 100, 0)]]]
#speeds = [100]
#operating_times = [1e6]
#
#bc_results = []
#for bearing_combination in BA.bearing_combinations:
#    li_bg_results = []
#    for bearing in bearing_combination.bearings:
#        li_bg_results.append(bearings.BearingSimulationResult())
#    bc_results.append(bearings.BearingCombinationSimulationResult(li_bg_results))
#BASR = bearings.BearingAssemblySimulationResult(bc_results, loads, speeds, operating_times)
#            
#BAS = bearings.BearingAssemblySimulation(BA, BASR)
#d = BAS.Dict()
#obj = bearings.BearingAssemblySimulation.DictToObject(d)
#d = obj.Dict()
#
#BA.Update([0, 0.1], [0.02, 0.025], [0, 0.1], 
#          [0.05, 0.07], [0.07, 0.04])
#BA.ShaftLoad([0, 0.1], BAS.bearing_assembly_simulation_result)
#BAS.bearing_assembly_simulation_result.axial_load_model.Plot(intensity_factor=1e-5)
#
##import json
##d = BA.Dict()
##print(json.dumps(d))
##obj = bearings.BearingAssembly.DictToObject(d)
##obj.Plot()
##obj = bearings.BearingAssembly.DictToObject(d)
##d = obj.Dict()
#
##print(d['bearings'])
##obj.Plot(typ=None, box=False)
##print(obj.PlotData(typ='Load'))
##BA.Plot(box = False, typ = 'Load')
#
##BA.PlotGraph()
#
##d = BA.bearings[3].Dict()
##obj = bearings.RadialBearing.DictToObject(d)
##obj.Plot()
##import json
###print(json.dumps(d))
##
##sol = bearings.BearingCombination.DictToObject(d)
##sol.Plot(typ='Load', box=False)
##
##bg = BA.bearings_solution[0].Plot(typ=None)
##export = BA.bearings[2].PlotData()
##
##sol = BA.PlotData()
###print(BA.PlotD3())
###print(json.dumps(export))
##print(json.dumps(BA.PlotData()))
##
###ax = bg.MPLPlot(style='-ob')
###li = []
###for item in bg.basis_primitives:
###    if 'Arc2D' in str(item.__class__):
###        li.extend(item.Discret())
###c=vm.Contour2D(li)
###c.MPLPlot(style='ob')