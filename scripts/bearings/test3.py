#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""
import sys as sys
#del sys.modules['mechanical_components.optimization']
import mechanical_components.bearings as bearings
from volmdlr import plot_data

import numpy as npy
import volmdlr as vm
import json
import pkg_resources
                 
b0 = bearings.TaperedRollerBearing(d = 0.02, D = 0.05, B = 0.02, i = 1, 
                                       Z = 20, alpha=0.2)
print(b0.Dw)
#d = b0.to_dict()
#jd = json.dumps(d)
#b0_bis = bearings.RadialBallBearing.dict_to_object(d)
#b0_bis = DessiaObject.dict_to_object(d)

plots = b0.plot_data(pos=0, direction=-1, quote=False, constructor=False)
pdg = plot_data.plot_d3(plots)

#with pkg_resources.resource_stream(pkg_resources.Requirement('mechanical_components'),
#                           'mechanical_components/catalogs/schaeffler_v2.json') as schaeffler_json:
#    schaeffler_catalog = bearings.BearingCatalog.LoadFromFile(schaeffler_json)

#schaeffler_catalog.SaveToFile('essai')

#%%
b0_bis = bearings.RadialBearing.DictToObject(d)
b0_bis.Plot(typ=None)
b0_bis.Graph()
b0_bis.Plot(typ=None)

b1 = bearings.RadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0)
d = b1.Dict()
jd = json.dumps(d)
b1_bis = bearings.RadialBallBearing.DictToObject(d)
b1_bis = bearings.RadialBearing.DictToObject(d)
b1_bis.Plot()
export = b1_bis.PlotData(quote = True)
print(json.dumps(export))

b2 = bearings.RadialRollerBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005)
d = b2.Dict()
jd = json.dumps(d)
b2_bis = bearings.AngularBallBearing.DictToObject(d)
b2_bis = bearings.RadialBearing.DictToObject(d)
b2_bis.Plot()
b2_bis.PlotData()

b3 = bearings.TaperedRollerBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.004, alpha = 0.2)
d = b3.Dict()
jd = json.dumps(d)
b3_bis = bearings.AngularBallBearing.DictToObject(d)
b3_bis = bearings.RadialBearing.DictToObject(d)
b3_bis.Plot()
