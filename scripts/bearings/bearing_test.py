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
import dessia_common as dc

import numpy as npy
import volmdlr as vm
import json
import pkg_resources
                 
b0 = bearings.TaperedRollerBearing(d = 0.02, D = 0.05, B = 0.02, i = 1, 
                                       Z = 20, alpha=0.2)

plots = b0.plot_data(pos=0, direction=-1, quote=False, constructor=False)
#pdg = plot_data.plot_d3(plots)

b0_bis = dc.dict_to_object(b0.to_dict())
if not b0_bis == b0:
    raise KeyError('Non equal bearing object with dict_to_object')

b1 = bearings.RadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005)
b1_bis = dc.dict_to_object(b1.to_dict())
if not b1_bis == b1:
    raise KeyError('Non equal bearing object with dict_to_object')
    
d = b1.plot_data()
#plot_data.plot_d3(d)

#plots = b1.plot_data(pos=0, direction=-1, quote=False, constructor=False)
#pdg = plot_data.plot_d3(plots)
#b0_bis = dc.dict_to_object(b1.to_dict())

# d = b1.Dict()
# jd = json.dumps(d)
# b1_bis = bearings.RadialBallBearing.dict_to_object(d)
# b1_bis = bearings.RadialBearing.dict_to_object(d)
# b1_bis.Plot()
# export = b1_bis.PlotData(quote = True)
# print(json.dumps(export))

# b2 = bearings.RadialRollerBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
#                                        Z = 20, Dw = 0.005)
# d = b2.Dict()
# jd = json.dumps(d)
# b2_bis = bearings.AngularBallBearing.dict_to_object(d)
# b2_bis = bearings.RadialBearing.dict_to_object(d)
# b2_bis.Plot()
# b2_bis.PlotData()

# b3 = bearings.TaperedRollerBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
#                                        Z = 20, Dw = 0.004, alpha = 0.2)
# d = b3.Dict()
# jd = json.dumps(d)
# b3_bis = bearings.AngularBallBearing.dict_to_object(d)
# b3_bis = bearings.RadialBearing.dict_to_object(d)
# b3_bis.Plot()
