#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""
import sys as sys
#del sys.modules['mechanical_components.optimization']
import mechanical_components.bearings as bearings
import plot_data
import dessia_common as dc

import numpy as npy
import volmdlr as vm
import json
import pkg_resources
                 
b0 = bearings.TaperedRollerBearing(d = 0.02, D = 0.05, B = 0.02, i = 1,
                                        Z = 20, alpha=0.2)

plots = b0.plot_data(pos=0)
#pdg = plot_data.plot_d3(plots)

b0_bis = dc.dict_to_object(b0.to_dict())
if not b0_bis == b0:
    raise KeyError('Non equal bearing object with dict_to_object')

# b1 = bearings.RadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
#                                        Z = 20, Dw = 0.005)
# b1_bis = dc.dict_to_object(b1.to_dict())
# if not b1_bis == b1:
#     raise KeyError('Non equal bearing object with dict_to_object')
    
d = b0.plot_data(stroke_width=1)

plot_data.plot_canvas(d[0],debug_mode=False)
