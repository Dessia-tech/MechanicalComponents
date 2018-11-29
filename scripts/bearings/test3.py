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
import json

b0 = bearings.AngularBallBearing(d = 0.02, D = 0.04, B = 0.01, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0, direction = -1)
d = b0.Dict()
jd = json.dumps(d)
b0_bis = bearings.AngularBallBearing.DictToObject(d)
b0_bis = bearings.RadialBearing.DictToObject(d)
b0_bis.Plot()

b1 = bearings.RadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0)
d = b1.Dict()
jd = json.dumps(d)
b1_bis = bearings.AngularBallBearing.DictToObject(d)
b1_bis = bearings.RadialBearing.DictToObject(d)
b1_bis.Plot()

b2 = bearings.RadialRollerBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005)
d = b2.Dict()
jd = json.dumps(d)
b2_bis = bearings.AngularBallBearing.DictToObject(d)
b2_bis = bearings.RadialBearing.DictToObject(d)
b2_bis.Plot()

b3 = bearings.TaperedRollerBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.004, alpha = 0.2)
d = b3.Dict()
jd = json.dumps(d)
b3_bis = bearings.AngularBallBearing.DictToObject(d)
b3_bis = bearings.RadialBearing.DictToObject(d)
b3_bis.Plot()