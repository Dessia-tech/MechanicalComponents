#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 17:13:16 2018

@author: Pierrem
"""
import mechanical_components.optimization.bearings as bearings
import numpy as npy

B1 = bearings.RadialRollerBearing(d = 0.02, D = 0.04, B = 0.015, i = 1,
                                  Z = 20, Dw = 0.005, typ = 'NF')

B1.FreeCADExport(fcstd_filepath = 'bearing1', python_path = '/Applications/FreeCAD.app/Contents/MacOS/FreeCADCmd',
            path_lib_freecad = '/Applications/FreeCAD.app/Contents/lib')


B2 = bearings.AngularBallBearing(d = 0.02, D = 0.04, B = 0.01, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0, direction = -1)
B2.FreeCADExport(fcstd_filepath = 'bearing2', python_path = '/Applications/FreeCAD.app/Contents/MacOS/FreeCADCmd',
            path_lib_freecad = '/Applications/FreeCAD.app/Contents/lib')

B3 = bearings.RadialBallBearing(d = 0.02, D = 0.04, B = 0.015, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0)
B3.FreeCADExport(fcstd_filepath = 'bearing3', python_path = '/Applications/FreeCAD.app/Contents/MacOS/FreeCADCmd',
            path_lib_freecad = '/Applications/FreeCAD.app/Contents/lib')

B4 = bearings.TaperedRollerBearing(d = 0.02, D = 0.04, B = 0.01, i = 1, 
                                       Z = 20, Dw = 0.005, alpha = 0.2, direction = -1)
B4.FreeCADExport(fcstd_filepath = 'bearing4', python_path = '/Applications/FreeCAD.app/Contents/MacOS/FreeCADCmd',
            path_lib_freecad = '/Applications/FreeCAD.app/Contents/lib')