#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 17:13:16 2018

@author: Pierrem
"""

import mechanical_components.bearing as bearing

C1=bearing.BearingCombination()
input_dict=(
          {'type':'d','min':0.01,'max':0.03},
          {'type':'B','min':0.01,'max':0.03},
          {'type':'D','min':0.03,'max':0.05},
          {'type':'grade','nom':'Gr_gn'},
          {'type':'Fr','nom':15000},
          {'type':'Fa','nom':0},
          {'type':'n','nom':100},
          {'type':'S','nom':0.97},
          {'type':'T','nom':40},
          {'type':'oil_name','nom':'iso_vg_100'},
          {'type':'nb_sol','nom':5},
          {'type':'mini','nom':'C0r'},
          {'type':'typ','nom':'NF'}
          )
C1.OptimizerBearing(input_dict)
for i,j in enumerate(C1.solution):
    C1.solution[i].FreeCADExport('GearAssembly_{}'.format(i),['fcstd','stl'])

print(C1.solution)