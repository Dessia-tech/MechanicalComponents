#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 17:13:16 2018

@author: Pierrem
"""

import mechanical_components.bearings as bearings

C1=bearings.BearingCombination()
input_dict=(
                  {'type':'d','min':0.06,'max':0.15},
                  {'type':'B','min':0.02,'max':0.06},
                  {'type':'D','min':0.07,'max':0.15},
                  {'type':'L10','min':10, 'max':2100000000},
                  {'type':'grade','nom':'Gr_gn'},
                  {'type':'Fr','nom':3000},
                  {'type':'Fa','nom':10},
                  {'type':'n','nom':1500},
                  {'type':'S','nom':0.90},
                  {'type':'T','nom':40},
                  {'type':'oil_name','nom':'iso_vg_100'},
                  {'type':'nb_sol','nom':1},
                  {'type':'mini','nom':'mass'},
                  {'type':'typ','nom':'NF'}
                  )

C1.OptimizerBearing(input_dict)
for i,b in enumerate(C1.solution):
    print(b)
    b.FreeCADExport('Bearing_{}'.format(i),['fcstd','stl'])

print(C1.solution)