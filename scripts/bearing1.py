#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 17:13:16 2018

@author: Pierrem
"""

import mechanical_components.bearings as bearings
import numpy as npy

C1=bearings.BearingCombination()
#input_dict=(
#                  {'type':'d','min':0.06,'max':0.15},
#                  {'type':'B','min':0.02,'max':0.06},
#                  {'type':'D','min':0.07,'max':0.15},
#                  {'type':'L10','min':10, 'max':2100000000},
#                  {'type':'grade','nom':'Gr_gn'},
#                  {'type':'Fr','nom':3000},
#                  {'type':'Fa','nom':10},
#                  {'type':'n','nom':1500},
#                  {'type':'S','nom':0.90},
#                  {'type':'T','nom':40},
#                  {'type':'oil_name','nom':'iso_vg_100'},
#                  {'type':'nb_sol','nom':1},
#                  {'type':'mini','nom':'mass'},
#                  {'type':'typ','nom':'NF'}
#                  )

Fa=1340.1788731883905
Fr=3373.492621666226
N=249.21560590286492 
L10=356.97506017574415

input_dict=(
                  {'type':'d','min':0.04,'max':0.10},
                  {'type':'B','min':0.01,'max':0.08},
                  {'type':'D','min':0.04,'max':0.15},
                  {'type':'L10','min':L10, 'max':1.20*L10},
                  {'type':'grade','nom':'Gr_gn'},
                  {'type':'Fr','nom':Fr},
                  {'type':'Fa','nom':Fa},
                  {'type':'n','nom':N},
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

    v=b.VolumeModel(npy.random.random(3),npy.random.random(3))
    v.FreeCADExport('python','Bearing_{}'.format(i),'/usr/lib/freecad/lib')
#    b.FreeCADExport('Bearing_{}'.format(i),['fcstd','stl'])

print(C1.solution)