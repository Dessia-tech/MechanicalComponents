#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 17:13:16 2018

@author: Pierrem
"""

import mechanical_components.bearing as bearing

C1=bearing.BearingCombination()
input_dict=(
          {'type':'d','min':0.01,'max':0.02},
          {'type':'B','min':0.01,'max':0.02},
          {'type':'D','min':0,'max':0.04},
          {'type':'L10','min':150,'max':200},
          {'type':'C0r','min':0,'max':15000},
          {'type':'grade','nom':'Gr_gn'},
          {'type':'Fr','nom':1500},
          {'type':'Fa','nom':0},
          {'type':'n','nom':100},
          {'type':'S','nom':0.99},
          {'type':'T','nom':40},
          {'type':'oil_name','nom':'iso_vg_1000'},
          {'type':'nb_sol','nom':5},
          {'type':'mini','nom':'C0r'}
          )
C1.OptimizerBearing(input_dict)

print(C1.solution)