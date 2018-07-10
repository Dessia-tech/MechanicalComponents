#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 17:13:16 2018

@author: Pierrem
"""

import mechanical_components.bearings as bearings
import numpy as npy

C1=bearings.BearingCombination()

#Test 1
Fa=1340.1788731883905
Fa=3000
Fr=3373.492621666226
N=249.21560590286492 
L10=10
C1.OptimizerBearing(d={'min':0.04,'max':0.10},D={'min':0.04,'max':0.15},B={'min':0.01,'max':0.08},
                    L10={'min':L10,'max':1e10*L10},
                    Fr=Fr,Fa=Fa,n=N,mini=['D'],typ='NF')
for i,b in enumerate(C1.solution):
    print(b)
    v=b.VolumeModel(npy.random.random(3),npy.random.random(3))
    v.FreeCADExport('python','Bearing_{}'.format(i),'/usr/lib/freecad/lib')
    
#Test 2
Fa=1340.1788731883905
Fa=3000
Fr=3373.492621666226
N=249.21560590286492 
L10=10
C1.OptimizerBearing(d={'nom':0.06,'err':0.03},D={'nom':0.1},B={'min':0.01,'max':0.08},
                    L10={'min':L10,'max':1e10*L10},
                    Fr=Fr,Fa=Fa,n=N,mini=['D'],typ='NF')
for i,b in enumerate(C1.solution):
    print(b)
    v=b.VolumeModel(npy.random.random(3),npy.random.random(3))
    v.FreeCADExport('python','Bearing_{}'.format(i),'/usr/lib/freecad/lib')
    
#Test 3
Fa=1340.1788731883905
Fa=3000
Fr=3373.492621666226
N=249.21560590286492 
L10=10
C1.OptimizerBearing(d={'nom':0.06,'err':0.3},D={'nom':0.1,'err':0.3},B={'min':0.01,'max':0.08},
                    L10={'min':L10,'max':1e10*L10},
                    Fr=Fr,Fa=Fa,n=N,mini=['mass'],typ='NF')
for i,b in enumerate(C1.solution):
    print(b)
    v=b.VolumeModel(npy.random.random(3),npy.random.random(3))
    v.FreeCADExport('python','Bearing_{}'.format(i),'/usr/lib/freecad/lib')