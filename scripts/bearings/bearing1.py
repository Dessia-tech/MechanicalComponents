#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 17:13:16 2018

@author: Pierrem
"""
import sys
#del sys.modules['mechanical_components.optimization']
import mechanical_components.optimization.bearings as bearings
import numpy as npy

C1=bearings.RollerBearingContinuousOptimizer()

#Test 1
Fa = [2000]
Fr = [8230.476546051736]
N = [30.]
t = [250000/30*3600]
T = [70]
L10 = 1500

C1.Optimization(d = {'min':0.02,'max':0.1},
                    D = {'min':0.04,'max':0.15},
                    B = {'min':0.01,'max':0.1},
                    Lnm = {'min':L10,'max':npy.inf},
#                    L10={'min':L10,'max':npy.inf},
                    Fr = Fr, Fa = Fa, N = N, t = t, T = T, typ = 'NF', 
                    nb_sol = 10,
                    verbose = True)


for i,b in enumerate(C1.bearing_assemblies):
#    v=b.VolumeModel()
    b.FreeCADExport(fcstd_filepath = 'Bearing_{}'.format(i), python_path = '/Applications/FreeCAD.app/Contents/MacOS/FreeCADCmd',
            path_lib_freecad = '/Applications/FreeCAD.app/Contents/lib')
#    print(b.CheckFNRRules(Fr, Fa, N))
    
##Test 2
#Fa=1340.1788731883905
#Fa=3000
#Fr=3373.492621666226
#N=249.21560590286492 
#L10=10
#C1.OptimizerBearing(d={'nom':0.06,'err':0.03},D={'nom':0.1},B={'min':0.01,'max':0.08},
#                    L10={'min':L10,'max':1e10*L10},
#                    Fr=Fr,Fa=Fa,n=N,mini=['D'],typ='NF')
#for i,b in enumerate(C1.solution):
#    print(b)
#    v=b.VolumeModel(npy.random.random(3),npy.random.random(3))
#    v.FreeCADExport('python','Bearing_{}'.format(i),'/usr/lib/freecad/lib')
#    
##Test 3
#Fa=1340.1788731883905
#Fa=3000
#Fr=3373.492621666226
#N=249.21560590286492 
#L10=10
#C1.OptimizerBearing(d={'nom':0.06,'err':0.3},D={'nom':0.1,'err':0.3},B={'min':0.01,'max':0.08},
#                    L10={'min':L10,'max':1e10*L10},
#                    Fr=Fr,Fa=Fa,n=N,mini=['mass'],typ='NF')
#for i,b in enumerate(C1.solution):
#    print(b)
#    v=b.VolumeModel(npy.random.random(3),npy.random.random(3))
#    v.FreeCADExport('python','Bearing_{}'.format(i),'/usr/lib/freecad/lib')