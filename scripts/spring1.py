#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 18:01:54 2018

@author: jezequel
"""
import math
import mechanical_components.springs as springs
import numpy as npy
import matplotlib.pyplot as plt
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D

n_spires = npy.linspace(0.5, 9.5, 10)
results = []
for n in n_spires:
    spring = springs.Spring()
    so = springs.SpringOptimizer(spring, {'displacement' : 0.003,
                                          'D' : (0.001, 0.3),
                                          'd' : (0.0001, 0.1),
                                          'F0' : 100,
                                          'F1' : 1000,
                                          'n' : n})
    
    res = so.Optimize()
    print(res)
    if res['success']:
        results.append(so.spring)
    print('---------------------------------------------\n',
          'Resultat', n, '\n\n',
          res)
    print('D : ', so.spring.D,
          '\nd : ', so.spring.d,
          '\nn : ', so.spring.n,
          '\nh : ', so.spring.displacement,
          '\nTargetStiffness : ', so.spring.TargetStiffness(),
          '\nStiffness : ', so.spring.Stiffness(),
          '\nEngagedResultingForce : ', so.spring.EngagedResultingForce(),
          '\nLoseResultingForce : ', so.spring.LoseResultingForce(),
          '\nLengths : ', so.spring.Lengths(),
          '\n---------------------------------------------',
          '\n\n\n')
    
a = [i.d**4/(i.D**3) for i in results]
b = [i.D/i.d for i in results]
n_spires = [i.n for i in results]
print(a)

plt.plot(n_spires, a, 'ro')