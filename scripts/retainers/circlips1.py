#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:36:30 2019

@author: ringhausen
"""

import mechanical_components.retainers as retainers
import matplotlib.pyplot as plt
import math

#circlip = retainers.Circlips()
#circlips = retainers.ExternalCirclips(0.050,0.0069,0.0051,0.002,0.0025)
#print(circlips.ExternalContour())
#print(circlips.InternalContour())
#print(circlips.Plot())
#circlips.CADExport('circlipsss')


#circlips = retainers.InternalCirclips(0.05,0.0065,0.0046,0.002,0.0025)
#print(circlips.Plot())



#x = []
#y = []
#x = [2*math.pi*i/100 for i in range(101)]
#for i in range(101):      
#    y.append(math.atan(x[i]))
#
#plt.figure(figsize = (20,10))
#plt.plot(x, y)



#for circlip in retainers.external_circlips_catalog.circlips[0:20]:
#    print('d = ', circlip.d)
#    print(circlip.Plot())

D, F, w_max = 0.05, 10000, 0.003
shaft = [D, F, w_max, retainers.groove_materials[0]]
a = retainers.external_circlips_catalog.ExternalShaftToCirclipGroove(D, F, w_max)
print('Circlip, groove : ', a)
if a[0][1] != None:
    print('Groove', a[0][1].D, a[0][1].x, a[0][1].w)
    
    
print("\n")
    
D, F, w_max = 0.05, 10000, 0.003
shaft = [D, F, w_max, retainers.groove_materials[0]]
a = retainers.internal_circlips_catalog.InternalShaftToCirclipGroove(D, F, w_max)
print('Circlip, groove : ', a)
if a[0][1] != None:
    print('Groove', a[0][1].D, a[0][1].x, a[0][1].w)