#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 19:17:38 2018

Network is grid randomly routed
"""

#import mechanical_components.wires as wires
import mechanical_components.optimization.wires as wires_opt
import numpy as npy
import random

import volmdlr as vm


n_wpts = (12, 4, 3)# Length, width, heightvm.Line2D(waypoints[i],waypoints[i+1]).DirectionVector(unit=True).Dot(vm.Line2D(waypoints[i+1], waypoints[i+2]).DirectionVector(unit=True))==1
grid_size = (0.25, 0.15, 0.08)
min_length_paths = n_wpts[0] + 2
n_wires = 5

waypoints = []
for i in range(n_wpts[0]):
    for j in range(n_wpts[1]):
        for k in range(n_wpts[2]):
            grid_point = vm.Point3D((i*grid_size[0], j*grid_size[1], k*grid_size[2]))
#            random_deviation = vm.Point3D(0.2 * (npy.dot((npy.random.random(3)-0.5), grid_size)))
            waypoints.append(grid_point)

routes = []            
for i in range(n_wpts[0]):
    for j in range(n_wpts[1]):
        for k in range(n_wpts[2]-1):
            if random.random()< 0.7:
                routes.append((waypoints[i*n_wpts[1]*n_wpts[2] + j*n_wpts[2] +k],
                               waypoints[i*n_wpts[1]*n_wpts[2] + j*n_wpts[2] +k+1]))

for i in range(n_wpts[0]):
    for k in range(n_wpts[2]):
        for j in range(n_wpts[1]-1):
            if random.random()< 0.7:
                routes.append((waypoints[i*n_wpts[1]*n_wpts[2] + j*n_wpts[2] +k],
                               waypoints[i*n_wpts[1]*n_wpts[2] + (j+1)*n_wpts[2] +k]))
                
for j in range(n_wpts[1]):
    for k in range(n_wpts[2]):
        for i in range(n_wpts[0]-1):
            if random.random()< 0.7:
                routes.append((waypoints[i*n_wpts[1]*n_wpts[2] + j*n_wpts[2] +k],
                               waypoints[(i+1)*n_wpts[1]*n_wpts[2] + j*n_wpts[2] +k]))

wires_specs = []
for i in range(n_wires):
    source = random.choice(waypoints[:n_wpts[0]])
    destination = random.choice(waypoints[-n_wpts[0]:])
    
    wires_specs.append({'source': source,
                        'destination': destination,
                        'diameter': 0.005 + 0.005*random.random()})
    

wo = wires_opt.WiringOptimizer(waypoints, routes)

wiring = wo.Optimize(wires_specs)

wiring.CADExport('harness')