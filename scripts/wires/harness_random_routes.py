#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 19:17:38 2018

Network is grid randomly routed
"""

import mechanical_components.wires as wires
import mechanical_components.optimization.wires as wires_opt
import numpy as npy
import random

import volmdlr as vm


n_wpts = (20, 7, 4)# Length, width, heightvm.Line2D(waypoints[i],waypoints[i+1]).DirectionVector(unit=True).Dot(vm.Line2D(waypoints[i+1], waypoints[i+2]).DirectionVector(unit=True))==1
grid_size = (0.16, 0.15, 0.08)
min_length_paths = n_wpts[0] + 2
n_wires = 5
connection_probability = 0.7

waypoints = []
for i in range(n_wpts[0]):
    for j in range(n_wpts[1]):
        for k in range(n_wpts[2]):
            grid_point = vm.Point3D(i*grid_size[0], j*grid_size[1], k*grid_size[2])
            waypoints.append(grid_point)

routes = []            
for i in range(n_wpts[0]):
    for j in range(n_wpts[1]):
        for k in range(n_wpts[2]-1):
            if random.random() < connection_probability:
                routes.append((waypoints[i*n_wpts[1]*n_wpts[2] + j*n_wpts[2] +k],
                               waypoints[i*n_wpts[1]*n_wpts[2] + j*n_wpts[2] +k+1]))

for i in range(n_wpts[0]):
    for k in range(n_wpts[2]):
        for j in range(n_wpts[1]-1):
            if random.random() < connection_probability:
                routes.append((waypoints[i*n_wpts[1]*n_wpts[2] + j*n_wpts[2] +k],
                               waypoints[i*n_wpts[1]*n_wpts[2] + (j+1)*n_wpts[2] +k]))
                
for j in range(n_wpts[1]):
    for k in range(n_wpts[2]):
        for i in range(n_wpts[0]-1):
            if random.random() < connection_probability:
                routes.append((waypoints[i*n_wpts[1]*n_wpts[2] + j*n_wpts[2] +k],
                               waypoints[(i+1)*n_wpts[1]*n_wpts[2] + j*n_wpts[2] +k]))


wo = wires_opt.WiringOptimizer(routes)


wires_specs = []
connected_sources = []
for i in range(n_wires):
    source = random.choice(waypoints[:n_wpts[0]])
    destination = random.choice(waypoints[-n_wpts[0]:])
    
    wires_specs.append(wires.RoutingSpec(source=source,
                                         destination=destination,
                                         diameter=0.005 + 0.005*random.random()))
    

# Connecting unconnected subnetworks


# wo.plot_graph()

wiring = wo.route(wires_specs)

wiring.babylonjs()

# wiring.Draw(vm.X3D, vm.Y3D)
# wiring.Draw(vm.Y3D, vm.Z3D)
# wiring.Draw(vm.Z3D, vm.X3D)
