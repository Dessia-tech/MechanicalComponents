#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 15:49:08 2018

"""

import mechanical_components.wires as wires
import mechanical_components.optimization.wires_protected as wires_opt

import volmdlr as vm

p1 = vm.Point3D(0, 0, 0)
p2 = vm.Point3D(0.5, 0.2, 0)
p3 = vm.Point3D(0, 0.23, 0)
p4 = vm.Point3D(0.3, 0.36, 0)
p5 = vm.Point3D(0, 0.76, -0.1)
p6 = vm.Point3D(-.1, 0.23, 0)
p7 = vm.Point3D(0, 0, 0.32)

waypoints = [p1, p2, p3, p4, p5, p6, p7]
routes = [(p1, p2), (p2, p4), (p2, p4), (p4, p6), (p3, p6), (p6, p7)]

wires_specs = [wires.RoutingSpec(source=p1, destination=p7, diameter=0.008, name='wire1'),
               wires.RoutingSpec(source=p3, destination=p7, diameter=0.006, name='wire2')]

wo = wires_opt.WiringOptimizer(waypoints, routes)
wiring = wo.route(wires_specs)

wiring.Draw()
wiring.babylonjs()

