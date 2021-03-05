#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
import mechanical_components.wires as wires
import mechanical_components.optimization.wires_protected as wires_opt

import volmdlr as vm

p1 = vm.Point3D(0, 0, 0)
p2 = vm.Point3D(0, 0.2, 0)
p3 = vm.Point3D(0.3, 0.23, 0)
p4 = vm.Point3D(0.45, 0.15, 0)
p5 = vm.Point3D(0.56, 0.28, 0)
p6 = vm.Point3D(0.56, 0.12, 0)
p7 = vm.Point3D(-0.10, 0.25, 0)

waypoints = [p1, p2, p3, p4, p5, p6, p7]
routes = [(p1, p2), (p7, p2), (p2, p3), (p3, p4), (p4, p5), (p4, p6)]

wires_specs = [wires.RoutingSpec(source=p1, destination=p4, diameter=0.004),
               wires.RoutingSpec(source=p7, destination=p5, diameter=0.006),
               wires.RoutingSpec(source=p2, destination=p6, diameter=0.003)
              ]
wo = wires_opt.WiringOptimizer(waypoints, routes)
wiring = wo.route(wires_specs)

wiring.spaced_wires()

# wiring.Draw()
# wiring.babylonjs()
