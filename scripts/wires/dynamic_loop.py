#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 11:53:29 2021

@author: steven
"""


import mechanical_components.wires as wires

import volmdlr as vm

p1 = vm.Point3D(0, 0, 0)
p2 = vm.Point3D(0.05, 0.03, 0.12)
tan1 = vm.Vector3D(0, 1, 0)
tan2 = vm.Vector3D(0, -1, 0)

length = 1.3 * p1.point_distance(p2)

wire = wires.JunctionWire(p1, tan1, p2, tan2, length, diameter=5e-3)
ax = wire.path.plot()
# wire.babylonjs()
# wire.Draw()

wire2 = wires.JunctionWire.curvature_radius(point1=p1, tangeancy1=tan1, 
                                            point2=p2, tangeancy2=tan2, 
                                            targeted_curv=40e-3, 
                                            length_min=length*0.5, length_max=length,
                                            diameter=5e-3)
wire2.path.plot(ax=ax, color='r')

print(wire2.path.minimum_radius(), wire.path.minimum_radius() ,40e-3)