#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 12:11:20 2020

@author: launay
"""

import volmdlr as vm
from volmdlr.primitives2D import ClosedRoundedLineSegments2D
from volmdlr.primitives3D import RevolvedProfile

p1 = vm.Point2D((-0.08, 0.1755))
p2 = vm.Point2D((2, 0.1755))
p3 = vm.Point2D((2, 0.33))
p4 = vm.Point2D((-0.08, 0.33))
points1 = [p1, p2, p3, p4 ]

c1 = vm.Polygon2D(points1)
c2 = ClosedRoundedLineSegments2D(points1, {0: 0.002, 1: 0.002, 3: 0.002})

c1.MPLPlot(plot_points=True)
c2.MPLPlot(plot_points=True)

p21 = vm.Point2D((0.1, 0.))
p22 = vm.Point2D((0.1, 0.1))
p23 = vm.Point2D((0.08, 0.1))
p24 = vm.Point2D((0.08, 0.2))
p25 = vm.Point2D((0., 0.2))
p26 = vm.Point2D((0., 0.05))
p27 = vm.Point2D((0.04, 0.05))
p28 = vm.Point2D((0.04, 0.))

points2 = [p21, p22, p23, p24, p25, p26, p27, p28]
c3 = vm.Polygon2D(points2)
c4 = ClosedRoundedLineSegments2D(points2, {0: 0.002, 1: 0.002, 3: 0.002})

c3.MPLPlot(plot_points=True)
c4.MPLPlot(plot_points=True)
vector_1=vm.Vector3D((0,2,0))
profile1 = RevolvedProfile(vector_1, vm.Y3D, vm.Z3D, c1, vm.O3D, vm.Y3D)
# profile2 = RevolvedProfile(0.5*vm.X3D, vm.Y3D, vm.Z3D, c2, 0.5*vm.X3D, vm.Y3D)
# profile3 = RevolvedProfile(-0.5*vm.X3D, vm.Y3D, vm.Z3D, c1, -0.5*vm.X3D, vm.Y3D)
# profile4 = RevolvedProfile(vm.X3D, vm.Y3D, vm.Z3D, c2, vm.X3D, vm.Y3D)

model = vm.VolumeModel([profile1])
model.babylonjs()
