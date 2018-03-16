import math
import mechanical_components.clutches as clutches
import numpy as np
import matplotlib.pyplot as plt
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 17:29:05 2018

@author: jezequel
"""

verin = clutches.HydraulicCylinder()
emb = clutches.Clutch(verin, separator_tooth_type = 'outer')

regime = np.linspace(0, 2500*math.pi/30, 100)
test = emb.DragTorque(regime, 0)
#print(test)

"""
plt.plot(regime*30/math.pi, test)
plt.show()
"""

co=clutches.ClutchOptimizer(emb, {'plate_inner_radius' : 0.060, 'plate_height' : (0.007, 0.020)})
res=co.Optimize()

## Export Freecad vérin
primitives = []
primitives.extend(co.clutch.hydraulic_cylinder.chamber_volume)
primitives.extend(co.clutch.hydraulic_cylinder.piston_volume)
#model_verin = vm.VolumeModel(primitives)
#resp=model_verin.FreeCADExport('python','cylinder','/usr/lib/freecad/lib/',['stl','fcstd'])

## Export Freecad vérin
#primitives = []
primitives.extend(co.clutch.hydraulic_cylinder.spring_volume)
#model_ressort = vm.VolumeModel(primitives)
#resp=model_ressort.FreeCADExport('python','spring','/usr/lib/freecad/lib/',['stl','fcstd'])
#print(resp)

# Export Freecad embrayage
#primitives = []
primitives.extend(co.clutch.separator_plate_volume)
primitives.extend(co.clutch.friction_plate_volume)
model_emb = vm.VolumeModel(primitives)
resp=model_emb.FreeCADExport('python','plate','/usr/lib/freecad/lib/',['stl','fcstd'])
print(resp)

#p0 = vm.Point2D((0, 0))
#p1 = vm.Point2D((0.09, 0))
#p2 = vm.Point2D((0, 0.09))
#p3 = vm.Point2D((-0.09, 0))
#p4 = vm.Point2D((0, -0.09))
#
#l0 = vm.Arc2D(p1, p2, p3)
#l1 = vm.Arc2D(p3, p4, p1)
#
#primitives = [l0, l1]
#c = vm.Contour2D(primitives)

#p0 = vm.Point2D((5, 0.5))
#p1 = vm.Point2D((-5, 0.5))
#p2 = vm.Point2D((-5, -0.5))
#p3 = vm.Point2D((5, -0.5))
#
#l0 = vm.Line2D(p0, p1)
#l1 = vm.Line2D(p1, p2)
#l2 = vm.Line2D(p2, p3)
#l3 = vm.Line2D(p3, p0)
#
#primitives = [l0, l1, l2, l3]
#c = vm.Contour2D(primitives)

