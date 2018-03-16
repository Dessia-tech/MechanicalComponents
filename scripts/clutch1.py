import math
import mechanical_components.clutches as clutches
import numpy as npy
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
#regime = npy.linspace(0, 2500*math.pi/30, 100)
#test = emb.DragTorque(regime, 0)
#print(test)
#
#
#plt.plot(regime*30/math.pi, test)
#plt.show()

n_friction_plates = [1, 2, 3, 4, 5, 6, 7]
results = []
for i, n_plates in enumerate(n_friction_plates):
    verin = clutches.HydraulicCylinder()
    emb = clutches.Clutch(verin, separator_tooth_type = 'outer')
    co=clutches.ClutchOptimizer(emb, {'plate_inner_radius' : 0.060,
                                      'plate_outer_radius' : (0.100, 0.120),
                                      'separator_plate_width' : (0.0010, 0.0020),
                                      'friction_plate_width' : (0.0006, 0.0010),
                                      'friction_paper_width' : (0.0004, 0.0008),
                                      'clearance' : (0.00005, 0.00040),
                                      'n_friction_plates' : n_plates,
                                      'max_pressure' : 5000000,
                                      'max_time' : 0.2,
                                      'max_drag_torque' : 50})
    res=co.Optimize()
    print(res)
    [print(val*(co.bounds[i][1] - co.bounds[i][0]) + co.bounds[i][0]) for i, val in enumerate(res.x)]
    results.append(co.clutch)
    
masses = [i.Mass() for i in results]
n_friction_plates = [i.n_friction_plates for i in results]
transferred_torques = [i.ClosedTransferredTorque() for i in results]
inner_radius = [i.plate_inner_radius for i in results]
outer_radius = [i.plate_outer_radius for i in results]
#engagement_time = [i.]

a = results[0]
a.CADExport()

plt.plot(n_friction_plates, masses)
plt.show()

#resp =co.clutch.CADExport()
#print(resp)