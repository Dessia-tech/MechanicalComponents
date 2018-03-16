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
regime = npy.linspace(0, 2500*math.pi/30, 100)

n_friction_plates = [2, 3]
results = []
for i, n_plates in enumerate(n_friction_plates):
    verin = clutches.HydraulicCylinder()
    emb = clutches.Clutch(verin, separator_tooth_type = 'outer')
    co=clutches.ClutchOptimizer(emb, {'plate_inner_radius' : 0.060,
                                      'plate_height' : (10**-6, 0.060),
                                      'separator_plate_width' : (0.0010, 0.0020),
                                      'friction_plate_width' : (0.0006, 0.0010),
                                      'friction_paper_width' : (0.0004, 0.0008),
                                      'clearance' : (0.00005, 0.00040),
                                      'n_friction_plates' : n_plates,
                                      'hydraulic_cylinder.inner_radius' : (0.010, 0.050),
                                      'hydraulic_cylinder.height' : (0.001, 0.050),
                                      'hydraulic_cylinder.chamber_width' : 0.100,
                                      'hydraulic_cylinder.thickness' : (0.0005, 0.0010), 
                                      'hydraulic_cylinder.engaged_chamber_pressure' : 500000,
                                      'hydraulic_cylinder.spring_wire_diameter' : (0.0005, 0.0015), 
                                      'hydraulic_cylinder.spring_outer_diameter' : (0.005, 0.015), 
                                      'hydraulic_cylinder.spring_free_length' : 0.01,
                                      'hydraulic_cylinder.spring_final_length' : 0.005,
                                      'max_pressure' : 5000000,
                                      'max_time' : 0.2,
                                      'max_drag_torque' : 50})
    res = co.Optimize()
    print(res)
    print('\nr1 : ', co.clutch.plate_inner_radius,
          '\nr2 : ', co.clutch.plate_inner_radius + co.clutch.plate_height,
          '\ne_splate : ', co.clutch.separator_plate_width,
          '\ne_fplate : ', co.clutch.friction_plate_width,
          '\ne_fpaper : ', co.clutch.friction_paper_width,
          '\nh : ', co.clutch.clearance,
          '\nr1_piston : ', co.clutch.hydraulic_cylinder.inner_radius,
          '\nr2_piston : ', co.clutch.hydraulic_cylinder.outer_radius,
          '\nL_chamber : ', co.clutch.hydraulic_cylinder.chamber_width,
          '\ne_chamber : ', co.clutch.hydraulic_cylinder.thickness,
          '\nwire_d : ', co.clutch.hydraulic_cylinder.spring_wire_diameter,
          '\nspring_d : ', co.clutch.hydraulic_cylinder.spring_outer_diameter,
          '\nl0 : ', co.clutch.hydraulic_cylinder.spring_free_length,
          '\nl : ', co.clutch.hydraulic_cylinder.spring_final_length,
          '\nt : ', co.clutch.EngagementTime(),
          '\nTd : ', max(co.clutch.DragTorque(regime, 0)),
          '\n')

    results.append(co.clutch)
    
masses = [i.Mass() for i in results]
n_friction_plates = [i.n_friction_plates for i in results]
transferred_torques = [i.EngagedTransferredTorque() for i in results]
inner_radius = [i.plate_inner_radius for i in results]
outer_radius = [i.plate_outer_radius for i in results]
engagement_time = [i.EngagementTime() for i in results]

a = results[1]
a.CADExport()

plt.plot(engagement_time, masses, 'ro')
plt.show()

test = co.clutch.DragTorque(regime, 0)

plt.plot(regime*30/math.pi, test)
plt.show()