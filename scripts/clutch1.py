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

n_friction_plates = [1, 2, 3, 4, 5, 6]
results = []
for i, n_plates in enumerate(n_friction_plates):
    verin = clutches.HydraulicCylinder()
    emb = clutches.Clutch(verin, separator_tooth_type = 'outer')
    co=clutches.ClutchOptimizer(emb, {'plate_inner_radius' : 0.060,
                                      'plate_height' : (0.010, 0.060),
                                      'separator_plate_width' : 0.0018,
                                      'friction_plate_width' : 0.0008,
                                      'friction_paper_width' : 0.0004,
                                      'clearance' : 0.0005,
                                      'n_friction_plates' : n_plates,
                                      'hydraulic_cylinder.inner_radius' : (0.010, 0.050),
                                      'hydraulic_cylinder.height' : (0.010, 0.050),
                                      'hydraulic_cylinder.chamber_width' : 0.100,
                                      'hydraulic_cylinder.thickness' : 0.0010,
                                      'hydraulic_cylinder.engaged_chamber_pressure' : 500000,
                                      'hydraulic_cylinder.spring_wire_diameter' : (0.0005, 0.0015),
                                      'hydraulic_cylinder.spring_outer_diameter' : (0.001, 0.005),
                                      'hydraulic_cylinder.spring_free_length' : 0.01,
                                      'hydraulic_cylinder.spring_final_length' : 0.005,
                                      'max_pressure' : 5000000,
                                      'max_time' : 0.2,
                                      'max_drag_torque' : 50,
                                      'min_torque' : 30})
    res = co.Optimize()
    print('---------------------------------------------\n',
          'Resultat', i+1, '\n\n',
          res)
    print('r1 : ', co.clutch.plate_inner_radius,
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
          '\nt : ', max(co.clutch.EngagementTime(regime)),
          '\nTd : ', max(co.clutch.DragTorque(regime, 0)),
          '\nT : ', co.clutch.EngagedTransferredTorque(),
          '\n---------------------------------------------',
          '\n\n\n')

    if res['success']:
        results.append(co.clutch)
    
masses = [i.Mass() for i in results]
piston_masses = [i.hydraulic_cylinder.Mass() for i in results]
n_friction_plates = [i.n_friction_plates for i in results]
J = [i.J_sep + i.J_fric + i.J_paper for i in results]
transferred_torques = [i.EngagedTransferredTorque() for i in results]
inner_radius = [i.plate_inner_radius for i in results]
outer_radius = [i.plate_outer_radius for i in results]
engagement_time = [max(i.EngagementTime(regime)) for i in results]
pressure = [i.PlatePressure() for i in results]
piston_force = [i.hydraulic_cylinder.PistonForce() for i in results]
spring_force = [i.hydraulic_cylinder.SpringResultingForce() for i in results]
spring_stiffness = [i.hydraulic_cylinder.spring_stiffness for i in results]
piston_height = [i.hydraulic_cylinder.outer_radius - i.hydraulic_cylinder.inner_radius for i in results]
ratio = [i.hydraulic_cylinder.spring_outer_diameter/i.hydraulic_cylinder.spring_wire_diameter for i in results]


G = 210*10**9/(2*(1+0.33))
n = 5
D = [i.hydraulic_cylinder.spring_outer_diameter for i in results]
d = [i.hydraulic_cylinder.spring_wire_diameter for i in results]
stiffness = [G*d[i]**4/(8*n*D[i]**3) for i in range(len(D))]


a = results[1]
a.CADExport()

plt.plot(piston_height, piston_force, 'ro')

test = co.clutch.DragTorque(regime, 0)
#plt.plot(regime*30/math.pi, test)