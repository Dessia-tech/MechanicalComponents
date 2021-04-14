#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
import math
import mechanical_components.parking_pawl as mcpp

# wheel = Wheel(0.020, 0.080, 0.060, 9, 0.25, 0.3, 0.035)
# wheel.outer_contour().plot()
# wheel.babylonjs()
# pawl = Pawl(vm.Point2D(-0.06, 0.042), 0.03, 0.025, 0.030,finger_height=0.012,
#             finger_angle=math.radians(20), finger_width=0.008, slope_height_start=0.015,
#             slope_angle=math.radians(12), slope_offset=0.005, slope_length=0.035,
#             width=0.030)

# parking_pawl.mpl_plot()
w = -3/3.6/(0.73/2)*12
C = 0.5*1950*9.81/0.73*math.sin(math.atan(0.3))

# parking_pawl = ParkingPawl(wheel, pawl)
locking_mechanism = mcpp.RollerLockingMechanism(roller_diameter=0.010,
                                                roller_width = 0.025,
                                                spring_stiffness=23000,
                                                spring_active_length=0.011
                                                )

slope_angle = math.radians(40)
parking_pawl = mcpp.ParkingPawl(wheel_inner_diameter=0.030,
                                wheel_lower_tooth_diameter=0.050,
                                wheel_outer_diameter=0.075,
                                teeth_number=7,
                                lower_tooth_ratio=0.60,
                                pressure_angle=math.radians(20.),
                                contact_diameter=0.070,
                                width = 0.025,
                                pawl_offset = 0.04,
                                axis_inner_diameter=0.020,
                                axis_outer_diameter=0.030,
                                finger_height=0.012,
                                roller_rest_length=0.6*locking_mechanism.roller_diameter,
                                finger_width=0.004,
                                slope_start_height=0.015,
                                slope_angle=slope_angle,
                                slope_offset=0.005, slope_length=0.011/math.cos(slope_angle),
                                pawl_spring_stiffness=20,
                                locking_mechanism=locking_mechanism)
# parking_pawl.pawl.outer_contour().plot()
# parking_pawl.wheel.outer_contour().plot()
# plot_data.plot_canvas(plot_data_object=parking_pawl.wheel.plot_data()[0], debug_mode=True)
# parking_pawl.babylonjs()
# parking_pawl.plot()#pawl_angle=parking_pawl.up_pawl_angle, wheel_angle=0)
# parking_pawl.wheel.mpl_plot()
# parking_pawl.pawl.mpl_plot()
# simulation = parking_pawl.static_locking_simulation()
# simulation.babylonjs()
parking_pawl.pawl.size_torsion_spring(3*9.81)
parking_pawl.size_width(wheel_torque=C)
print('rest margin',parking_pawl.rest_margin()*1000, 'mm')
print('engaged_slack', math.degrees(parking_pawl.engaged_slack()), '°')

simulation = parking_pawl.locking_simulation(wheel_speed=w)
print(parking_pawl.check())
simulation.babylonjs()
simulation.plot()
#
# t = abs((parking_pawl.wheel.junction_angle+parking_pawl.wheel.lower_tooth_angle)/w)
#
# print('engaging time margin', t-simulation.time[-1])
# print('speed begin engagment', simulation.max_locking_speed/12*0.73/2*3.6, 'km/h')
#
# simulation._check_platform()