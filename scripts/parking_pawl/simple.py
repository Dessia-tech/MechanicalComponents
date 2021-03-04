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

# parking_pawl = ParkingPawl(wheel, pawl)
locking_mechanism = mcpp.RollerLockingMechanism(roller_diameter=0.035,
                                                # center_distance=0.062,
                                                width = 0.025)
parking_pawl = mcpp.ParkingPawl(wheel_inner_diameter=0.030,
                                wheel_lower_tooth_diameter=0.060,
                                wheel_outer_diameter=0.080,
                                teeth_number=9,
                                lower_tooth_ratio=0.60,
                                basis_diameter=0.063,
                                contact_diameter=0.070,
                                width = 0.035,
                                pawl_offset = 0.06,
                                axis_inner_diameter=0.025,
                                axis_outer_diameter=0.030,
                                finger_height=0.012,
                                # finger_angle=math.radians(20),
                                finger_width=0.018,
                                slope_start_height=0.015,
                                slope_angle=math.radians(26),
                                slope_offset=0.005, slope_length=0.035,
                                locking_mechanism=locking_mechanism)
# parking_pawl.pawl.outer_contour().plot()
# parking_pawl.wheel.outer_contour().plot()
# plot_data.plot_canvas(plot_data_object=parking_pawl.wheel.plot_data()[0], debug_mode=True)
# parking_pawl.babylonjs()
parking_pawl.mpl_plot()#pawl_angle=parking_pawl.up_pawl_angle, wheel_angle=0)
# parking_pawl.wheel.mpl_plot()
# parking_pawl.pawl.mpl_plot()
simulation = parking_pawl.static_locking_simulation()
print(parking_pawl.check())
simulation.babylonjs()