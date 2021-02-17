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
parking_pawl = mcpp.ParkingPawl(0.030, 0.060, 0.080,
                           9, 0.25, 0.3,
                           basis_diameter=0.030,
                           contact_diameter=0.070,
                           width = 0.035,
                           pawl_offset = 0.06,
                           axis_inner_diameter=0.025,
                           axis_outer_diameter=0.030,
                           finger_height=0.012,
                           finger_angle=math.radians(20), finger_width=0.008,
                           slope_start_height=0.015,
                           slope_angle=math.radians(12),
                           slope_offset=0.005, slope_length=0.035)
# parking_pawl.pawl.outer_contour().plot()
# parking_pawl.wheel.outer_contour().plot()
# plot_data.plot_canvas(plot_data_object=parking_pawl.wheel.plot_data()[0], debug_mode=True)
# parking_pawl.babylonjs()
parking_pawl.mpl_plot()
print(parking_pawl.check())