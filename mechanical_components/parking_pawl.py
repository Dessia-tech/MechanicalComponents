#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import math
import volmdlr as vm
import volmdlr.edges as vme
import volmdlr.wires as vmw
import volmdlr.primitives2d as p2d
import volmdlr.primitives3d as p3d
import dessia_common as dc
import plot_data
# import dessia_common.typings as dct


class Wheel(dc.DessiaObject):
    def __init__(self, inner_diameter:float,
                 lower_tooth_diameter:float,
                 outer_diameter:float,
                 number_teeth:int,
                 upper_tooth_ratio:float, lower_tooth_ratio:float, 
                 width:float):
        self.inner_diameter = inner_diameter
        self.lower_tooth_diameter = lower_tooth_diameter
        self.outer_diameter = outer_diameter
        self.upper_tooth_ratio = upper_tooth_ratio
        self.lower_tooth_ratio = lower_tooth_ratio
        self.number_teeth = number_teeth
        self.width = width
        
        # Computed attributes
        self.tooth_angle = vm.TWO_PI/self.number_teeth
        self.upper_tooth_angle = self.tooth_angle*self.upper_tooth_ratio
        self.lower_tooth_angle = self.tooth_angle*self.lower_tooth_ratio
        self.junction_angle = 0.5*(self.tooth_angle-self.upper_tooth_angle-self.lower_tooth_angle)
        
        
    def outer_contour(self):

        
        start_upper = vm.Point2D(0, 0.5*self.outer_diameter)
        interior_upper = start_upper.rotation(vm.O2D, 0.5*self.upper_tooth_angle)
        end_upper = start_upper.rotation(vm.O2D, self.upper_tooth_angle)
        arc_upper = vme.Arc2D(start_upper, interior_upper, end_upper)
        
        start_lower = vm.Point2D(0, 0.5*self.lower_tooth_diameter)
        start_lower.rotation(vm.O2D, self.upper_tooth_angle+self.junction_angle, copy=False)
        interior_lower = start_lower.rotation(vm.O2D, 0.5*self.lower_tooth_angle)
        end_lower = start_lower.rotation(vm.O2D, self.lower_tooth_angle)
        arc_lower = vme.Arc2D(start_lower, interior_lower, end_lower)
        
        junction1 = vme.LineSegment2D(arc_upper.end, arc_lower.start)
        next_tooth_start = start_upper.rotation(vm.O2D, self.tooth_angle)
        junction2 = vme.LineSegment2D(arc_lower.end, next_tooth_start)
        
        primitives = [arc_upper, junction1, arc_lower, junction2]
        for i in range(1, self.number_teeth):
            primitives.extend([arc_upper.rotation(vm.O2D, i*self.tooth_angle),
                                junction1.rotation(vm.O2D, i*self.tooth_angle),
                                arc_lower.rotation(vm.O2D, i*self.tooth_angle),
                                junction2.rotation(vm.O2D, i*self.tooth_angle)])
        return vmw.Contour2D(primitives)
    
    def plot_data(self):
        return [plot_data.PrimitiveGroup([self.outer_contour().plot_data()])]
        
    def volmdlr_primitives(self, frame=vm.OXYZ):
        inner_contour = vmw.Circle2D(vm.O2D, 0.5*self.inner_diameter)
        return [p3d.ExtrudedProfile(frame.origin, frame.v, frame.w,
                                    self.outer_contour(), [inner_contour],
                                    frame.u*self.width)]
        
        
class Pawl(dc.DessiaObject):
    def __init__(self, axis_position:vm.Point2D, contact_radius:float,
                 axis_inner_diameter:float, axis_outer_diameter:float,
                 finger_height:float, finger_angle:float,
                 finger_width:float,
                 slope_height_start, slope_length,
                 slope_offset:float, slope_angle:float,
                 width:float):
        self.axis_position = axis_position
        self.contact_radius = contact_radius
        self.finger_height = finger_height
        self.finger_angle = finger_angle
        self.finger_width = finger_width
        self.slope_height_start = slope_height_start
        self.slope_length = slope_length
        self.slope_offset = slope_offset
        self.slope_angle = slope_angle
        self.axis_inner_diameter = axis_inner_diameter
        self.axis_outer_diameter = axis_outer_diameter
        self.width = width
        
    def outer_contour(self):
        
        p1 = vm.Point2D(-0.5*self.finger_width, self.contact_radius)
        p2 = vm.Point2D(0.5*self.finger_width, self.contact_radius)
        
        p3 = p2 + vm.Point2D(math.sin(self.finger_angle)*self.finger_height,
                             self.finger_height)

        p0 = p1 + vm.Point2D(-math.sin(self.finger_angle)*self.finger_height,
                             self.finger_height)
        
        p4 = p3 + vm.Point2D(self.slope_offset, 0)
        p5 = vm.Point2D(p4.x, self.contact_radius + self.slope_height_start)
        p6 = p5 - vm.Point2D(self.slope_length, 0.).rotation(vm.O2D, -self.slope_angle)
        
        pa1 = self.axis_position + vm.Point2D(0., 0.5*self.axis_outer_diameter)
        pa2 = self.axis_position + vm.Point2D(-0.5*self.axis_outer_diameter, 0.)
        pa3 = self.axis_position + vm.Point2D(0., -0.5*self.axis_outer_diameter)
        pa3.rotation(self.axis_position, math.radians(60), copy=False)
        
        primitives = ([vme.Arc2D(pa1, pa2, pa3)] + 
                      p2d.OpenedRoundedLineSegments2D([pa3, p0, p1, p2, p3, p4, p5, p6, pa1],
                                                      {}).primitives)
        
        return vmw.Contour2D(primitives)
        
    def volmdlr_primitives(self, frame=vm.OXYZ):
        inner_contour = vmw.Circle2D(self.axis_position, 0.5*self.axis_inner_diameter)
        return [p3d.ExtrudedProfile(frame.origin, frame.v, frame.w,
                                    self.outer_contour(), [inner_contour],
                                    frame.u*self.width)]
    
class ParkingPawl(dc.DessiaObject):
    def __init__(self,
                 wheel_inner_diameter:float, wheel_lower_tooth_diameter:float,
                 outer_diameter:float,
                 number_teeth:int,
                 upper_tooth_ratio:float, lower_tooth_ratio:float, 
                 width:float,
                 axis_position:vm.Point2D,
                 axis_inner_diameter:float, axis_outer_diameter:float,
                 finger_height:float, finger_angle:float,
                 finger_width:float,
                 slope_height_start, slope_length,
                 slope_offset:float, slope_angle:float,
                 ):
        self.wheel = Wheel(wheel_inner_diameter, wheel_lower_tooth_diameter,
                           outer_diameter, number_teeth,
                           upper_tooth_ratio, lower_tooth_ratio,
                           width)
        self.pawl = Pawl(axis_position,
                         axis_inner_diameter,
                         wheel_lower_tooth_diameter,
                         axis_outer_diameter,
                         finger_height, finger_angle, finger_width,
                         slope_height_start=slope_height_start,
                         slope_length=slope_length,
                         slope_offset=slope_offset,
                         slope_angle=slope_angle,
                         width=width
                         )
        
        self.engaged_wheel_angle = self.wheel.junction_angle+0.5*self.wheel.lower_tooth_angle
    
    def volmdlr_primitives(self, frame=vm.OXYZ):
        wheel_frame = frame.rotation(frame.origin, frame.u, self.engaged_wheel_angle)
        return (self.wheel.volmdlr_primitives(frame=wheel_frame)
                + self.pawl.volmdlr_primitives(frame=frame))


# wheel = Wheel(0.020, 0.080, 0.060, 9, 0.25, 0.3, 0.035)
# wheel.outer_contour().plot()
# wheel.babylonjs()
# pawl = Pawl(vm.Point2D(-0.06, 0.042), 0.03, 0.025, 0.030,finger_height=0.012,
#             finger_angle=math.radians(20), finger_width=0.008, slope_height_start=0.015,
#             slope_angle=math.radians(12), slope_offset=0.005, slope_length=0.035,
#             width=0.030)

# parking_pawl = ParkingPawl(wheel, pawl)
parking_pawl = ParkingPawl(0.020, 0.060, 0.080,
                           9, 0.25, 0.3, width = 0.035,
                           axis_position = vm.Point2D(-0.06, 0.042),
                           axis_inner_diameter=0.025,
                           axis_outer_diameter=0.030,
                           finger_height=0.012,
                           finger_angle=math.radians(20), finger_width=0.008,
                           slope_height_start=0.015,
                           slope_angle=math.radians(12),
                           slope_offset=0.005, slope_length=0.035)
parking_pawl.pawl.outer_contour().plot()
parking_pawl.wheel.outer_contour().plot()
plot_data.plot_canvas(plot_data_object=parking_pawl.wheel.plot_data()[0], debug_mode=True)
parking_pawl.babylonjs()