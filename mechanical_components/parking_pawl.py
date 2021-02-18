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
                 teeth_number:int,
                 upper_tooth_ratio:float, lower_tooth_ratio:float,
                 basis_diameter:float,
                 contact_diameter:float,
                 width:float):
        self.inner_diameter = inner_diameter
        self.lower_tooth_diameter = lower_tooth_diameter
        self.outer_diameter = outer_diameter
        self.upper_tooth_ratio = upper_tooth_ratio
        self.lower_tooth_ratio = lower_tooth_ratio
        self.teeth_number = teeth_number
        self.basis_diameter = basis_diameter
        self.contact_diameter = contact_diameter
        self.width = width
        
        # Computed attributes
        self.tooth_angle = vm.TWO_PI/self.teeth_number
        self.upper_tooth_angle = self.tooth_angle*self.upper_tooth_ratio
        self.lower_tooth_angle = self.tooth_angle*self.lower_tooth_ratio
        self.junction_angle = 0.5*(self.tooth_angle-self.upper_tooth_angle-self.lower_tooth_angle)

        action_point1 = vm.Point2D(0., 0.5 * self.contact_diameter)
        print(math.degrees(math.asin(self.basis_diameter/self.contact_diameter)))
        action_point2 = vm.Point2D(0., 0.5 * self.basis_diameter)\
                            .rotation(vm.O2D, -math.asin(self.basis_diameter/self.contact_diameter))

        self.action_line = vme.Line2D(action_point1, action_point2)
        
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
        for i in range(1, self.teeth_number):
            primitives.extend([arc_upper.rotation(vm.O2D, i*self.tooth_angle),
                                junction1.rotation(vm.O2D, i*self.tooth_angle),
                                arc_lower.rotation(vm.O2D, i*self.tooth_angle),
                                junction2.rotation(vm.O2D, i*self.tooth_angle)])
        return vmw.Contour2D(primitives)

    def inner_contour(self):
        return vmw.Circle2D(vm.O2D, 0.5*self.inner_diameter)

    def plot_data(self):

        primitives = [self.outer_contour().plot_data(),
                      self.inner_contour().plot_data(),
                      self.action_line.plot_data(color='red')]
        return [plot_data.PrimitiveGroup(primitives)]
        
    def volmdlr_primitives(self, frame=vm.OXYZ):
        return [p3d.ExtrudedProfile(frame.origin, frame.v, frame.w,
                                    self.outer_contour(), [inner_contour],
                                    frame.u*self.width)]

        
class Pawl(dc.DessiaObject):
    def __init__(self, axis_position:vm.Point2D,
                 wheel_lower_tooth_diameter:float,
                 axis_inner_diameter:float, axis_outer_diameter:float,
                 finger_height:float, finger_angle:float,
                 finger_width:float,
                 slope_start_height, slope_length,
                 slope_offset:float, slope_angle:float,
                 width:float):
        self.axis_position = axis_position
        self.wheel_lower_tooth_diameter = wheel_lower_tooth_diameter
        self.finger_height = finger_height
        self.finger_angle = finger_angle
        self.finger_width = finger_width
        self.slope_start_height = slope_start_height
        self.slope_length = slope_length
        self.slope_offset = slope_offset
        self.slope_angle = slope_angle
        self.axis_inner_diameter = axis_inner_diameter
        self.axis_outer_diameter = axis_outer_diameter
        self.width = width

    def outer_contour(self):
        
        p1 = vm.Point2D(-0.5*self.finger_width, 0.5*self.wheel_lower_tooth_diameter)
        p2 = vm.Point2D(0.5*self.finger_width, 0.5*self.wheel_lower_tooth_diameter)
        
        p3 = p2 + vm.Point2D(math.sin(self.finger_angle)*self.finger_height,
                             self.finger_height)

        p0 = p1 + vm.Point2D(-math.sin(self.finger_angle)*self.finger_height,
                             self.finger_height)
        
        p4 = p3 + vm.Point2D(self.slope_offset, 0)
        p5 = vm.Point2D(p4.x,
                        0.5*self.wheel_lower_tooth_diameter + self.slope_start_height)
        p6 = p5 - vm.Point2D(self.slope_length, 0.).rotation(vm.O2D, -self.slope_angle)
        
        pa1 = self.axis_position + vm.Point2D(0., 0.5*self.axis_outer_diameter)
        pa2 = self.axis_position + vm.Point2D(-0.5*self.axis_outer_diameter, 0.)
        pa3 = self.axis_position + vm.Point2D(0., -0.5*self.axis_outer_diameter)
        pa3.rotation(self.axis_position, math.radians(60), copy=False)
        
        primitives = ([vme.Arc2D(pa1, pa2, pa3)] + 
                      p2d.OpenedRoundedLineSegments2D([pa3, p0, p1, p2, p3, p4, p5, p6, pa1],
                                                      {}).primitives)
        
        return vmw.Contour2D(primitives)

    def inner_contour(self):
        return vmw.Circle2D(self.axis_position, 0.5*self.axis_inner_diameter)

    def plot_data(self):
        return [plot_data.PrimitiveGroup([self.outer_contour().plot_data(),
                                          self.inner_contour().plot_data()])]

    def volmdlr_primitives(self, frame=vm.OXYZ):
        return [p3d.ExtrudedProfile(frame.origin, frame.v, frame.w,
                                    self.outer_contour(), [self.inner_contour()],
                                    frame.u*self.width)]
    
class ParkingPawl(dc.DessiaObject):
    def __init__(self,
                 wheel_inner_diameter:float, wheel_lower_tooth_diameter:float,
                 outer_diameter:float,
                 teeth_number:int,
                 upper_tooth_ratio:float, lower_tooth_ratio:float,
                 basis_diameter:float,
                 contact_diameter:float,
                 width:float,
                 pawl_offset:float,
                 axis_inner_diameter:float, axis_outer_diameter:float,
                 finger_height:float, finger_angle:float,
                 finger_width:float,
                 slope_start_height, slope_length,
                 slope_offset:float, slope_angle:float,
                 ):


        self.wheel = Wheel(inner_diameter=wheel_inner_diameter,
                           lower_tooth_diameter=wheel_lower_tooth_diameter,
                           outer_diameter=outer_diameter,
                           teeth_number=teeth_number,
                           upper_tooth_ratio=upper_tooth_ratio,
                           lower_tooth_ratio=lower_tooth_ratio,
                           basis_diameter=basis_diameter,
                           contact_diameter=contact_diameter,
                           width=width)
        self.slope_end_height = 0.5*wheel_lower_tooth_diameter + slope_start_height + slope_length*math.sin(slope_angle)

        self.pawl = Pawl(axis_position=vm.Point2D(-pawl_offset, self.slope_end_height-0.5*axis_outer_diameter),
                         wheel_lower_tooth_diameter=wheel_lower_tooth_diameter,
                         axis_inner_diameter=axis_inner_diameter,
                         axis_outer_diameter=axis_outer_diameter,
                         finger_height=finger_height,
                         finger_angle=finger_angle,
                         finger_width=finger_width,
                         slope_start_height=slope_start_height,
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

    def engaged_slack(self):
        return self.wheel.lower_tooth_angle * 0.5 * self.wheel.lower_tooth_diameter - self.pawl.finger_width

    def check(self):
        if self.engaged_slack() < 0:
            return False
        return True

    def plot_data(self):
        red_dashed = plot_data.EdgeStyle(color_stroke=plot_data.colors.RED,
                                         dashline=[3,1,1,3])

        wheel_contour = self.wheel.outer_contour().rotation(vm.O2D, self.engaged_wheel_angle)
        action_line = self.wheel.action_line.rotation(vm.O2D, self.engaged_wheel_angle)
        # basis_circle = plot_data.Circle2D(0, 0, 0.5*self.wheel.basis_diameter, edge_style=red_dashed)
        return [plot_data.PrimitiveGroup([wheel_contour.plot_data(),
                                         self.wheel.inner_contour().plot_data(),
                                         action_line.plot_data(edge_style=red_dashed),
                                         # basis_circle,
                                         ]
                                         +self.pawl.plot_data()[0].primitives)]

