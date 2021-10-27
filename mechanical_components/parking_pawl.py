#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import math
import random
import numpy as npy
from typing import List
import volmdlr as vm
import volmdlr.edges as vme
import volmdlr.wires as vmw
import volmdlr.faces as vmf
import volmdlr.primitives2d as p2d
import volmdlr.primitives3d as p3d
import dessia_common as dc
import dessia_common.utils as dc_utils
import dessia_common.optimization as dc_opt
import plot_data
import scipy.optimize
import cma
# import dessia_common.typings as dct

class Material(dc.DessiaObject):
    def __init__(self, elastic_limit:float, ultimate_limit:float, shear_elasticity_limit:float, name:str=''):
        self.elastic_limit = elastic_limit
        self.ultimate_limit = ultimate_limit
        self.shear_elasticity_limit = shear_elasticity_limit
        self.name = name

    def bending_stress(self, force:float, width:float, breadth:float, force_lever:float):
        moment = force_lever*force 
        return 6*moment/width/breadth**2

    def shear_stress(self, force:float, width:float, breadth:float):
        return force/width/breadth

    def width_from_shear_force(self, force:float, breadth:float)->float:
        return abs(force)/self.shear_elasticity_limit/breadth

    def width_from_bending_force(self, force:float, breadth:float, force_lever:float)->float:
        # print('force', force)
        # print('breadth', breadth)
        # print('force_lever', force_lever)
        moment = abs(force_lever*force)
        return 6*moment/self.elastic_limit/(breadth**2)



ASTM_A36_STEEL = Material(250e6, 400e6, 0.5*250e6, 'ASTM A36 steel')
ASTM_A514_STEEL = Material(690e6, 760e6, 0.7*690e6, 'ASTM A514 steel')


class Wheel(dc.DessiaObject):
    def __init__(self, inner_diameter:float,
                 lower_tooth_diameter:float,
                 outer_diameter:float,
                 teeth_number:int,
                 # upper_tooth_ratio:float,
                 lower_tooth_ratio:float,
                 # basis_diameter:float,
                 pressure_angle:float,
                 contact_diameter:float,
                 width:float,
                 name:str=''):
        self.inner_diameter = inner_diameter
        self.lower_tooth_diameter = lower_tooth_diameter
        self.outer_diameter = outer_diameter
        # self.upper_tooth_ratio = upper_tooth_ratio
        self.lower_tooth_ratio = lower_tooth_ratio
        self.teeth_number = teeth_number
        # self.basis_diameter = basis_diameter
        self.pressure_angle = pressure_angle
        self.pressure_angle = pressure_angle
        self.contact_diameter = contact_diameter
        self.width = width
        self.name = name



        # Computed attributes
        self.tooth_angle = vm.TWO_PI/self.teeth_number
        # print('tooth_angle', math.degrees(self.tooth_angle))
        # self.upper_tooth_angle = self.tooth_angle*self.upper_tooth_ratio

        # self.junction_angle = 0.5*(self.tooth_angle-self.upper_tooth_angle-self.lower_tooth_angle)

        # self.contact_height_ratio = (self.contact_diameter - self.lower_tooth_diameter)/(self.outer_diameter - self.lower_tooth_diameter)

        self.contact_point1 = vm.Point2D(0, 0.5*self.contact_diameter)

        self.basis_diameter = self.contact_diameter * math.cos(
            self.pressure_angle)

        action_point12 = vm.O2D.rotation(self.contact_point1,
                                          math.asin(self.basis_diameter/self.contact_diameter))

        self.action_line1 = vme.Line2D(self.contact_point1, action_point12)

        self.basis_circle = vmw.Circle2D(vm.O2D, 0.5 * self.basis_diameter)
        self.side_tooth_radius = 4*self.contact_diameter

        side1_center = self.contact_point1 + self.side_tooth_radius * self.action_line1.unit_direction_vector()
        side1_circle = vmw.Circle2D(side1_center, self.side_tooth_radius)

        interior_upper = vm.Point2D(0, 0.5 * self.outer_diameter)
        # end_upper.rotation(vm.O2D, -self.tooth_angle, copy=False)

        start_upper = interior_upper.rotation(vm.O2D,
                                              -self.tooth_angle)
        end_upper = interior_upper.rotation(vm.O2D, self.tooth_angle)
        arc_upper = vme.Arc2D(start_upper, interior_upper, end_upper)

        start_lower = vm.Point2D(0, 0.5 * self.lower_tooth_diameter)
        interior_lower = start_lower.rotation(vm.O2D,
                                              0.5 * self.tooth_angle)
        end_lower = interior_lower.rotation(vm.O2D, self.tooth_angle)
        arc_lower = vme.Arc2D(start_lower, interior_lower, end_lower)



        corrected_upper_end = side1_circle.arc_intersections(arc_upper)[0]
        corrected_lower_start = side1_circle.arc_intersections(arc_lower)[0]

        self.junction_angle = (math.atan2(corrected_lower_start.y,
                                         corrected_lower_start.x)
                               - math.atan2(corrected_upper_end.y,
                                           corrected_upper_end.x))

        self.lower_junction_angle = (math.atan2(corrected_lower_start.y,
                                          corrected_lower_start.x)
                               - math.atan2(self.contact_point1.y,
                                            self.contact_point1.x))




        remaining_angle = self.tooth_angle-2*self.junction_angle
        # print('remaining_angle', math.degrees(remaining_angle))
        if remaining_angle < 0:
            raise ValueError('Negative remaining space: {}'.format(remaining_angle))

        self.lower_tooth_angle = remaining_angle * self.lower_tooth_ratio
        self.upper_tooth_angle = remaining_angle * (1-self.lower_tooth_ratio)


        self.arc_side1 = vme.Arc2D(corrected_upper_end, self.contact_point1,
                              corrected_lower_start)

        self.contact_angle = (2 * self.lower_junction_angle
                                      + self.lower_tooth_angle)

        self.contact_point2 = vm.Point2D(0, 0.5 * self.contact_diameter)
        self.contact_point2.rotation(vm.O2D,
                                     self.contact_angle,
                                     copy=False)

        action_point22 = vm.O2D.rotation(self.contact_point2,
                      -math.asin(self.basis_diameter / self.contact_diameter)
                       # + 2 * self.contact_height_ratio * self.junction_angle)
                                         )

        self.action_line2 = vme.Line2D(self.contact_point2, action_point22)

        side2_center = self.contact_point2 + self.side_tooth_radius * self.action_line2.unit_direction_vector()
        side2_circle = vmw.Circle2D(side2_center, self.side_tooth_radius)


        # next_arc_upper.plot(ax=a, color='g')

        corrected_lower_end = side2_circle.arc_intersections(arc_lower)[0]
        corrected_next_tooth_start = \
            side2_circle.arc_intersections(arc_upper)[0]
        corrected_upper_start = corrected_next_tooth_start.rotation(vm.O2D,
                                                                    -self.tooth_angle,
                                                                    )


        arc_upper = arc_upper.split(corrected_upper_start)[1]
        self.arc_upper = arc_upper.split(corrected_upper_end)[0]

        arc_lower = arc_lower.split(corrected_lower_start)[1]
        self.arc_lower = arc_lower.split(corrected_lower_end)[0]

        self.arc_side2 = vme.Arc2D(corrected_lower_end, self.contact_point2,
                                   corrected_next_tooth_start)


        # a = side1_circle.plot()
        # side2_circle.plot(ax=a, color='grey')
        # self.action_line1.plot(ax=a, color='k')
        # action_point12.plot(ax=a, color='grey')
        # self.action_line2.plot(ax=a, color='grey')
        # arc_upper.plot(ax=a, color='r')
        # arc_lower.plot(ax=a, color='b')
        # self.contact_point1.plot(ax=a)
        # self.contact_point2.plot(ax=a)
        # corrected_next_tooth_start.plot(ax=a, color='grey')
        # corrected_upper_start.plot(ax=a, color='g')
        # corrected_upper_end.plot(ax=a, color='r')
        # corrected_lower_start.plot(ax=a, color='r')
        # corrected_lower_end.plot(ax=a, color='g')
        # self.basis_circle.plot(ax=a, color='grey')
        # a.set_aspect('equal')

    def outer_contour(self):
        # Starting from contact point
        primitives = [self.arc_upper, self.arc_side1, self.arc_lower, self.arc_side2]
        for i in range(1, self.teeth_number):
            primitives.extend([self.arc_upper.rotation(vm.O2D, i*self.tooth_angle),
                               self.arc_side1.rotation(vm.O2D, i*self.tooth_angle),
                               self.arc_lower.rotation(vm.O2D, i*self.tooth_angle),
                               self.arc_side2.rotation(vm.O2D, i*self.tooth_angle)])
        return vmw.Contour2D(primitives)

    def inner_contour(self):
        return vmw.Circle2D(vm.O2D, 0.5*self.inner_diameter)

    def surface2d(self):
        return vmf.Surface2D(self.outer_contour(), [self.inner_contour()])

    def mass(self):
        area = self.surface2d().area()
        return 7800*area*self.width

    def plot_data(self, angle=0.,x=0,y=0):
        line_style = plot_data.EdgeStyle(color_stroke=plot_data.colors.WATERCRESS,
                                         dashline=[5, 5, 20, 5])
        point_translated=vm.Point2D(x,y)
        outer_contour_rotation=self.outer_contour().rotation(vm.O2D, angle)
        inner_contour_rotation= self.inner_contour().rotation(vm.O2D, angle)
        action_line1_rotation=self.action_line1.rotation(vm.O2D, angle)
        action_line2_rotation=self.action_line2.rotation(vm.O2D, angle)

        primitives = [self.basis_circle.translation(point_translated).plot_data(edge_style=line_style),
                      outer_contour_rotation.translation(point_translated).plot_data(),
                      inner_contour_rotation.translation(point_translated).plot_data(),
                      action_line1_rotation.translation(point_translated).plot_data(edge_style=line_style),
                      action_line2_rotation.translation(point_translated).plot_data(edge_style=line_style)]
        return [plot_data.PrimitiveGroup(primitives)]
        
    def volmdlr_primitives(self, frame=vm.OXYZ):
        return [p3d.ExtrudedProfile(frame.origin, frame.v, frame.w,
                                    self.outer_contour(), [self.inner_contour()],
                                    frame.u*self.width,
                                    name='Wheel {}'.format(self.name))]

    def minimum_width(self, wheel_torque:float, safety_factor:float=2):
        material = ASTM_A514_STEEL

        force = wheel_torque/0.5/self.contact_diameter
        l = (self.upper_tooth_angle+self.junction_angle) * 0.5 * self.lower_tooth_diameter
        h = 0.5*(self.contact_diameter - self.lower_tooth_diameter)
        return max(material.width_from_bending_force(force, l, h),
                   material.width_from_shear_force(force, l))*safety_factor

    def max_stress(self, wheel_torque:float):
        material = ASTM_A514_STEEL
        force = wheel_torque*0.5*self.contact_diameter
        l = self.upper_tooth_angle * 0.5 * self.lower_tooth_diameter
        h = 0.5*(self.contact_diameter - self.lower_tooth_diameter)
        return max(material.bending_stress(force, self.width, l, h),
                   material.shear_stress(force, self.width, l))

class TorsionSpring(dc.DessiaObject):
    def __init__(self, torque:float, name:str=''):
        self.torque = torque
        self.name = name
    
class Pawl(dc.DessiaObject):
    def __init__(self, axis_position:vm.Point2D,
                 wheel_lower_tooth_diameter:float,
                 contact_diameter:float,
                 # basis_diameter:float,
                 pressure_angle:float,
                 axis_inner_diameter:float, axis_outer_diameter:float,
                 finger_height:float,
                 # finger_angle:float,
                 finger_width:float,
                 slope_start_height:float, slope_length:float,
                 roller_rest_length:float,
                 slope_offset:float, slope_angle:float,
                 width:float,
                 pawl_spring_stiffness:float,
                 name:str=''):
        self.axis_position = axis_position
        self.wheel_lower_tooth_diameter = wheel_lower_tooth_diameter
        self.contact_diameter = contact_diameter
        # self.basis_diameter = basis_diameter
        self.pressure_angle = pressure_angle
        self.finger_height = finger_height
        # self.finger_angle = finger_angle
        self.finger_width = finger_width
        self.slope_start_height = slope_start_height
        self.roller_rest_length=roller_rest_length
        self.slope_length = slope_length
        self.slope_offset = slope_offset
        self.slope_angle = slope_angle
        self.axis_inner_diameter = axis_inner_diameter
        self.axis_outer_diameter = axis_outer_diameter
        self.width = width
        self.pawl_spring_stiffness = pawl_spring_stiffness
        self.name = name

        self.side_tooth_radius = 10*self.contact_diameter


        pa1 = self.axis_position + vm.Point2D(0., 0.5*self.axis_outer_diameter)
        pa2 = self.axis_position + vm.Point2D(-0.5*self.axis_outer_diameter, 0.)
        pa3 = self.axis_position + vm.Point2D(0., -0.5*self.axis_outer_diameter)
        pa3.rotation(self.axis_position, math.radians(60), copy=False)

        self.axis_arc = vme.Arc2D(pa1, pa2, pa3)


        contact_middle_point = vm.Point2D(0, 0.5*self.contact_diameter)
        # Find minimal contact width
        # self.action_angle = math.asin(self.basis_diameter/self.contact_diameter)
        # self.pressure_angle = math.acos(self.basis_diameter/self.contact_diameter)
        self.basis_diameter = self.contact_diameter * math.cos(self.pressure_angle)
        # print('pressure angle', math.degrees(self.pressure_angle))

        self.min_finger_width = (self.contact_diameter - self.wheel_lower_tooth_diameter)*math.tan(self.pressure_angle)
        # print('self.min_finger_width', self.min_finger_width)
        self.finger_lower_angle =  math.atan(self.finger_width/self.wheel_lower_tooth_diameter)
        # print('finger_lower_angle', math.degrees(self.finger_lower_angle))
        self.finger_side_angle = math.asin((self.contact_diameter-self.wheel_lower_tooth_diameter)/self.wheel_lower_tooth_diameter*math.tan(self.pressure_angle))
        # print('finger_side_angle', math.degrees(self.finger_side_angle))
        # self.contact_angle = math.atan(self.finger_contact_width/self.contact_diameter)
        self.contact_angle = (self.finger_lower_angle + self.finger_side_angle)
        self.finger_contact_width = self.contact_diameter * math.sin(self.contact_angle)
        self.finger_upper_width = 2*(0.5*self.wheel_lower_tooth_diameter+self.finger_height) * math.sin(self.contact_angle)
        # print('contact angle', math.degrees(self.contact_angle))
        self.contact_point1 = contact_middle_point.rotation(vm.O2D, -self.contact_angle)
        self.contact_point2 = contact_middle_point.rotation(vm.O2D, self.contact_angle)
        action_point12 = vm.O2D.rotation(self.contact_point1, 0.5*math.pi - self.pressure_angle)
        action_point22 = vm.O2D.rotation(self.contact_point2, -(0.5*math.pi - self.pressure_angle))

        self.action_line1 = vme.Line2D(self.contact_point1, action_point12)
        self.action_line2 = vme.Line2D(self.contact_point2, action_point22)

        side1_center = self.contact_point1 - self.side_tooth_radius * self.action_line1.unit_direction_vector()
        side1_circle = vmw.Circle2D(side1_center, self.side_tooth_radius)

        side2_center = self.contact_point2 - self.side_tooth_radius * self.action_line2.unit_direction_vector()
        side2_circle = vmw.Circle2D(side2_center, self.side_tooth_radius)


        lower_finger_line = vme.Line2D(vm.Point2D(0, 0.5*self.wheel_lower_tooth_diameter),
                                       vm.Point2D(1, 0.5*self.wheel_lower_tooth_diameter))
        upper_finger_line = vme.Line2D(vm.Point2D(0, 0.5*self.wheel_lower_tooth_diameter+self.finger_height),
                                       vm.Point2D(1, 0.5*self.wheel_lower_tooth_diameter+self.finger_height))


        # ax = side1_circle.plot(color='grey')
        # vmw.Circle2D(vm.O2D, 0.5*self.basis_diameter).plot(color='grey', ax=ax)
        # side2_circle.plot(ax=ax, color='grey')
        # self.action_line1.plot(ax=ax, color='grey')
        # self.action_line2.plot(ax=ax, color='grey')
        # action_point12.plot(ax=ax, color='g')
        # action_point22.plot(ax=ax, color='g')
        # self.contact_point1.plot(ax=ax)
        # self.contact_point2.plot(ax=ax)
        # lower_finger_line.plot(ax=ax)
        # upper_finger_line.plot(ax=ax)
        #
        # cp11 = contact_middle_point.rotation(vm.O2D, -self.finger_lower_angle)
        # cp11.plot(ax=ax, color='y')
        # cp12 = contact_middle_point.rotation(vm.O2D, -(self.finger_lower_angle+self.finger_side_angle))
        # cp12.plot(ax=ax, color='y')
        # dl1 = vme.LineSegment2D(vm.O2D, cp11)
        # dl1.plot(ax=ax, color='y', alpha=0.5)
        # dl2 = vme.LineSegment2D(vm.O2D, cp12)
        # dl2.plot(ax=ax, color='y', alpha=0.5)


        # ax.set_xlim(-0.5*self.wheel_lower_tooth_diameter, 0.5*self.wheel_lower_tooth_diameter)
        # ax.set_ylim(0.35*self.wheel_lower_tooth_diameter, 1.5*self.wheel_lower_tooth_diameter)
        # ax.set_aspect('equal')

        side1_start  = sorted(side1_circle.line_intersections(lower_finger_line), key=lambda p:p.x)[1]
        side1_end  = sorted(side1_circle.line_intersections(upper_finger_line), key=lambda p:p.x)[1]

        side2_end  = sorted(side2_circle.line_intersections(lower_finger_line), key=lambda p:p.x)[0]
        side2_start  = sorted(side2_circle.line_intersections(upper_finger_line), key=lambda p:p.x)[0]

        # side1_start.plot(ax=ax, color='r')
        # side1_end.plot(ax=ax, color='g')
        # side2_start.plot(ax=ax, color='g')
        # side2_end.plot(ax=ax, color='r')

        self.arc_side1 = vme.Arc2D(side1_start, self.contact_point1, side1_end)
        self.arc_side2 = vme.Arc2D(side2_start, self.contact_point2, side2_end)

        self.lower_finger_line = vme.LineSegment2D(side2_end, side1_start)

        self.junction1 = vme.LineSegment2D(pa3, side2_start)

        finger_lower_end1 = side1_end + vm.Point2D(self.roller_rest_length+self.slope_offset, 0.)
        finger_lower_end2 = vm.Point2D(finger_lower_end1.x,
                                       self.slope_start_height+0.5*self.wheel_lower_tooth_diameter)
        slope_start = vm.Point2D(side1_end.x+self.slope_offset, self.slope_start_height+0.5*self.wheel_lower_tooth_diameter)
        self.junction2 = vme.LineSegment2D(side1_end, finger_lower_end1)
        self.junction3 = vme.LineSegment2D(finger_lower_end1, finger_lower_end2)
        self.roller_rest = vme.LineSegment2D(finger_lower_end2, slope_start)

        slope_end = slope_start - vm.Point2D(self.slope_length, 0.).rotation(vm.O2D, -self.slope_angle)
        self.slope = vme.LineSegment2D(slope_start, slope_end)
        self.junction4 = vme.LineSegment2D(slope_end, pa1)

        # ax = self.roller_rest.plot()
        # self.junction2.plot(ax=ax, color='grey')
        # self.junction3.plot(ax=ax, color='g')
        # self.slope.plot(ax=ax, color='r')
        # self.junction4.plot(ax=ax, color='b')

        self.torsion_spring = None

    def profile(self):
        radius = 0.2*self.slope_length
        profile = p2d.OpenedRoundedLineSegments2D([self.roller_rest.start,
                                                   self.slope.start,
                                                   self.junction4.start,
                                                   self.junction4.end],
                                                  {1:radius, 2:radius}
                                                  )
        return profile

    def outer_contour(self):
        
        # p1 = vm.Point2D(-0.5*self.finger_width, 0.5*self.wheel_lower_tooth_diameter)
        # p2 = vm.Point2D(0.5*self.finger_width, 0.5*self.wheel_lower_tooth_diameter)
        
        # p3 = p2 + vm.Point2D(math.sin(self.finger_angle)*self.finger_height,
        #                      self.finger_height)

        # p0 = p1 + vm.Point2D(-math.sin(self.finger_angle)*self.finger_height,
        #                      self.finger_height)
        
        # p4 = p3 + vm.Point2D(self.slope_offset, 0)
        # p5 = vm.Point2D(p4.x,
        #                 0.5*self.wheel_lower_tooth_diameter + self.slope_start_height)
        # p6 = p5 - vm.Point2D(self.slope_length, 0.).rotation(vm.O2D, -self.slope_angle)

        profile = self.profile()
        
        primitives = [self.axis_arc, self.junction1, self.arc_side2,
                      self.lower_finger_line, self.arc_side1, self.junction2,
                      self.junction3,# self.slope, self.junction4]
                      ]+profile.primitives
        return vmw.Contour2D(primitives)

    def inertia(self):
        Ix, Iy, _ = self.surface2d().second_moment_area(self.axis_position)
        return 7800*(Ix+Iy)*self.width

    def mass(self):
        area = self.surface2d().area()
        return 7800*area*self.width


    def inner_contour(self):
        return vmw.Circle2D(self.axis_position, 0.5*self.axis_inner_diameter)

    def surface2d(self):
        return vmf.Surface2D(self.outer_contour(), [self.inner_contour()])

    def plot_data(self, angle=0.,x=0,y=0):
        line_style = plot_data.EdgeStyle(color_stroke=plot_data.colors.FERN,
                                         dashline=[5, 5, 20, 5])
        point_translated=vm.Point2D(x,y)
        outer_contour_rotation=self.outer_contour().rotation(self.axis_position, angle)
        inner_contour_rotation= self.inner_contour().rotation(self.axis_position, angle)
        action_line1_rotation=self.action_line1.rotation(self.axis_position, angle)
        action_line2_rotation=self.action_line2.rotation(self.axis_position, angle)
        
        primitives = [outer_contour_rotation.translation(point_translated).plot_data(),
                     inner_contour_rotation.translation(point_translated).plot_data(),
                     action_line1_rotation.translation(point_translated).plot_data(edge_style=line_style),
                     action_line2_rotation.translation(point_translated).plot_data(edge_style=line_style)]
        
        

        return [plot_data.PrimitiveGroup(primitives)]

    def volmdlr_primitives(self, frame=vm.OXYZ):
        return [p3d.ExtrudedProfile(frame.origin, frame.v, frame.w,
                                    self.outer_contour(), [self.inner_contour()],
                                    frame.u*self.width,
                                    name='Pawl {}'.format(self.name))]

    def size_torsion_spring(self, max_acceleration:float):
        cog = self.surface2d().center_of_mass()
        max_lever = cog.point_distance(self.axis_position)
        max_force = self.mass() * max_acceleration
        torque_up = max_force * max_lever


        self.torsion_spring = TorsionSpring(torque_up)

    def minimum_width(self, wheel_torque:float, safety_factor:float=2.):
        material = ASTM_A514_STEEL
        force = wheel_torque/0.5/self.contact_diameter
        l = self.finger_contact_width
        h = 0.5*(self.wheel_lower_tooth_diameter+self.finger_height-self.contact_diameter)
        return max(material.width_from_bending_force(force, l, h),
                   material.width_from_shear_force(force, l)) * safety_factor

    def max_stress(self, wheel_torque:float):
        material = ASTM_A514_STEEL
        force = wheel_torque*0.5*self.contact_diameter
        l = self.finger_contact_width
        h = 0.5*(self.wheel_lower_tooth_diameter+self.finger_height-self.contact_diameter)

        return max(material.bending_stress(force, self.width, l, h),
                   material.shear_stress(force, self.width, l))

class RollerLockingMechanism(dc.DessiaObject):
    _standalone_in_db = True

    def __init__(self, roller_diameter:float,
                 roller_width:float,
                 spring_stiffness:float,
                 spring_active_length:float,
                 name:str=''):
        # self.start_position = start_position
        # self.end_position = end_position
        # self.center_distance = center_distance
        self.roller_diameter = roller_diameter
        self.roller_width = roller_width
        self.spring_stiffness = spring_stiffness
        self.spring_active_length = spring_active_length
        self.name = name

    def contacts_from_profile(self, profile:vmw.Wire2D, center_distance:float):
        # ax = slope.plot()
        locking_mechanism_line = vme.Line2D(vm.Point2D(0, center_distance),
                                            vm.Point2D(0.1, center_distance))


        offset_profile = profile.offset(-0.5*self.roller_diameter)

        # ax = offset_profile.plot(color='grey')
        # profile.plot(ax=ax)
        # locking_mechanism_line.plot(ax=ax, color='grey')

        inters = sorted(offset_profile.line_intersections(locking_mechanism_line),
                                  key=lambda p:p[0].x)

        contacts = []
        for roller_center, contact_offset_edge in inters:

            contact_edge = profile.primitives[offset_profile.primitives.index(contact_offset_edge)]

            contact_normal = contact_offset_edge.normal_vector(contact_offset_edge.abscissa(roller_center))
            contact_point = roller_center + 0.5*self.roller_diameter*contact_normal

            # contact_point.plot(ax=ax, color='r')
            # roller_center.plot(ax=ax, color='b')

            contact_line = vme.LineSegment2D(contact_point, contact_point)
            contacts.append([roller_center.x, contact_point, contact_normal, contact_edge])

        return contacts

    def spring_force(self, position):
        if position > self.spring_active_length:
            return 0.
        else:

            return self.spring_stiffness*(self.spring_active_length - position)

    def outer_contour(self):
        return vmw.Circle2D(vm.Point2D(0., 0.),
                            0.5 * self.roller_diameter)
    def inner_contour(self):
        return vmw.Circle2D(vm.Point2D(0., 0.),
                            0.25 * self.roller_diameter)

    def plot_data(self, center_distance=0., position=0.,x=0,y=0):
        # line_style = plot_data.EdgeStyle(color_stroke=plot_data.colors.FERN,
        #                                  dashline=[5, 5, 20, 5])
        primitives = [self.outer_contour().translation(vm.Point2D(position+x, center_distance+y)).plot_data(),
                      self.inner_contour().translation(vm.Point2D(position+x, center_distance+y)).plot_data(),
                      self.outer_contour().translation(
                          vm.Point2D(position+x, center_distance+self.roller_diameter+y)).plot_data(),
                      self.inner_contour().translation(
                          vm.Point2D(position+x, center_distance+self.roller_diameter+y)).plot_data(),

                      ]
        
        

        return [plot_data.PrimitiveGroup(primitives)]

    def volmdlr_primitives(self, frame=vm.OXYZ):

        roller1 = p3d.HollowCylinder(position=frame.origin+0.5*self.roller_width*frame.u,
                                     axis=frame.u,
                                     inner_radius=0.25*self.roller_diameter,
                                     outer_radius=0.5*self.roller_diameter,
                                     length=self.roller_width,
                                     name='Roller 1 {}'.format(self.name))
        roller2 = p3d.HollowCylinder(position=(frame.origin
                                               + 0.5*self.roller_width*frame.u
                                               + self.roller_diameter*frame.w),
                                     axis=frame.u,
                                     inner_radius=0.25*self.roller_diameter,
                                     outer_radius=0.5*self.roller_diameter,
                                     length=self.roller_width,
                                     name='Roller 2 {}'.format(self.name))
        return [roller1, roller2]


class ParkingPawl(dc.DessiaObject):
    _standalone_in_db = True
    _allowed_methods = ['static_locking_simulation', 'locking_simulation']

    def __init__(self,
                 wheel_inner_diameter:float,
                 wheel_lower_tooth_diameter:float,
                 wheel_outer_diameter:float,
                 teeth_number:int,
                 lower_tooth_ratio:float,
                 contact_diameter:float,
                 pressure_angle:float,
                 width:float,
                 pawl_offset:float,
                 axis_inner_diameter:float, axis_outer_diameter:float,
                 finger_height:float,
                 finger_width:float,
                 roller_rest_length:float,
                 slope_start_height:float,
                 slope_length:float,
                 slope_offset:float,
                 slope_angle:float,
                 pawl_spring_stiffness:float,
                 locking_mechanism:RollerLockingMechanism,
                 open_clearance:float=0.002,
                 name:str=''
                 ):

        self.wheel_inner_diameter = wheel_inner_diameter
        self.wheel_lower_tooth_diameter = wheel_lower_tooth_diameter
        self.wheel_outer_diameter = wheel_outer_diameter
        self.teeth_number = teeth_number
        self.lower_tooth_ratio = lower_tooth_ratio
        self.pressure_angle = pressure_angle
        self.contact_diameter = contact_diameter
        self.width = width
        self.pawl_offset = pawl_offset
        self.axis_inner_diameter = axis_inner_diameter
        self.axis_outer_diameter = axis_outer_diameter
        self.finger_height = finger_height
        self.finger_width = finger_width
        self.roller_rest_length = roller_rest_length
        self.slope_start_height = slope_start_height
        self.slope_length = slope_length
        self.slope_offset = slope_offset
        self.slope_angle = slope_angle
        self.pawl_spring_stiffness = pawl_spring_stiffness

        self.locking_mechanism = locking_mechanism
        self.open_clearance = open_clearance
        self.name = name

        self._utd_locking_contact_results = False

        self.wheel = Wheel(inner_diameter=wheel_inner_diameter,
                           lower_tooth_diameter=wheel_lower_tooth_diameter,
                           outer_diameter=wheel_outer_diameter,
                           teeth_number=teeth_number,
                           # upper_tooth_ratio=upper_tooth_ratio,
                           lower_tooth_ratio=lower_tooth_ratio,
                           pressure_angle=pressure_angle,
                           contact_diameter=contact_diameter,
                           width=width)
        self.slope_end_height = 0.5*wheel_lower_tooth_diameter + slope_start_height + slope_length*math.sin(slope_angle)

        self.pawl = Pawl(axis_position=vm.Point2D(-pawl_offset, self.slope_end_height-0.5*axis_outer_diameter),
                         wheel_lower_tooth_diameter=wheel_lower_tooth_diameter,
                         contact_diameter=contact_diameter,
                         pressure_angle=pressure_angle,
                         axis_inner_diameter=axis_inner_diameter,
                         axis_outer_diameter=axis_outer_diameter,
                         finger_height=finger_height,
                         roller_rest_length=roller_rest_length,
                         finger_width=finger_width,
                         slope_start_height=slope_start_height,
                         slope_length=slope_length,
                         slope_offset=slope_offset,
                         slope_angle=slope_angle,
                         width=width,
                         pawl_spring_stiffness=pawl_spring_stiffness
                         )



        self.contact1_wheel_angle = -0.5*self.pawl.contact_angle
        self.contact2_wheel_angle = -(self.wheel.lower_tooth_angle
                                     +2*self.wheel.lower_junction_angle
                                     -0.5*self.pawl.contact_angle)

        R = self.pawl.arc_side2.end.point_distance(self.pawl.axis_position)
        x, y = self.pawl.arc_side2.end-self.pawl.axis_position
        alpha = math.atan2(y, x)

        self.minimal_up_pawl_angle = -(math.asin(math.sin(alpha)
                               -(0.5*(wheel_outer_diameter
                                      -wheel_lower_tooth_diameter)
                                 +open_clearance)/R)
                               - alpha)

        # Finding locking mech center distance
        roller_rest_contact_point = self.pawl.roller_rest.point_at_abscissa(self.pawl.roller_rest.length()-0.7*self.locking_mechanism.roller_diameter)
        y, z = self.pawl.roller_rest.start.rotation(self.pawl.axis_position,
                                                    self.minimal_up_pawl_angle)
        self.locking_mechanism_center_distance = self.pawl.slope.end.y + 0.5*self.locking_mechanism.roller_diameter
        self.locking_mechanism_start_position = self.locking_contact_results[1][-1]
        self.locking_mechanism_end_position = self.pawl.slope.end.x - 0.002
        
        dc.DessiaObject.__init__(self, name=name)

    def mass(self):
        return abs(self.wheel.mass() + self.pawl.mass())

    def minimum_width(self, wheel_torque:float, safety_factor:float=2):
        wheel_min_width = self.wheel.minimum_width(wheel_torque)
        pawl_min_width = self.pawl.minimum_width(wheel_torque)
        return max(wheel_min_width,
                   pawl_min_width)

    def ejection_levers(self):
        # action_line1 = self.wheel.action_line1.rotation(vm.O2D, self.contact1_wheel_angle)
        # action_line2 = self.wheel.action_line2.rotation(vm.O2D, self.contact2_wheel_angle)
        u1 = self.pawl.action_line1.unit_direction_vector()
        u2 = self.pawl.action_line1.unit_direction_vector()
        
        
        # ax = self.pawl.outer_contour().plot()
        # self.wheel.outer_contour().plot(ax=ax)
        # action_line1.plot(ax=ax, color='r')
        # action_line2.plot(ax=ax, color='r')
        # self.pawl.axis_position.plot(color='b', ax=ax)
        
        lever1 = (self.pawl.axis_position - self.pawl.contact_point1).cross(u1)
        lever2 = (self.pawl.axis_position - self.pawl.contact_point2).cross(u2)
        return lever1, lever2


    def functional_footprint(self):
        xmin = self.pawl.axis_position.x - 0.5* self.pawl.axis_outer_diameter
        ymax = self.locking_mechanism_center_distance + 1.5*self.locking_mechanism.roller_diameter
        ymin = -0.5 * self.wheel.outer_diameter
        xmax = self.locking_mechanism_start_position + 0.5*self.locking_mechanism.roller_diameter 

        return xmin, xmax, ymin, ymax

    def volmdlr_primitives(self, frame=vm.OXYZ):
        wheel_frame = frame.rotation(frame.origin, frame.u, self.contact1_wheel_angle)
        locking_mech_frame = frame.copy()
        locking_mech_frame.origin.z = locking_mech_frame.origin.z+self.locking_mechanism_center_distance
        
        return (
                self.wheel.volmdlr_primitives(frame=wheel_frame)
                + self.pawl.volmdlr_primitives(frame=frame)
                + self.locking_mechanism.volmdlr_primitives(frame=locking_mech_frame))

    def engaged_slack(self):
        """
        in radians
        """
        return self.wheel.contact_angle - 2*self.pawl.contact_angle

    def axis_wheel_clearance(self):
        center_distance = self.pawl.axis_position.point_distance(vm.O2D)
        return center_distance - 0.5*self.wheel_outer_diameter - 0.5*self.pawl.axis_outer_diameter

    def rest_margin(self):
        return (self.pawl.slope_length*math.sin(self.pawl.slope_angle-self.minimal_up_pawl_angle)
                - self.pawl.junction4.length()*math.sin(self.minimal_up_pawl_angle))

    def check(self, clearance:float=0.003):
        if self.engaged_slack() < 0:
            return False
        if self.pawl.slope_angle < self.minimal_up_pawl_angle:
            return False

        if self.axis_wheel_clearance() < clearance:
            return False

        if self.rest_margin() < 0:
            return False

        if min(self.ejection_levers()) < 0:
            return False

        return True

    def plot_data(self, pawl_angle=None, wheel_angle=None,x=0,y=0):

        if pawl_angle is None and wheel_angle is None:
            primitives_p1 = self.locking_mechanism.plot_data(position=self.locking_mechanism_end_position,
                                                             center_distance=self.locking_mechanism_center_distance,
                                                             x=x,y=y)
            primitives_p1 += self.wheel.plot_data(angle=self.contact1_wheel_angle,x=x,y=y)[0].primitives
            primitives_p1 += self.pawl.plot_data(angle=0.,x=x,y=y)[0].primitives
            
            primitives_p2 = self.locking_mechanism.plot_data(position=self.locking_mechanism_end_position,
                                                             center_distance=self.locking_mechanism_center_distance,
                                                             x=x,y=y)
            primitives_p2 += self.wheel.plot_data(angle=self.contact2_wheel_angle,x=x,y=y)[0].primitives
            primitives_p2 += self.pawl.plot_data(angle=0.,x=x,y=y)[0].primitives

            primitives_p3 = self.locking_mechanism.plot_data(position=self.locking_mechanism_start_position,
                                                             center_distance=self.locking_mechanism_center_distance,
                                                             x=x,y=y)
            primitives_p3 += self.wheel.plot_data(angle=0.,x=x,y=y)[0].primitives
            primitives_p3 += self.pawl.plot_data(angle=self.minimal_up_pawl_angle,x=x,y=y)[0].primitives

            # An example of Graph2D instantiation. It draws one or several datasets on
            # one canvas and is useful for displaying numerical functions

            to_disp_attribute_names = ['lock position', 'pawl angle']
            tooltip = plot_data.Tooltip(to_disp_attribute_names=to_disp_attribute_names)
            pawl_angles, travels, _, _ = self.locking_contact_results
            elements = []
            for pawl_angle, travel in zip(pawl_angles, travels):
                elements.append({'pawl angle': pawl_angle,
                                 'lock position': travel})

            # The previous line instantiates a dataset with limited arguments but
            # several customizations are available
            point_style = plot_data.PointStyle(color_fill=plot_data.colors.RED,
                                               color_stroke=plot_data.colors.BLACK)
            edge_style = plot_data.EdgeStyle(color_stroke=plot_data.colors.BLUE,
                                             dashline=[10, 5])

            dataset = plot_data.Dataset(elements=elements,
                                        name='travel pawl',
                                        tooltip=tooltip,
                                        point_style=point_style,
                                        edge_style=edge_style)

            return [plot_data.PrimitiveGroup(primitives_p1),
                    plot_data.PrimitiveGroup(primitives_p2),
                    plot_data.PrimitiveGroup(primitives_p3),
                    plot_data.Graph2D(graphs=[dataset],
                                      to_disp_attribute_names=to_disp_attribute_names)]
        else:
            if pawl_angle is None:
                pawl_angle = 0.
            if wheel_angle is None:
                wheel_angle = 0.
            primitives_p1 = self.locking_mechanism.plot_data()
            primitives_p1 += self.pawl.plot_data(angle=pawl_angle)[0].primitives
            primitives_p1 += self.wheel.plot_data(angle=wheel_angle)[
                0].primitives
            return [plot_data.PrimitiveGroup(primitives_p1)]

    def _solve_locking_contact(self, angular_resolution:float = 0.01, max_steps=40):
        max_angle = min(3*self.minimal_up_pawl_angle, math.radians(80.))
        number_steps = min(math.ceil((max_angle/angular_resolution)), max_steps)

        contacts = [(0, self.pawl.profile().primitives[-1].start.x,
                     self.pawl.profile().primitives[-1].start,
                     -vm.Y2D,
                     self.pawl.profile().primitives[-1])]
        for i in range(number_steps):
            pawl_angle = (i+1) * max_angle / number_steps
            #
            # slope = self.pawl.slope.rotation(self.pawl.axis_position, pawl_angle)
            pawl_profile = self.pawl.profile().rotation(self.pawl.axis_position, pawl_angle)
            step_contacts = self.locking_mechanism.contacts_from_profile(pawl_profile,self.locking_mechanism_center_distance)
            for step_contact in step_contacts:
                step_contact.insert(0, pawl_angle)
                contacts.append(step_contact)

        sorted_contacts = sorted(contacts, key=lambda c:c[1])
        last_travel = None
        pawl_angles = []
        locking_travels =[]
        contact_points = []
        contact_normals = []
        for pawl_angle, travel, contact_point, contact_normal, contact_edge in sorted_contacts:
            if travel != last_travel and travel > self.pawl.axis_position.x:
                pawl_angles.append(pawl_angle)
                locking_travels.append(travel)
                contact_normals.append(contact_normal)
                contact_points.append(contact_point)
                last_travel = travel
        return pawl_angles, locking_travels, contact_points, contact_normals

    @property
    def locking_contact_results(self):
        if not self._utd_locking_contact_results:
            self._locking_contact_results = self._solve_locking_contact()
            self._utd_locking_contact_results = True
        return self._locking_contact_results

    def pawl_angle_from_locking_position(self, locking_position:float):
        if locking_position < self.locking_contact_results[1][0]:
            return 0.
        pawl_angles, locking_positions, _, _ = self.locking_contact_results
        istep = dc_utils.istep_from_value_on_list(locking_positions, locking_position)
        return dc_utils.interpolate_from_istep(pawl_angles, istep)

    def locking_position_from_pawl_angle(self, pawl_angle:float):
        pawl_angles, locking_positions, _, _ = self.locking_contact_results
        istep = dc_utils.istep_from_value_on_list(pawl_angles, pawl_angle, extrapolate=True)
        return dc_utils.interpolate_from_istep(locking_positions, istep)


    def contact_point_from_locking_position(self, locking_position:float):
        if locking_position < self.locking_contact_results[1][0]:
            return self.locking_contact_results[2][0]
        pawl_angles, locking_positions, contact_points, _ = self.locking_contact_results
        istep = dc_utils.istep_from_value_on_list(locking_positions, locking_position)
        return dc_utils.interpolate_from_istep(contact_points, istep)


    def contact_normal_from_locking_position(self, locking_position:float):
        if locking_position < self.locking_contact_results[1][0]:
            return self.locking_contact_results[3][0]

        pawl_angles, locking_positions, _, contact_normals = self.locking_contact_results
        istep = dc_utils.istep_from_value_on_list(locking_positions, locking_position)
        return dc_utils.interpolate_from_istep(contact_normals, istep)


    def static_locking_simulation(self, distance_step=0.002):
        wheel_angles = []
        pawl_angles = []
        locking_positions = []
        time = []
        locking_forces = []

        travel = self.locking_mechanism_end_position - self.locking_mechanism_start_position
        number_steps = round(abs(travel)/distance_step) + 1
        distance_step = travel/distance_step


        for i in range(number_steps+1):
            time.append(i)
            locking_position = self.locking_mechanism_start_position + i/(number_steps)*travel
            wheel_angles.append(0.)
            locking_positions.append(locking_position)
            pawl_angles.append(self.pawl_angle_from_locking_position(locking_position))
            locking_forces.append(0.)


        return ParkingPawlSimulation(self, time, wheel_angles, pawl_angles,
                                     locking_positions, locking_forces,
                                     name='Locking simulation')

    def locking_simulation(self, wheel_speed=0., initial_time_step=0.001,
                           min_step_number=10)->'ParkingPawlSimulation':

        locking_position = self.locking_mechanism_start_position

        if wheel_speed > 0:
            wheel_angle = -(self.wheel.junction_angle+self.wheel.lower_junction_angle+self.wheel.lower_tooth_angle)
        elif wheel_speed < 0:
            wheel_angle = self.contact2_wheel_angle-(self.wheel.junction_angle+self.wheel.lower_junction_angle+self.wheel.lower_tooth_angle)
        else:
            wheel_angle = 0.

        # print('wheel_angle', wheel_angle)

        wheel_angles = [wheel_angle]
        initial_angle = self.pawl_angle_from_locking_position(locking_position)
        pawl_angles = [initial_angle]
        locking_positions = [locking_position]
        locking_forces = [0.]
        time = [0.]

        pawl_angle = initial_angle

        # distance_step = travel/distance_step
        # print(travel)
        t = 0.
        inertia = self.pawl.inertia()
        # print('inertia', inertia)
        angular_speed = 0.
        
        time_step = initial_time_step
        max_delta_angle = self.minimal_up_pawl_angle / min_step_number
        
        n_step = 1
        last_step_step_change = 0
        if not self.pawl.torsion_spring:
            raise NotImplementedError('No torsion spring defined')
        # ax = self.pawl.surface2d().plot()
        # ax = None
        while (n_step < min_step_number*10) and (locking_position > self.locking_mechanism_end_position):
            # print('====\nStep n°{} time step: {}s @{}s'.format(n_step, time_step, t))
            # print('locking_position', locking_position, locking_position > self.locking_mechanism_end_position)

            contact_point = self.contact_point_from_locking_position(locking_position)
            contact_normal = self.contact_normal_from_locking_position(locking_position)

            spring_force = self.locking_mechanism.spring_force(
                self.locking_mechanism_start_position-locking_position)
            locking_force = -spring_force*vm.X2D.dot(contact_normal)
            # print('force transmission ratio', vm.X2D.dot(contact_normal))
            # print('normal', contact_normal)
            contact_line = vme.LineSegment2D(contact_point, contact_point+0.05*contact_normal)
            # ax = contact_line.plot(ax= ax)
            # print('contact_normal', contact_normal, contact_point)
            # contact_normal.plot(ax=ax)
            # contact_point.plot(ax=ax)

            # print('torsion spring torque: ', self.pawl.torsion_spring.torque)
            lever = (contact_point - self.pawl.axis_position).cross(contact_normal)
            locking_moment = locking_force * lever

            sum_moments =  locking_moment + self.pawl.torsion_spring.torque
            # print('sum_moments', sum_moments)
            # print('lever', lever)
            if sum_moments > 0 and (sum_moments * lever < 0):
                # Contact orientation is not well positioned to push
                # Roller will slide
                # print('pushing roller')
                locking_position -= 0.0005
                t += 0.0005
                if locking_position < self.locking_mechanism_end_position:
                    locking_position = self.locking_mechanism_end_position
                    break
                pawl_angle = self.pawl_angle_from_locking_position(locking_position)
                # print('pawl_angle', pawl_angle)

                wheel_angles.append(wheel_angle)
                locking_forces.append(locking_force)
                locking_positions.append(locking_position)
                pawl_angles.append(pawl_angle)
                time.append(t)

                continue

            angular_acceleration = sum_moments/inertia
            delta_pawl_angle = 0.5*angular_acceleration*time_step**2 + angular_speed*time_step
            angular_speed += angular_acceleration*time_step

            # print('angle ratio: ', abs(delta_pawl_angle)/max_delta_angle)

            if last_step_step_change != n_step:
                if abs(delta_pawl_angle) > max_delta_angle:
                    time_step /= (2*abs(delta_pawl_angle)/max_delta_angle)
                    # print('reducing time step to ', time_step)
                    last_step_step_change = n_step
                    delta_pawl_angle = 0.5 * angular_acceleration * time_step ** 2 + angular_speed * time_step
                    angular_speed += angular_acceleration * time_step
                    continue

                if abs(delta_pawl_angle) < 0.2 * max_delta_angle:
                    time_step *= 1.5
                    # print('increasing time step to ', time_step)
                    last_step_step_change = n_step
                    delta_pawl_angle = 0.5 * angular_acceleration * time_step ** 2 + angular_speed * time_step
                    angular_speed += angular_acceleration * time_step
                    continue

            # if delta_pawl_angle < 0.01*max_delta_angle:
            #     time_step *= 1.5
            #     print('increasing time step to ', time_step)
            #     continue

            n_step += 1

            pawl_angle += delta_pawl_angle
            angular_speed += angular_speed * time_step
            wheel_angle += wheel_speed*time_step
            t += time_step
            
            locking_position = self.locking_position_from_pawl_angle(pawl_angle)
            # locking_position = self.locking_mechanism_start_position + i/(number_steps)*travel# CHANGE!

            # print('spring_force', spring_force)
            # print('locking_force', locking_force)
            # print('sum_moments', sum_moments)
            # print('angular_acceleration', angular_acceleration)
            # print('angular_speed', angular_speed)
            # print('pawl_angle', pawl_angle)
            # print('locking_position', locking_position)


            # Saving results

            wheel_angles.append(wheel_angle)
            locking_forces.append(locking_force)
            locking_positions.append(locking_position)
            pawl_angles.append(pawl_angle)
            time.append(t)


            if pawl_angle < 0:
                # print('finishing simulation as angle<0')
                break

        locking_forces[0] = locking_forces[1]

        return ParkingPawlSimulation(self, time, wheel_angles, pawl_angles,
                                     locking_positions, locking_forces,
                                     name='Locking simulation')

    def max_stress(self, wheel_torque:float):
        return max(self.pawl.max_stress(wheel_torque),
                   self.wheel.max_stress(wheel_torque))

    def size_width(self, wheel_torque:float, safety_factor:float=2):
        min_width = self.minimum_width(wheel_torque, safety_factor=safety_factor)
        self.width = min_width
        self.locking_mechanism.roller_width = min_width
        self.pawl.width = min_width
        self.wheel.width = min_width
        # print('Sizing width: ', min_width)


    def to_markdown(self):
        s = '''## Parking pawl datasheet\n
Engaged slack: {} °\n
Axis/wheel clearance: {} mm\n
ejections levers: {}mm / {}mm\n\n'''.format(round(math.degrees(self.engaged_slack()),3),
                                           round(1000*self.axis_wheel_clearance(), 3),
                                           round(1000*self.ejection_levers()[0], 3),
                                           round(1000*self.ejection_levers()[1], 3)
                                           )

        return s

class ParkingPawlSimulation(dc.DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = ['wheel', 'pawl']

    def __init__(self, parking_pawl:ParkingPawl,
                 time:List[float],
                 wheel_angles:List[float],
                 pawl_angles:List[float],
                 locking_mechanism_positions:List[float],
                 locking_mechanism_forces:List[float],
                 name:str=''):

        max_locking_speed = (parking_pawl.wheel.junction_angle
                             + parking_pawl.wheel.lower_tooth_angle)/time[-1]

        dc.DessiaObject.__init__(self, parking_pawl=parking_pawl,
                                 time=time,
                                 wheel_angles=wheel_angles,
                                 pawl_angles=pawl_angles,
                                 locking_mechanism_positions=locking_mechanism_positions,
                                 locking_mechanism_forces=locking_mechanism_forces,
                                 max_locking_speed=max_locking_speed,
                                 name=name)
        
        self.mass = parking_pawl.mass()
        self.axis_wheel_clearance = parking_pawl.axis_wheel_clearance()
        self.engaged_slack = parking_pawl.engaged_slack()
        self.rest_margin = parking_pawl.rest_margin()
        self.width = parking_pawl.width
        
        
        
    def volmdlr_primitives(self):
        wheel = self.parking_pawl.wheel.volmdlr_primitives()[0]
        pawl = self.parking_pawl.pawl.volmdlr_primitives()[0]
        lock_primitives = self.parking_pawl.locking_mechanism.volmdlr_primitives()

        pawl = pawl.translation(vm.Point3D(0,
                                    -self.parking_pawl.pawl.axis_position.x,
                                    -self.parking_pawl.pawl.axis_position.y))
        return [wheel, pawl] + lock_primitives
    
    def volmdlr_primitives_step_frames(self):
        frames = []
        number_lock_parts = len(self.volmdlr_primitives()) - 2

        for wheel_angle, pawl_angle, locking_position in zip(self.wheel_angles,
                                                             self.pawl_angles,
                                                             self.locking_mechanism_positions):
            wheel_frame = vm.OXYZ.rotation(vm.O3D, vm.X3D, wheel_angle)
            pawl_frame = vm.OXYZ.rotation(vm.O3D, vm.X3D, pawl_angle)
            pawl_frame.origin.y, pawl_frame.origin.z = self.parking_pawl.pawl.axis_position
            locking_frame = vm.Frame3D(vm.Point3D(0, locking_position,
                                                  self.parking_pawl.locking_mechanism_center_distance),
                                       vm.X3D, vm.Y3D, vm.Z3D)
            frames.append([wheel_frame, pawl_frame] + [locking_frame]*number_lock_parts)
            
        return frames

    def plot_data(self):
        # to_disp_attribute_names = ['time', 'locking_force', 'pawl_angle']
        # tooltip = plot_data.Tooltip(to_disp_attribute_names=to_disp_attribute_names)
        elements1 = []
        elements2 = []
        elements3 = []
        for t, locking_force, pawl_angle, locking_position in zip(self.time,
                                                self.locking_mechanism_forces,
                                                self.pawl_angles,
                                                self.locking_mechanism_positions):
            elements1.append({'time': t,
                             'locking_force': locking_force})
            elements2.append({'time': t,
                              'pawl_angle': pawl_angle})
            elements3.append({'time': t,
                              'locking_position': locking_position})

        # The previous line instantiates a dataset with limited arguments but
        # several customizations are available
        point_style = plot_data.PointStyle(color_fill=plot_data.colors.BLUE,
                                            color_stroke=plot_data.colors.BLACK)
        edge_style = plot_data.EdgeStyle(color_stroke=plot_data.colors.BLUE,
                                          )

        dataset1 = plot_data.Dataset(elements=elements1, edge_style=edge_style,
                                     point_style=point_style,
                                    name='locking_force')
        dataset2 = plot_data.Dataset(elements=elements2, edge_style=edge_style,
                                     point_style=point_style,
                                     name='pawl_angle')
        dataset3 = plot_data.Dataset(elements=elements3, edge_style=edge_style,
                                     point_style=point_style,
                                     name='locking_position')
        # print(elements1)
        # print(elements2)
        return [plot_data.Graph2D(graphs=[dataset1], to_disp_attribute_names=['time', 'locking_force']),
                plot_data.Graph2D(graphs=[dataset2], to_disp_attribute_names=['time', 'pawl_angle']),
                plot_data.Graph2D(graphs=[dataset3],
                                  to_disp_attribute_names=['time',
                                                           'locking_position'])
                ]

class ParkingPawSimulationlList(dc.DessiaObject):
    _standalone_in_db = True
    def __init__(self, parking_pawl_simulations:List[ParkingPawlSimulation], name:str=''):
        self.parking_pawl_simulations = parking_pawl_simulations
        self.name = name
        
    def plot_data(self):
        
        dataset = []
        for ipls,  parking_pawl_simulation in enumerate(self.parking_pawl_simulations):
            xmin, xmax, ymin, ymax = parking_pawl_simulation.parking_pawl.functional_footprint()
            width = xmax - xmin
            height = ymax - ymin
            dataset.append({'mass': parking_pawl_simulation.parking_pawl.mass(),
                            'engaged_slack': parking_pawl_simulation.parking_pawl.engaged_slack(),
                            'axis_wheel_clearance': parking_pawl_simulation.parking_pawl.mass(),
                            'rest_margin': parking_pawl_simulation.parking_pawl.rest_margin(),
                            # 'wheel_speed': parking_pawl_simulation.wheel_angles[1]/parking_pawl_simulation.time[1],
                            'width': parking_pawl_simulation.parking_pawl.width,
                            'footprint_width': width,
                            'footprint_height': height
                            })
        
        objects = []
        
        to_disp_attribute_names = ['mass',
                                   'engaged_slack',
                                   'axis_wheel_clearance', 'rest_margin',
                                   'width', 'footprint_width', 'footprint_height']
        
        
        parallel_plot = plot_data.ParallelPlot(disposition='horizontal',
                                               axes=to_disp_attribute_names)
        objects.append(parallel_plot)
        
        # Scatter        
        scatter1 = plot_data.Scatter(elements=dataset,
                                     tooltip = plot_data.Tooltip(to_disp_attribute_names=to_disp_attribute_names),
                                     x_variable='footprint_width', y_variable='mass')
        objects.append(scatter1)
        
        
        scatter2 = plot_data.Scatter(elements=dataset,
                                     tooltip = plot_data.Tooltip(to_disp_attribute_names=to_disp_attribute_names),
                                     x_variable='footprint_width', y_variable='footprint_height')
        objects.append(scatter2)
        
        coords = [(0, 0), (400, 0), (600, 600)]
        sizes = [plot_data.Window(width=860, height=300),
                 plot_data.Window(width=360, height=300),
                 plot_data.Window(width=360, height=300),
                 ]
        
        
        multipleplots = plot_data.MultiplePlots(elements=dataset, plots=objects,
                                                sizes=sizes, coords=coords,
                                                initial_view_on=True)
        
        return [multipleplots]   
        
class ParkingPawlOptimizer(dc_opt.InstantiatingModelOptimizer):
    _standalone_in_db = True
    _allowed_methods = ['optimize_gradient', 'optimize_cma', 'objective_from_model']
    
    def __init__(self, wheel_locking_speed:float,
                 wheel_torque:float,
                 locking_mechanism_travel:float,
                 max_engaged_slack:float,
                 fixed_parameters:List[dc_opt.FixedAttributeValue],
                 optimization_bounds:List[dc_opt.BoundedAttributeValue],
                 name:str=''):
        self.wheel_locking_speed = wheel_locking_speed
        self.wheel_torque = wheel_torque
        self.locking_mechanism_travel = locking_mechanism_travel
        self.max_engaged_slack = max_engaged_slack
        self.fixed_parameters = fixed_parameters
        self.optimization_bounds = optimization_bounds
        self.name = name

        self.number_parameters = len(self.optimization_bounds)
    

    
    
    def instantiate_model(self, attributes_values):
        # print('attributes_values')
        # print(attributes_values)
        wheel_lower_tooth_diameter = attributes_values['wheel_lower_tooth_diameter']
        # basis_diameter = wheel_lower_tooth_diameter + attributes_values['relative_basis_diameter']
        contact_diameter = wheel_lower_tooth_diameter + attributes_values['relative_contact_diameter']
        wheel_outer_diameter = contact_diameter + attributes_values['relative_wheel_outer_diameter']
        finger_height = 1.1*0.5 * (wheel_outer_diameter - wheel_lower_tooth_diameter)
        slope_start_height = finger_height + attributes_values['slope_start_finger_overheight']
        
        locking_mechanism = RollerLockingMechanism(roller_diameter=attributes_values['roller_diameter'],
                                                   roller_width=0.025,
                                                   spring_stiffness=23000,
                                                   spring_active_length=self.locking_mechanism_travel
                                                  )

        parking_pawl = ParkingPawl(wheel_inner_diameter=0.5*wheel_lower_tooth_diameter,
                                   wheel_lower_tooth_diameter=wheel_lower_tooth_diameter,
                                   wheel_outer_diameter=wheel_outer_diameter,
                                   teeth_number=attributes_values['teeth_number'],
                                   lower_tooth_ratio=attributes_values['lower_tooth_ratio'],
                                   pressure_angle=attributes_values['pressure_angle'],
                                   contact_diameter=contact_diameter,
                                   width=0.025,
                                   pawl_offset=attributes_values['pawl_offset'],
                                   axis_inner_diameter=0.5*attributes_values['axis_outer_diameter'],
                                   axis_outer_diameter=attributes_values['axis_outer_diameter'],
                                   finger_height=finger_height,
                                   roller_rest_length=0.6 * locking_mechanism.roller_diameter,
                                   finger_width=attributes_values['finger_width'],
                                   slope_start_height=slope_start_height,
                                   slope_angle=attributes_values['slope_angle'],
                                   slope_offset=0.005,
                                   slope_length=self.locking_mechanism_travel / math.cos(
                                        attributes_values['slope_angle']),
                                   pawl_spring_stiffness=20,
                                   locking_mechanism=locking_mechanism)
        parking_pawl.size_width(self.wheel_torque)
        return parking_pawl

    def dimensionless_vector_to_vector(self, dl_vector):
        return [bound.dimensionless_to_value(dl_xi) for dl_xi, bound in zip(dl_vector, self.optimization_bounds)]

    def vector_to_attributes_values(self, vector:List[float]):
        attributes = {fp.attribute_name: fp.value for fp in self.fixed_parameters}
        
        for bound, xi in zip(self.optimization_bounds, vector):
            attributes[bound.attribute_name] = xi
        return attributes

    def objective_from_dimensionless_vector(self, dl_vector):
        attributes_values = self.vector_to_attributes_values(self.dimensionless_vector_to_vector(dl_vector))

        try:
            model = self.instantiate_model(attributes_values)
        except Exception as err:
            # print('Error in instanciating model')
            return 1000000+100*random.random()
        return self.objective_from_model(model)
        
    def objective_from_model(self, model:ParkingPawl, clearance:float=0.003):

        objective = model.mass()

        if model.axis_wheel_clearance() < clearance:
            objective += (clearance - model.axis_wheel_clearance())*10000

        model_slack = model.engaged_slack()
        if  model_slack < 0.:
            objective -= (model_slack)*100000
        elif model_slack > self.max_engaged_slack:
            objective += (model_slack - self.max_engaged_slack)*100000


        if model.rest_margin() < 0:
            objective -= 100000*model.rest_margin()
            
        lever_limit = 0.01*model.wheel.outer_diameter
        for lever in model.ejection_levers():
            if lever < lever_limit:
                # print('lever', lever, lever_limit)
                objective += (lever_limit - lever)*1000000

        model.pawl.size_torsion_spring(10 * 9.81)
        try:
            simulation = model.locking_simulation(wheel_speed=self.wheel_locking_speed)
        except:
            # print('Something went wrong in simulation')
            objective += 1000000
            return objective
       
        tooth_time = abs((model.wheel.junction_angle+model.wheel.lower_tooth_angle)/self.wheel_locking_speed)
        locking_time_ratio = (tooth_time-simulation.time[-1])/tooth_time
        if locking_time_ratio < 0:
            objective += 10000*abs(locking_time_ratio)
            
        # print('rest margin/locking time ratio', model.rest_margin(), locking_time_ratio)
        # print('\tmodel objective', objective)
        
        return objective

