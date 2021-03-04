#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import math
from typing import List
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
                 # upper_tooth_ratio:float,
                 lower_tooth_ratio:float,
                 basis_diameter:float,
                 contact_diameter:float,
                 width:float):
        self.inner_diameter = inner_diameter
        self.lower_tooth_diameter = lower_tooth_diameter
        self.outer_diameter = outer_diameter
        # self.upper_tooth_ratio = upper_tooth_ratio
        self.lower_tooth_ratio = lower_tooth_ratio
        self.teeth_number = teeth_number
        self.basis_diameter = basis_diameter
        self.contact_diameter = contact_diameter
        self.width = width
        
        # Computed attributes
        self.tooth_angle = vm.TWO_PI/self.teeth_number
        # print('tooth_angle', math.degrees(self.tooth_angle))
        # self.upper_tooth_angle = self.tooth_angle*self.upper_tooth_ratio

        # self.junction_angle = 0.5*(self.tooth_angle-self.upper_tooth_angle-self.lower_tooth_angle)

        # self.contact_height_ratio = (self.contact_diameter - self.lower_tooth_diameter)/(self.outer_diameter - self.lower_tooth_diameter)

        self.contact_point1 = vm.Point2D(0, 0.5*self.contact_diameter)

        action_point12 = vm.O2D.rotation(self.contact_point1,
                                          math.asin(self.basis_diameter/self.contact_diameter))

        self.action_line1 = vme.Line2D(self.contact_point1, action_point12)

        self.basis_circle = vmw.Circle2D(vm.O2D, 0.5 * self.basis_diameter)
        self.side_tooth_radius = self.contact_diameter

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

        # print('junction angle: ', math.degrees(self.junction_angle))

        remaining_angle = self.tooth_angle-2*self.junction_angle
        # print('remaining_angle', math.degrees(remaining_angle))
        if remaining_angle < 0:
            raise ValueError('Negative remaining space')

        self.lower_tooth_angle = remaining_angle * self.lower_tooth_ratio
        self.upper_tooth_angle = remaining_angle * (1-self.lower_tooth_ratio)
        # print('lower_tooth_angle', math.degrees(self.lower_tooth_angle))
        # print('upper_tooth_angle', math.degrees(self.upper_tooth_angle))


        self.arc_side1 = vme.Arc2D(corrected_upper_end, self.contact_point1,
                              corrected_lower_start)

        self.contact_point2 = vm.Point2D(0, 0.5 * self.contact_diameter)
        self.contact_point2.rotation(vm.O2D,
                                     (2 * self.lower_junction_angle
                                      + self.lower_tooth_angle),
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
        # a = side1_circle.plot()
        # # side2_circle.plot(ax=a, color='grey')
        # self.action_line1.plot(ax=a, color='k')
        # action_point12.plot(ax=a, color='grey')
        # # self.action_line2.plot(ax=a, color='grey')
        # arc_upper.plot(ax=a, color='r')
        # arc_lower.plot(ax=a, color='b')
        # self.contact_point1.plot(ax=a)
        # self.contact_point2.plot(ax=a)
        # corrected_next_tooth_start.plot(ax=a, color='grey')
        # corrected_upper_start.plot(ax=a, color='g')
        # corrected_upper_end.plot(ax=a, color='r')
        # corrected_lower_start.plot(ax=a, color='r')
        # corrected_lower_end.plot(ax=a, color='r')
        # self.basis_circle.plot(ax=a, color='grey')
        # a.set_aspect('equal')

        arc_upper = arc_upper.split(corrected_upper_start)[1]
        self.arc_upper = arc_upper.split(corrected_upper_end)[0]

        arc_lower = arc_lower.split(corrected_lower_start)[1]
        self.arc_lower = arc_lower.split(corrected_lower_end)[0]

        self.arc_side2 = vme.Arc2D(corrected_lower_end, self.contact_point2,
                                   corrected_next_tooth_start)

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

    def plot_data(self, angle=0.):
        line_style = plot_data.EdgeStyle(color_stroke=plot_data.colors.WATERCRESS,
                                         dashline=[5, 5, 20, 5])

        primitives = [self.outer_contour().rotation(vm.O2D, angle).plot_data(),
                      self.inner_contour().rotation(vm.O2D, angle).plot_data(),
                      self.basis_circle.plot_data(edge_style=line_style),
                      self.action_line1.rotation(vm.O2D, angle).plot_data(edge_style=line_style),
                      self.action_line2.rotation(vm.O2D, angle).plot_data(edge_style=line_style)]
        return [plot_data.PrimitiveGroup(primitives)]
        
    def volmdlr_primitives(self, frame=vm.OXYZ):
        return [p3d.ExtrudedProfile(frame.origin, frame.v, frame.w,
                                    self.outer_contour(), [self.inner_contour()],
                                    frame.u*self.width)]

        
class Pawl(dc.DessiaObject):
    def __init__(self, axis_position:vm.Point2D,
                 wheel_lower_tooth_diameter:float,
                 contact_diameter:float,
                 basis_diameter:float,
                 axis_inner_diameter:float, axis_outer_diameter:float,
                 finger_height:float,
                 # finger_angle:float,
                 finger_width:float,
                 slope_start_height, slope_length,
                 slope_offset:float, slope_angle:float,
                 width:float):
        self.axis_position = axis_position
        self.wheel_lower_tooth_diameter = wheel_lower_tooth_diameter
        self.contact_diameter = contact_diameter
        self.basis_diameter = basis_diameter

        self.finger_height = finger_height
        # self.finger_angle = finger_angle
        self.finger_width = finger_width
        self.slope_start_height = slope_start_height
        self.slope_length = slope_length
        self.slope_offset = slope_offset
        self.slope_angle = slope_angle
        self.axis_inner_diameter = axis_inner_diameter
        self.axis_outer_diameter = axis_outer_diameter
        self.width = width

        self.side_tooth_radius = self.contact_diameter


        pa1 = self.axis_position + vm.Point2D(0., 0.5*self.axis_outer_diameter)
        pa2 = self.axis_position + vm.Point2D(-0.5*self.axis_outer_diameter, 0.)
        pa3 = self.axis_position + vm.Point2D(0., -0.5*self.axis_outer_diameter)
        pa3.rotation(self.axis_position, math.radians(60), copy=False)

        self.axis_arc = vme.Arc2D(pa1, pa2, pa3)

        contact_middle_point = vm.Point2D(0, 0.5*self.contact_diameter)
        self.contact_angle = math.asin(self.finger_width/self.contact_diameter)
        self.contact_point1 = contact_middle_point.rotation(vm.O2D, -0.5*self.contact_angle)
        self.contact_point2 = contact_middle_point.rotation(vm.O2D, 0.5 * self.contact_angle)
        action_point12 = vm.O2D.rotation(self.contact_point1,
                                          math.asin(self.basis_diameter/self.contact_diameter))
        action_point22 = vm.O2D.rotation(self.contact_point2,
                                          -math.asin(self.basis_diameter/self.contact_diameter))

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

        side1_start  = sorted(side1_circle.line_intersections(lower_finger_line), key=lambda p:p.x)[1]
        side1_end  = sorted(side1_circle.line_intersections(upper_finger_line), key=lambda p:p.x)[1]

        side2_end  = sorted(side2_circle.line_intersections(lower_finger_line), key=lambda p:p.x)[0]
        side2_start  = sorted(side2_circle.line_intersections(upper_finger_line), key=lambda p:p.x)[0]

        # ax = side1_circle.plot(color='grey')
        # side2_circle.plot(ax=ax, color='grey')
        # self.action_line1.plot(ax=ax, color='grey')
        # self.action_line2.plot(ax=ax, color='grey')
        # self.contact_point1.plot(ax=ax)
        # self.contact_point2.plot(ax=ax)
        # side1_start.plot(ax=ax, color='r')
        # side1_end.plot(ax=ax, color='g')
        # side2_start.plot(ax=ax, color='g')
        # side2_end.plot(ax=ax, color='r')

        self.arc_side1 = vme.Arc2D(side1_start, self.contact_point1, side1_end)
        self.arc_side2 = vme.Arc2D(side2_start, self.contact_point2, side2_end)

        self.lower_finger_line = vme.LineSegment2D(side2_end, side1_start)

        self.junction1 = vme.LineSegment2D(pa3, side2_start)

        finger_lower_end = side1_end + vm.Point2D(self.slope_offset, 0.)
        slope_start = vm.Point2D(finger_lower_end.x, self.slope_start_height+0.5*self.wheel_lower_tooth_diameter)
        self.junction2 = vme.LineSegment2D(side1_end, finger_lower_end)
        self.junction3 = vme.LineSegment2D(finger_lower_end, slope_start)

        slope_stop = slope_start - vm.Point2D(self.slope_length, 0.).rotation(vm.O2D, -self.slope_angle)
        self.slope = vme.LineSegment2D(slope_start, slope_stop)
        self.junction4 = vme.LineSegment2D(slope_stop, pa1)

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
        

        
        primitives = [self.axis_arc, self.junction1, self.arc_side2,
                      self.lower_finger_line, self.arc_side1, self.junction2,
                      self.junction3, self.slope, self.junction4]
        
        return vmw.Contour2D(primitives)

    def inner_contour(self):
        return vmw.Circle2D(self.axis_position, 0.5*self.axis_inner_diameter)

    def plot_data(self, angle=0.):
        line_style = plot_data.EdgeStyle(color_stroke=plot_data.colors.FERN,
                                         dashline=[5, 5, 20, 5])
        primitives = [self.outer_contour().rotation(self.axis_position, angle).plot_data(),
                     self.inner_contour().rotation(self.axis_position, angle).plot_data(),
                     self.action_line1.rotation(self.axis_position, angle).plot_data(edge_style=line_style),
                     self.action_line2.rotation(self.axis_position, angle).plot_data(edge_style=line_style)]


        return [plot_data.PrimitiveGroup(primitives)]

    def volmdlr_primitives(self, frame=vm.OXYZ):
        return [p3d.ExtrudedProfile(frame.origin, frame.v, frame.w,
                                    self.outer_contour(), [self.inner_contour()],
                                    frame.u*self.width)]

class RollerLockingMechanism(dc.DessiaObject):
    _standalone_in_db = True

    def __init__(self, roller_diameter:float,
                 # center_distance:float,
                 # start_position:float,
                 # end_position:float,
                 width:float,
                 name:str=''):
        # self.start_position = start_position
        # self.end_position = end_position
        # self.center_distance = center_distance
        self.roller_diameter = roller_diameter
        self.width = width
        self.name = name

    def outer_contour(self):
        return vmw.Circle2D(vm.Point2D(0., 0.),
                            0.5 * self.roller_diameter)
    def inner_contour(self):
        return vmw.Circle2D(vm.Point2D(0., 0.),
                            0.25 * self.roller_diameter)

    def plot_data(self):
        # line_style = plot_data.EdgeStyle(color_stroke=plot_data.colors.FERN,
        #                                  dashline=[5, 5, 20, 5])
        primitives = [self.outer_contour().plot_data(),
                     self.inner_contour().plot_data(),
                     ]

        return [plot_data.PrimitiveGroup(primitives)]

    def volmdlr_primitives(self, frame=vm.OXYZ):
        # roller = p3d.ExtrudedProfile(frame.origin, frame.v, frame.w,
        #                             self.outer_contour(), [],
        #                             frame.u*self.width)
        roller = p3d.HollowCylinder(position=0.5*self.width*vm.X3D,
                                    axis=vm.X3D,
                                    inner_radius=0.25*self.roller_diameter,
                                    outer_radius=0.5*self.roller_diameter,
                                    length=self.width)
        # roller = p3d.HollowCylinder(0.1 * vm.Y3D, vm.X3D,
        #                                         0.02, 0.06, 0.03,
        #                                         name='cylinder2')
        return [roller]


class ParkingPawl(dc.DessiaObject):
    _standalone_in_db = True
    def __init__(self,
                 wheel_inner_diameter:float,
                 wheel_lower_tooth_diameter:float,
                 wheel_outer_diameter:float,
                 teeth_number:int,
                 lower_tooth_ratio:float,
                 basis_diameter:float,
                 contact_diameter:float,
                 width:float,
                 pawl_offset:float,
                 axis_inner_diameter:float, axis_outer_diameter:float,
                 finger_height:float,
                 finger_width:float,
                 slope_start_height, slope_length,
                 slope_offset:float, slope_angle:float,
                 locking_mechanism:RollerLockingMechanism,
                 open_clearance:float=0.002,
                 name:str=''
                 ):
        self.locking_mechanism = locking_mechanism
        self.open_clearance = open_clearance
        self.name = name

        self.wheel = Wheel(inner_diameter=wheel_inner_diameter,
                           lower_tooth_diameter=wheel_lower_tooth_diameter,
                           outer_diameter=wheel_outer_diameter,
                           teeth_number=teeth_number,
                           # upper_tooth_ratio=upper_tooth_ratio,
                           lower_tooth_ratio=lower_tooth_ratio,
                           basis_diameter=basis_diameter,
                           contact_diameter=contact_diameter,
                           width=width)
        self.slope_end_height = 0.5*wheel_lower_tooth_diameter + slope_start_height + slope_length*math.sin(slope_angle)

        self.pawl = Pawl(axis_position=vm.Point2D(-pawl_offset, self.slope_end_height-0.5*axis_outer_diameter),
                         wheel_lower_tooth_diameter=wheel_lower_tooth_diameter,
                         contact_diameter=contact_diameter,
                         basis_diameter=basis_diameter,
                         axis_inner_diameter=axis_inner_diameter,
                         axis_outer_diameter=axis_outer_diameter,
                         finger_height=finger_height,
                         # finger_angle=finger_angle,
                         finger_width=finger_width,
                         slope_start_height=slope_start_height,
                         slope_length=slope_length,
                         slope_offset=slope_offset,
                         slope_angle=slope_angle,
                         width=width
                         )
        


        self.contact1_wheel_angle = -0.5*self.pawl.contact_angle
        self.contact2_wheel_angle = -(self.wheel.lower_tooth_angle
                                     +2*self.wheel.lower_junction_angle
                                     -0.5*self.pawl.contact_angle)

        R = self.pawl.arc_side2.end.point_distance(self.pawl.axis_position)
        x, y = self.pawl.arc_side2.end-self.pawl.axis_position
        alpha = math.atan2(y, x)
        self.up_pawl_angle = -(math.asin(math.sin(alpha)
                               -(0.5*(wheel_outer_diameter
                                      -wheel_lower_tooth_diameter)
                                 +open_clearance)/R)
                               - alpha)

        # Finding locking mech center distance
        y, z = self.pawl.slope.start.rotation(self.pawl.axis_position,
                                              self.up_pawl_angle)
        self.locking_mechanism_center_distance = z + 0.5*self.locking_mechanism.roller_diameter
        self.locking_mechanism_start_position = y
        
        self.locking_mechanism_end_position = self.pawl.slope.end.x
        

    def volmdlr_primitives(self, frame=vm.OXYZ):
        wheel_frame = frame.rotation(frame.origin, frame.u, self.contact1_wheel_angle)
        locking_mech_frame = vm.OXYZ.copy()
        locking_mech_frame.origin.z = self.locking_mechanism_center_distance
        return (
                self.wheel.volmdlr_primitives(frame=wheel_frame)
                + self.pawl.volmdlr_primitives()
                + self.locking_mechanism.volmdlr_primitives(frame=locking_mech_frame))

    def engaged_slack(self):
        return self.wheel.lower_tooth_angle * 0.5 * self.wheel.lower_tooth_diameter - self.pawl.finger_width

    def check(self):
        if self.engaged_slack() < 0:
            return False
        if self.pawl.slope_angle < self.up_pawl_angle:
            return False
        return True

    def plot_data(self, pawl_angle=None, wheel_angle=None):

        if pawl_angle is None and wheel_angle is None:
            primitives_p1 = self.locking_mechanism.plot_data()
            primitives_p1 += self.pawl.plot_data(angle=0.)[0].primitives
            primitives_p1 += self.wheel.plot_data(angle=self.contact1_wheel_angle)[0].primitives

            primitives_p2 = self.locking_mechanism.plot_data()
            primitives_p2 += self.pawl.plot_data(angle=0.)[0].primitives
            primitives_p2 += self.wheel.plot_data(angle=self.contact2_wheel_angle)[0].primitives

            primitives_p3 = self.locking_mechanism.plot_data()
            primitives_p3 += self.pawl.plot_data(angle=self.up_pawl_angle)[0].primitives
            primitives_p3 += self.wheel.plot_data(angle=0.)[0].primitives

            return [plot_data.PrimitiveGroup(primitives_p1),
                    plot_data.PrimitiveGroup(primitives_p2),
                    plot_data.PrimitiveGroup(primitives_p3)]
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
        
    def static_locking_simulation(self, distance_step=0.002):
        wheel_angles = []
        pawl_angles = []
        locking_positions = []
        
        travel = self.locking_mechanism_end_position - self.locking_mechanism_start_position
        number_steps = round(abs(travel)/distance_step) + 1
        distance_step = travel/distance_step
        print(travel)
        
        for i in range(number_steps+1):
            locking_position = self.locking_mechanism_start_position + i/(number_steps)*travel
            wheel_angles.append(0.)
            locking_positions.append(locking_position)
            pawl_angles.append(self.up_pawl_angle - i/(number_steps)*self.up_pawl_angle)
        
        return ParkingPawlSimulation(self, wheel_angles, pawl_angles, locking_positions,
                                     name='Locking simulation')

class ParkingPawlSimulation(dc.DessiaObject):
    _standalone_in_db = True

    def __init__(self, parking_pawl:ParkingPawl,
                 wheel_angles:List[float],
                 pawl_angles:List[float],
                 locking_mechanism_positions:List[float],
                 name:str=''):
        dc.DessiaObject.__init__(self, parking_pawl=parking_pawl,
                                 wheel_angles=wheel_angles,
                                 pawl_angles=pawl_angles,
                                 locking_mechanism_positions=locking_mechanism_positions,
                                 name=name)
        
    def volmdlr_primitives(self):
        wheel, pawl, lock = self.parking_pawl.volmdlr_primitives()
    
        pawl = pawl.translation(vm.Point3D(0,
                                    -self.parking_pawl.pawl.axis_position.x,
                                    -self.parking_pawl.pawl.axis_position.y))
        return [wheel, pawl, lock]
    
    def volmdlr_primitives_step_frames(self):
        frames = []
        for wheel_angle, pawl_angle, locking_position in zip(self.wheel_angles,
                                                             self.pawl_angles,
                                                             self.locking_mechanism_positions):
            wheel_frame = vm.OXYZ.rotation(vm.O3D, vm.X3D, wheel_angle)
            pawl_frame = vm.OXYZ.rotation(vm.O3D, vm.X3D, pawl_angle)
            pawl_frame.origin.y, pawl_frame.origin.z = self.parking_pawl.pawl.axis_position
            locking_frame = vm.Frame3D(vm.Point3D(0, locking_position,
                                                  self.parking_pawl.locking_mechanism_center_distance),
                                       vm.X3D, vm.Y3D, vm.Z3D)
            frames.append([wheel_frame, pawl_frame, locking_frame])
            
        return frames
