#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import math
from typing import List
import dessia_common as dc
import volmdlr as vm
import volmdlr.edges as vme
import volmdlr.wires as vmw
import volmdlr.primitives3d as p3d


class RollerChain(dc.DessiaObject):
    _standalone_in_db = True
    
    # slack_plate_ratio = 5
    def __init__(self, pitch:float, roller_diameter:float, roller_width:float,
                 inner_plate_width:float, outer_plate_width:float,
                 pin_length:float, inner_plate_height:float, 
                 outer_plate_height:float, pin_diameter:float,
                 number_rows:int=1,
                 name:str=''):
        dc.DessiaObject.__init__(self, pitch=pitch,
                                 pin_diameter=pin_diameter,
                                 roller_diameter=roller_diameter,
                                 roller_width=roller_width,
                                 outer_plate_width=outer_plate_width,
                                 inner_plate_width=inner_plate_width,
                                 pin_length=pin_length,
                                 inner_plate_height=inner_plate_height,
                                 outer_plate_height=outer_plate_height,
                                 number_rows=number_rows,
                                 name=name)
        self.bushing_diameter = 0.5*(self.roller_diameter+self.pin_diameter)
        self.slack = 0.5*(0.5*(self.pin_length - self.roller_width) - self.outer_plate_width-self.inner_plate_width)
        # self.plate_width = 0.5*(self.outer_plates_width - self.inner_plates_width)*self.slack_plate_ratio/(2*self.slack_plate_ratio+1)
        # self.slack = self.plate_width/self.slack_plate_ratio
        # self.roller_width = self.width
        # self.plate_diameter = 1.2*self.roller_diameter
        # self.outer_plates_distance = self.overall_width - 2*self.slack
        # self.inner_plates_distance = self.outer_plates_distance-self.plate_width
        
    def plate_outer_contour(self, plate_height):
        circle1 = vmw.Circle2D(vm.O2D, 0.5*plate_height)
        circle2 = vmw.Circle2D(self.pitch*vm.X2D, 0.5*plate_height)

        center3 = vm.Point2D(0.5*self.pitch, 2.5*plate_height)
        line1 = vme.Line2D(circle1.center, center3)
        p1 = sorted(circle1.line_intersections(line1), key=lambda p:p.x)[1]
        circle3 = vmw.Circle2D(center3, center3.point_distance(p1))

        center4 = vm.Point2D(0.5*self.pitch, -2.5*plate_height)

        circle4 = vmw.Circle2D(center4, circle3.radius)
        
        p2 = vm.Point2D(self.pitch-p1.x, p1.y)
        p3 = vm.Point2D(p2.x, -p1.y)
        p4 = vm.Point2D(p1.x, p3.y)

        _, arc1  = circle1.split(p4, p1)
        arc2 ,_  = circle3.split(p1, p2)

        arc3 , _ = circle2.split(p2, p3)
        arc4, _ = circle4.split(p3, p4)
        
        return vmw.Contour2D([arc1, arc2, arc3, arc4])
        
    def plate_inner_contours(self):
        return [vmw.Circle2D(vm.O2D, 0.5*self.pin_diameter),
                vmw.Circle2D(vm.X2D*self.pitch, 0.5*self.pin_diameter)]

        
    def volmdlr_primitives(self, points=None, frame=vm.OXYZ):
        if points is None:
            points = [vm.Point2D(0., 0.), vm.Point2D(self.pitch, 0.),
                      vm.Point2D(2*self.pitch, 0.)]
            
        
        primitives = []
        outer_plate_outer_contour = self.plate_outer_contour(self.outer_plate_height)
        inner_plate_outer_contour = self.plate_outer_contour(self.inner_plate_height)
        plate_inner_contours = self.plate_inner_contours()
        for point1_2d, point2_2d, point3_2d in zip(points[::2],
                                                   (points[1:]+[points[0]])[::2],
                                                   (points[2:]+points[:2])[::2]):
            point1_3d = point1_2d.to_3d(frame.origin, frame.u, frame.v)
            point2_3d = point2_2d.to_3d(frame.origin, frame.u, frame.v)
            point3_3d = point3_2d.to_3d(frame.origin, frame.u, frame.v)
            u1 = (point2_3d - point1_3d)
            u2 = (point3_3d - point2_3d)
            u1.normalize()
            u2.normalize()
            v1 = frame.w.cross(u1)
            v2 = frame.w.cross(u2)
            
            pin1 = p3d.Cylinder(point1_3d, frame.w, 0.5*self.pin_diameter,
                                self.pin_length,
                                name='pin 1')
            pin2 = p3d.Cylinder(point2_3d, frame.w,
                                0.5*self.pin_diameter, self.pin_length,
                                name='pin 2')
            delta_w = -((self.number_rows)/2)*(self.pin_length-2*self.slack)*frame.w
            primitives.extend([pin1, pin2])
            for row_index in range(self.number_rows):
                point1_row = point1_3d + row_index*delta_w
                point2_row = point2_3d + row_index*delta_w
                
                roller1 = p3d.HollowCylinder(point1_row, frame.w,
                                              0.5*self.bushing_diameter,
                                              0.5*self.roller_diameter,
                                              self.roller_width, name='roller 1')
                roller2 = p3d.HollowCylinder(point2_row, frame.w,
                                              0.5*self.bushing_diameter,
                                              0.5*self.roller_diameter,
                                              self.roller_width, name='roller 2')
                outer_plate1 = p3d.ExtrudedProfile(point1_row+(0.5*self.pin_length-self.slack)*frame.w,
                                                    u1, v1,
                                                    outer_plate_outer_contour,
                                                    plate_inner_contours,
                                                    -self.outer_plate_width*frame.w,
                                                    name='outer plate 1')
                outer_plate2 = p3d.ExtrudedProfile(point1_row-(0.5*self.pin_length-self.slack)*frame.w,
                                                    u1, v1,
                                                    outer_plate_outer_contour,
                                                    plate_inner_contours,
                                                    self.outer_plate_width*frame.w,
                                                    name='outer plate 2')
                inner_plate1_position = (point2_row
                                          + (0.5 * (self.roller_width)) * frame.w
                                          )
                inner_plate1 = p3d.ExtrudedProfile(inner_plate1_position, u2, v2,
                                                    inner_plate_outer_contour,
                                                    plate_inner_contours,
                                                    self.inner_plate_width*frame.w,
                                                    name='inner plate 1')
                inner_plate2_position = (point2_row
                                          - (0.5 * (self.roller_width)) * frame.w
                                          )
        
                inner_plate2 = p3d.ExtrudedProfile(inner_plate2_position, u2, v2,
                                                    inner_plate_outer_contour,
                                                    plate_inner_contours,
                                                    -self.inner_plate_width*frame.w,
                                                    name='inner plate 2')
                primitives += [outer_plate1, outer_plate2, inner_plate1, inner_plate2,
                               roller1, roller2]
            # primitives += [pin1, pin2, roller1, ro, outer_plate1]
        return primitives

class Sprocket(dc.DessiaObject):
    def __init__(self, pitch:float, number_teeth:int, roller_diameter, width:float, name:str=''):
        diameter = number_teeth*pitch/math.pi
        dc.DessiaObject.__init__(self,
                                 pitch=pitch,
                                 roller_diameter=roller_diameter,
                                 number_teeth=number_teeth,
                                 width=width,
                                 diameter=diameter,
                                 name=name)

    def outer_contour(self):
        c = vm.Point2D(0, 0.5*self.diameter)
        i = vm.Point2D(0, 0.5*(self.diameter - self.roller_diameter))
        s = i.rotation(c, 0.5*math.pi)
        e = i.rotation(c, -0.5*math.pi)
        arc = vme.Arc2D(s, i, e)
        pattern_angle = 2*math.pi/self.number_teeth
        roller_angle = math.atan(2*self.roller_diameter/self.diameter)
        teeth_angle = pattern_angle -roller_angle 
        crest_angle =  0.5*teeth_angle
        junction_angle = 0.5*(teeth_angle - crest_angle)
        crest_start = vm.Point2D(0, 0.5*(self.diameter + self.roller_diameter)).rotation(vm.O2D, 0.5*(roller_angle)+junction_angle)
        crest_interior = crest_start.rotation(vm.O2D, 0.5*crest_angle)
        crest_end = crest_start.rotation(vm.O2D, crest_angle)
        crest = vme.Arc2D(crest_start, crest_interior, crest_end)
        
        junction1 = vme.LineSegment2D(arc.end, crest.start)
        junction2 = vme.LineSegment2D(crest.end, s.rotation(vm.O2D, pattern_angle))
        edges = [arc, junction1, crest, junction2]
        for i in range(1, self.number_teeth):
            edges.append(arc.rotation(vm.O2D, i*pattern_angle))
            edges.append(junction1.rotation(vm.O2D, i*pattern_angle))
            edges.append(crest.rotation(vm.O2D, i*pattern_angle))
            edges.append(junction2.rotation(vm.O2D, i*pattern_angle))
            
        return vmw.Contour2D(edges)
    
    def inner_contour(self):
        return vmw.Circle2D(vm.O2D, 0.5*self.diameter-2*self.roller_diameter)

    def volmdlr_primitives(self, frame=vm.OXYZ):
        sprocket = p3d.Cylinder(frame.origin, frame.w,
                                0.5*self.pitch*self.number_teeth,
                                self.width)
        sprocket = p3d.ExtrudedProfile(frame.origin-0.5*self.width*frame.w, frame.u, frame.v,
                                       self.outer_contour(), [self.inner_contour()],
                                       self.width*frame.w)
        return [sprocket]#, outer_plate2, inner_plate1, inner_plate2]
    
    
class SprocketLayout(dc.DessiaObject):
    def __init__(self, sprocket:Sprocket, position: vm.Point2D,
                 winding_side:bool, name:str=''):
        dc.DessiaObject.__init__(self, sprocket=sprocket,
                                 position=position, winding_side=winding_side,
                                 name=name)
        
        
        
    def volmdlr_primitives(self, frame=vm.OXYZ):
        sprocket_origin = self.position.to_3d(frame.origin, frame.u, frame.v)
            
        frame_s = frame.copy()
        frame_s.origin = sprocket_origin
                         
        return self.sprocket.volmdlr_primitives(frame=frame_s)
        
    
class RollerChainLayout(dc.DessiaObject):
    def __init__(self, roller_chain:RollerChain,
                 sprocket_layouts:List[SprocketLayout],
                 frame=vm.OXYZ,
                 name:str=''):
        dc.DessiaObject.__init__(self, roller_chain=roller_chain,
                                 sprocket_layouts=sprocket_layouts,
                                 frame=frame, 
                                 name=name)
        
        
    def chain_contour(self):
        circles = []
        # ax=None
        circle_to_sprocket = {}
        for sprocket_layout in self.sprocket_layouts:
            circle = vmw.Circle2D(sprocket_layout.position, 0.5*sprocket_layout.sprocket.diameter)
            circles.append(circle)
            circle_to_sprocket[circle] = sprocket_layout
            # ax = circle.plot(ax=ax, color='grey')
            
        line_segments = []
        arcs = []
        circle_ends = {}
        circle_starts = {}
        
        for circle1, sprocket_layout1, circle2, sprocket_layout2 in\
            zip(circles, self.sprocket_layouts,
                circles[1:]+[circles[0]], self.sprocket_layouts[1:]+[self.sprocket_layouts[0]]):
            
            r, R = sorted([circle1.radius, circle2.radius])
            L = circle1.center.point_distance(circle2.center)
            if sprocket_layout1.winding_side == sprocket_layout2.winding_side:
                theta = math.asin((R-r)/L)
            else:
                theta = math.asin((R+r)/L)
            
            u = circle2.center - circle1.center
            u.normalize()
            v = u.normal_vector()

            if sprocket_layout1.winding_side:
                p1b = (circle1.center-circle1.radius*v)
                winding_side1 = 1
            else:
                p1b = (circle1.center+circle1.radius*v)
                winding_side1 = -1
            if sprocket_layout2.winding_side:
                winding_side2 = 1
                p2b = (circle2.center-circle2.radius*v)
            else:
                winding_side2 = -1
                p2b = (circle2.center+circle2.radius*v)
            
            
            if sprocket_layout1.winding_side == sprocket_layout2.winding_side:
                if circle1.radius < circle2.radius:
                    # Transition to small circle to big one
                    angle_side1 = -winding_side1
                    angle_side2 = -winding_side2
                else:
                    angle_side1 = winding_side1
                    angle_side2 = winding_side2
            else:
                angle_side1 = winding_side1
                angle_side2 = -winding_side2

            
            p1 = p1b.rotation(circle1.center, angle_side1*theta)
            p2 = p2b.rotation(circle2.center, angle_side2*theta)
            
            
            circle_starts[circle2] = p2
            circle_ends[circle1] = p1
            # p2b.plot(ax=ax, color='grey')
            # p1b.plot(ax=ax, color='grey')
            
            # p2.plot(ax=ax, color='r')
            # p1.plot(ax=ax, color='g')
            
            ls = vme.LineSegment2D(p1, p2)
            # ls.plot(ax=ax)
            line_segments.append(ls)
            
        for circle, sprocket_layout in zip(circles, self.sprocket_layouts):
            arc1, arc2 = circle.split(circle_starts[circle], circle_ends[circle])
            if sprocket_layout.winding_side:
                if arc1.is_trigo:
                    arcs.append(arc1)
                else:
                    arcs.append(arc2)
            else:
                if arc2.is_trigo:
                    arcs.append(arc1)
                else:
                    arcs.append(arc2)

        edges = []
        for arc, ls in zip(arcs, line_segments):
            edges.append(arc)
            edges.append(ls)
            
        return vmw.Contour2D(edges)
            
                    
    def chain_points(self):
        c = self.chain_contour()
        l = c.length()
        n_points = int(l/self.roller_chain.pitch)
        return [c.point_at_abscissa(i*self.roller_chain.pitch) for i in range(n_points)]
    
    def volmdlr_primitives(self):
        primitives = []
        for sprocket_layout in self.sprocket_layouts:
            # sprocket_origin = sprocket_layout.position.to_3d(self.frame.origin, self.frame.u, self.frame.v)
            
            primitives.extend(sprocket_layout.volmdlr_primitives(vm.Frame3D(vm.O3D,
                                                                            self.frame.u,
                                                                            self.frame.v,
                                                                            self.frame.w)))
        points = self.chain_points()
        primitives += self.roller_chain.volmdlr_primitives(points)
        return primitives