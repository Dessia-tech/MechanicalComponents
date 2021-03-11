#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

from typing import List
import dessia_common as dc
import volmdlr as vm
import volmdlr.edges as vme
import volmdlr.wires as vmw
import volmdlr.primitives3d as p3d


class RollerChain(dc.DessiaObject):
    slack_plate_ratio = 5
    def __init__(self, pitch:float, pin_diameter, roller_outer_diameter, inner_width, overall_width):
        dc.DessiaObject.__init__(self, pitch=pitch,
                                 pin_diameter=pin_diameter,
                                 roller_outer_diameter=roller_outer_diameter,
                                 inner_width=inner_width,
                                 overall_width=overall_width)
        self.bushing_diameter = 0.5*(self.roller_outer_diameter+self.pin_diameter)
        self.plate_width = 0.5*(self.overall_width - self.inner_width)*self.slack_plate_ratio/(2*self.slack_plate_ratio+1)
        self.slack = self.plate_width/self.slack_plate_ratio
        self.plate_diameter = 1.2*self.roller_outer_diameter
        self.outer_plates_distance = self.overall_width - 2*self.slack
        self.inner_plates_distance = self.outer_plates_distance-self.plate_width
        
    def plate_outer_contour(self):
        circle1 = vmw.Circle2D(vm.O2D, 0.5*self.plate_diameter)
        circle2 = vmw.Circle2D(self.pitch*vm.X2D, 0.5*self.plate_diameter)

        center3 = vm.Point2D(0.5*self.pitch, 2.5*self.plate_diameter)
        line1 = vme.Line2D(circle1.center, center3)
        p1 = sorted(circle1.line_intersections(line1), key=lambda p:p.x)[1]
        circle3 = vmw.Circle2D(center3, center3.point_distance(p1))

        center4 = vm.Point2D(0.5*self.pitch, -2.5*self.plate_diameter)

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

        
    def volmdlr_primitives(self, frame=vm.OXYZ):
        pin1 = p3d.Cylinder(frame.origin, frame.u, 0.5*self.pin_diameter,
                            self.overall_width,
                            name='pin 1')
        pin2 = p3d.Cylinder(frame.origin+self.pitch*vm.Y3D, frame.u,
                            0.5*self.pin_diameter, self.overall_width,
                            name='pin 2')
        roller1 = p3d.HollowCylinder(frame.origin, frame.u,
                                      0.5*self.bushing_diameter,
                                      0.5*self.roller_outer_diameter,
                                      self.inner_width, name='bushing 1')
        roller2 = p3d.HollowCylinder(frame.origin+self.pitch*vm.Y3D, frame.u,
                                      0.5*self.bushing_diameter,
                                      0.5*self.roller_outer_diameter,
                                      self.inner_width, name='bushing 2')
        outer_plate1 = p3d.ExtrudedProfile(frame.origin+(0.5*self.overall_width-self.slack)*frame.u,
                                           frame.v, frame.w,
                                           self.plate_outer_contour(),
                                           self.plate_inner_contours(),
                                           -self.plate_width*frame.u,
                                           name='outer plate 1')
        outer_plate2 = p3d.ExtrudedProfile(frame.origin-(0.5*self.overall_width-self.slack)*frame.u,
                                           frame.v, frame.w,
                                           self.plate_outer_contour(),
                                           self.plate_inner_contours(),
                                           self.plate_width*frame.u,
                                           name='outer plate 2')
        inner_plate1_position = (frame.origin
                                 + (0.5 * self.overall_width
                                    - self.plate_diameter
                                    - self.slack) * frame.u
                                 + self.pitch*frame.v)
        inner_plate1 = p3d.ExtrudedProfile(inner_plate1_position, frame.v, frame.w,
                                           self.plate_outer_contour(),
                                           self.plate_inner_contours(),
                                           -self.plate_width*frame.u,
                                           name='inner plate 1')
        inner_plate2_position = (frame.origin
                                 - (0.5 * self.overall_width
                                    - self.plate_diameter
                                    - self.slack) * frame.u
                                 + self.pitch*frame.v)

        inner_plate2 = p3d.ExtrudedProfile(inner_plate2_position, frame.v, frame.w,
                                           self.plate_outer_contour(),
                                           self.plate_inner_contours(),
                                           self.plate_width*frame.u,
                                           name='inner plate 2')
        return [outer_plate1, outer_plate2, inner_plate1, inner_plate2,
                pin1, pin2, roller1, roller2]

class Sprocket(dc.DessiaObject):
    def _init__(self, pitch:float, number_teeth:int, width:float, name:str=''):
        dc.DessiaObject.__init__(self, pitch=pitch,
                                 number_teeth=number_teeth,
                                 width=width,
                                 name=name)

    def volmdlr_primitives(self, frame=vm.OXYZ):
        sprocket = p3d.Cylinder(frame.origin, frame.u,
                                0.5*self.pitch*self.number_teeth,
                                self.overall_width)
        return [sprocket]#, outer_plate2, inner_plate1, inner_plate2]
    
class RollerChainLayout(dc.DessiaObject):
    def __init__(self, roller_chain:RollerChain, sprockets:List[Sprocket]):
        dc.DessiaObject.__init__(self, roller_chain=roller_chain,
                                 sprockets=sprockets)
    


