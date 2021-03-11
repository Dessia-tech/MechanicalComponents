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
    def __init__(self, pitch:float, pin_diameter, roller_outer_diameter, inner_width, overall_width):
        dc.DessiaObject.__init__(self, pitch=pitch,
                                 pin_diameter=pin_diameter,
                                 roller_outer_diameter=roller_outer_diameter,
                                 inner_width=inner_width,
                                 overall_width=overall_width)
        self.bushing_diameter = 0.5*(self.roller_outer_diameter+self.pin_diameter)
        self.plate_width = 0.35*(self.overall_width - self.inner_width)
        self.slack = 0.15*(self.overall_width - self.inner_width)
        self.plate_diameter = 1.2*self.roller_outer_diameter
        
    def plate_outer_contour(self):
        circle1 = vmw.Circle2D(vm.O2D, 0.5*self.plate_diameter)
        circle2 = vmw.Circle2D(self.pitch*vm.X2D, 0.5*self.plate_diameter)
        circle3 = vmw.Circle2D(vm.Point2D(0.5*self.pitch, 2.5*self.plate_diameter), 2.1*self.plate_diameter)
        circle4 = vmw.Circle2D(vm.Point2D(0.5*self.pitch, -2.5*self.plate_diameter), 2.1*self.plate_diameter)
        
        p1 =  sorted(circle1.circle_intersections(circle3), key=lambda p:p.x)[1]
        p2 =  sorted(circle2.circle_intersections(circle3), key=lambda p:p.x)[0]
        p3 =  sorted(circle2.circle_intersections(circle4), key=lambda p:p.x)[0]
        p4 =  sorted(circle1.circle_intersections(circle4), key=lambda p:p.x)[1]
        ax = circle1.plot()
        # circle2.plot(ax=ax)        
        # circle3.plot(ax=ax)        
        # circle4.plot(ax=ax)        
        # p1.plot(ax=ax, color='r')
        # p2.plot(ax=ax, color='g')
        # p3.plot(ax=ax, color='b')
        # p4.plot(ax=ax)
        arc1, _  = circle1.split(p4, p1)
        arc2, a = circle3.split(p1, p2)
        arc2.plot(color='r', ax=ax)
        arc2.interior.plot(ax=ax, color='r')
        a.plot(color='b', ax=ax)
        
        ax.set_aspect('equal')
        _, arc3 = circle2.split(p2, p3)
        arc4, _ = circle4.split(p3, p4)
        
        # arc4.plot(ax=ax, color='g')
        return vmw.Contour2D([arc1, arc2, arc3, arc4])
        
    def plate_inner_contours(self):
        return [vmw.Circle2D(vm.O2D, 0.5*self.pin_diameter),
                vmw.Circle2D(vm.X2D*self.pitch, 0.5*self.pin_diameter)]

        
    def volmdlr_primitives(self, frame=vm.OXYZ):
        pin1 = p3d.Cylinder(frame.origin, frame.u, 0.5*self.pin_diameter, self.overall_width)
        pin2 = p3d.Cylinder(frame.origin+self.pitch*vm.Y3D, frame.u,
                            0.5*self.pin_diameter, self.overall_width)
        outer_plate1 = p3d.ExtrudedProfile(frame.origin, frame.v, frame.w,
                                           self.plate_outer_contour(),
                                           self.plate_inner_contours(),
                                           vm.X3D)
        return [pin1, pin2, outer_plate1]#, outer_plate2, inner_plate1, inner_plate2]

class Sprocket(dc.DessiaObject):
    def _init__(self, pitch:float, number_teeth:int):
        dc.DessiaObject.__init__(self, pitch=pitch, number_teeth=number_teeth)    
    
    
class RollerChainLayout(dc.DessiaObject):
    def __init__(self, roller_chain:RollerChain, sprockets:List[Sprocket]):
        dc.DessiaObject.__init__(self, roller_chain=roller_chain,
                                 sprockets=sprockets)
    
# https://www.tridistribution.fr/chaines-a-rouleaux-norme-europeene-iso/6184-5453-rollerchain-simplex-european-series.html#/4017-ref-04b1/4018-p-6/17601-w-280/17602-o_d-400/17603-o_d-185/17604-c-830
iso_chains = [RollerChain(pitch=0.006, inner_width=0.002,
                          overall_width=0.0083, roller_outer_diameter=0.004,
                          pin_diameter=0.00185)]
        
iso_chains[0].plate_outer_contour().plot()
        