#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

from dessia_common.core import DessiaObject

import volmdlr as vm
import volmdlr.primitives2d as primitives2d
import volmdlr.primitives3d as primitives3d
#import matplotlib.pyplot as plt

class Section(DessiaObject):
    def __init__(self):
        pass
    
    
class ISection(Section):
    def __init__(self, h, b, tw, tf, r):
        self.h = h
        self.b = b
        self.tw = tw# Thickness web
        self.tf = tf# Thickness flange
        self.r = r
        
    def contour(self, x=vm.X2D, z=vm.Y2D):
        p1 = vm.Point2D(-0.5*self.b, -0.5*self.h)
        p2 = vm.Point2D(0.5*self.b, -0.5*self.h)
        p3 = vm.Point2D(0.5*self.b, -0.5*self.h+self.tf)
        p4 = vm.Point2D(0.5*self.tw, -0.5*self.h + self.tf)
        p5 = vm.Point2D(0.5*self.tw, 0.5*self.h - self.tf)
        p6 = vm.Point2D(0.5*self.b, 0.5*self.h - self.tf)
        p7 = vm.Point2D(0.5*self.b, 0.5*self.h)
        p8 = vm.Point2D(-0.5*self.b, 0.5*self.h)
        p9 = vm.Point2D(-0.5*self.b, 0.5*self.h-self.tf)
        p10 = vm.Point2D(-0.5*self.tw, 0.5*self.h-self.tf)
        p11 = vm.Point2D(-0.5*self.tw, -0.5*self.h+self.tf)
        p12 = vm.Point2D(-0.5*self.b, -0.5*self.h+self.tf)

        rl = primitives2d.ClosedRoundedLineSegments2D([p1, p2, p3, p4, p5, p6, p7, p8,
                                                       p9, p10, p11, p12],
                                                      radius={3: self.r, 4: self.r,
                                                              9: self.r, 10: self.r})
        return rl
    
    def plot_data(self):
        return [self.contour().plot_data()]
        



class Beam(DessiaObject):
    def __init__(self, section, length, name=''):
        self.section = section
        self.length = length
        self.name = name
        
    def volmdlr_primitives(self, position=vm.O3D, x=vm.X3D, y=vm.Y3D):
        """
        x and y define the plane of section in the beam axis
        """
        z = x.cross(y)
        c = self.section.contour()
        return [primitives3d.ExtrudedProfile(position-0.5*self.length*z, x, y, c, [], z*self.length)]
        
