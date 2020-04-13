#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

from dessia_common.core import DessiaObject

import volmdlr as vm
import volmdlr.primitives2D as primitives2D
import volmdlr.primitives3D as primitives3D
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
        p1 = vm.Point2D((-0.5*self.b, -0.5*self.h))
        p2 = vm.Point2D((0.5*self.b, -0.5*self.h))
        p3 = vm.Point2D((0.5*self.b, -0.5*self.h+self.tf))
        p4 = vm.Point2D((0.5*self.tw, -0.5*self.h + self.tf))
        p5 = vm.Point2D((0.5*self.tw, 0.5*self.h - self.tf))
        p6 = vm.Point2D((0.5*self.b, 0.5*self.h - self.tf))
        p7 = vm.Point2D((0.5*self.b, 0.5*self.h))
        p8 = vm.Point2D((-0.5*self.b, 0.5*self.h))
        p9 = vm.Point2D((-0.5*self.b, 0.5*self.h-self.tf))
        p10 = vm.Point2D((-0.5*self.tw, 0.5*self.h-self.tf))
        p11 = vm.Point2D((-0.5*self.tw, -0.5*self.h+self.tf))
        p12 = vm.Point2D((-0.5*self.b, -0.5*self.h+self.tf))

        rl = primitives2D.ClosedRoundedLineSegments2D([p1, p2, p3, p4, p5, p6, p7, p8,
                                                       p9, p10, p11, p12],
                                                      radius={3: self.r, 4: self.r,
                                                              9: self.r, 10: self.r})
        return rl
    
    def plot(self):
        self.contour().MPLPlot()
        
        
ipe_80 = ISection(0.08, 0.046, 0.0038, 0.0052, 0.005)
ipe_100 = ISection(0.1, 0.055, 0.0041, 0.0057, 0.007)
ipe_120 = ISection(0.12, 0.064, 0.0044, 0.0063, 0.007)
ipe_140 = ISection(0.14, 0.073, 0.0047, 0.0069, 0.007)
ipe_160 = ISection(0.16, 0.082, 0.005, 0.0074, 0.009)
ipe_180 = ISection(0.18, 0.091, 0.0053, 0.008, 0.009)
ipe_200 = ISection(0.20, 0.1, 0.0056, 0.0085, 0.012)
ipe_220 = ISection(0.22, 0.11, 0.0059, 0.0092, 0.012)
ipe_240 = ISection(0.24, 0.12, 0.0062, 0.0098, 0.015)
ipe_270 = ISection(0.27, 0.135, 0.0066, 0.0102, 0.015)
ipe_300 = ISection(0.30, 0.15, 0.0071, 0.0107, 0.015)
ipe_330 = ISection(0.33, 0.16, 0.0075, 0.0115, 0.018)
ipe_360 = ISection(0.36, 0.17, 0.008, 0.0127, 0.018)
ipe_400 = ISection(0.40, 0.18, 0.0086, 0.0135, 0.021)
ipe_450 = ISection(0.45, 0.19, 0.0094, 0.0146, 0.021)
ipe_500 = ISection(0.50, 0.20, 0.0102, 0.016, 0.021)
ipe_550 = ISection(0.55, 0.21, 0.0111, 0.0172, 0.024)
ipe_600 = ISection(0.60, 0.22, 0.012, 0.019, 0.024)
ipe_750_147 = ISection(0.753, 0.265, 0.0132, 0.017, 0.017)


#ipe_80.plot()
#ipe_100.plot()
#ipe_120.plot()
#ipe_140.plot()

class Beam(DessiaObject):
    def __init__(self, section, length, name=''):
        self.section = section
        self.length = length
        self.name = name
        
    def volmdlr_primitives(self, position=vm.O3D, x=vm.X3D, y=vm.Y3D):
        """
        x and y define the plane of section in the beam axis
        """
        z = x.Cross(y)
        c = self.section.contour()
        return [primitives3D.ExtrudedProfile(position-0.5*self.length*z, x, y, c, [], z*self.length)]
        
#beam = Beam(ipe_120, 2.56)
#beam.cad_volumes()