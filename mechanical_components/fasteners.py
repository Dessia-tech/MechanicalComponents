#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import math
import volmdlr as vm
import volmdlr.edges as vme
import volmdlr.wires as vmw
import volmdlr.primitives3d
from typing import List

from dessia_common.core import DessiaObject

class HexagonNut(DessiaObject):
    _standalone_in_db = True

    def __init__(self, d:float, t:float, h:float, name:str=''):
        self.d = d
        self.t = t
        self.e = 2*self.t / math.sqrt(3)
        self.h = h
        DessiaObject.__init__(self, name=name)

    def __hash__(self):
        return round(1000*self.d + 2123*self.t + 782*self.h)

    def __eq__(self, other_nut):
        return self.d == other_nut.d and self.t == other_nut.t\
            and self.h == other_nut.h

    def outer_contour(self):
        p1 = vm.Point2D(0., -0.5*self.e)
        p2 = vm.Point2D(0.5*self.t, -0.25*self.e)
        l1 = vme.LineSegment2D(p1, p2)
        p3 = vm.Point2D(0.5*self.t, 0.25*self.e)
        l2 = vme.LineSegment2D(p2, p3)
        p4 = vm.Point2D(0., 0.5*self.e)
        l3 = vme.LineSegment2D(p3, p4)
        p5 = vm.Point2D(-0.5*self.t, 0.25*self.e)
        l4 = vme.LineSegment2D(p4, p5)
        p6 = vm.Point2D(-0.5*self.t, -0.25*self.e)
        l5 = vme.LineSegment2D(p5, p6)
        l6 = vme.LineSegment2D(p6, p1)
        return vmw.Contour2D([l1, l2, l3, l4, l5, l6])

    def inner_contour(self):
        return vmw.Circle2D(vm.O2D, 0.5*self.d)

    def volmdlr_primitives(self, center=vm.O3D, x=vm.X3D, y=vm.Y3D, z=vm.Z3D):
        extrusion = volmdlr.primitives3d.ExtrudedProfile(center, x, y,
                                                         self.outer_contour(),
                                                         [self.inner_contour()],
                                                         z*self.h,
                                                         name=self.name)
        return [extrusion]


class HexagonScrew(DessiaObject):
    _standalone_in_db = True

    def __init__(self, d:float, L:float, a:float, s:float, t:float, name:str=''):
        self.d = d
        self.L = L
        self.a = a
        self.s = s
        self.t = t
        self.e = 2*self.t / math.sqrt(3)
        DessiaObject.__init__(self, name=name)

    def __hash__(self):
        return round(1000*self.d + 2123*self.t + 782*self.L + 2839*self.s + 3829*self.a)

    def __eq__(self, other_screw):
        return self.d == other_screw.d and self.L == other_screw.L\
            and self.a == other_screw.a and self.s == other_screw.s\
            and self.t == other_screw.t


    def head_outer_contour(self):
        p1 = vm.Point2D(0., -0.5*self.e)
        p2 = vm.Point2D(0.5*self.t, -0.25*self.e)
        l1 = vme.LineSegment2D(p1, p2)
        p3 = vm.Point2D(0.5*self.t, 0.25*self.e)
        l2 = vme.LineSegment2D(p2, p3)
        p4 = vm.Point2D(0., 0.5*self.e)
        l3 = vme.LineSegment2D(p3, p4)
        p5 = vm.Point2D(-0.5*self.t, 0.25*self.e)
        l4 = vme.LineSegment2D(p4, p5)
        p6 = vm.Point2D(-0.5*self.t, -0.25*self.e)
        l5 = vme.LineSegment2D(p5, p6)
        l6 = vme.LineSegment2D(p6, p1)
        return vmw.Contour2D([l1, l2, l3, l4, l5, l6])

    def body_outer_contour(self):
        return vmw.Circle2D(vm.O2D, 0.5*self.d)

    def volmdlr_primitives(self, center=vm.O3D, x=vm.X3D, y=vm.Y3D, z=vm.Z3D):

        body_with_thread = volmdlr.primitives3d.ExtrudedProfile(center, x, y,
                                                                self.body_outer_contour(),
                                                                [],
                                                                z*(self.L-self.a),
                                                                name='thread '+self.name)
        body_without_thread = volmdlr.primitives3d.ExtrudedProfile(center+z*(self.L-self.a), x, y,
                                                                  self.body_outer_contour(),
                                                                  [],
                                                                  z*self.s,
                                                                  name='body '+self.name)
        head = volmdlr.primitives3d.ExtrudedProfile(center+z*(self.L-self.a+self.s), x, y,
                                                    self.head_outer_contour(),
                                                    [],
                                                    z*self.a,
                                                    name='head '+self.name)

        return [head, body_without_thread, body_with_thread]


class FlatWasher(DessiaObject):
    _standalone_in_db = True

    def __init__(self, D:float, A:float, e1:float, name=''):
        self.D = D
        self.A = A
        self.e1 = e1
        DessiaObject.__init__(self, name=name)


    def __hash__(self):
        return round(1000*self.A + 2123*self.D + 782*self.e1)

    def __eq__(self, other_nut):
        return self.A == other_nut.A and self.D == other_nut.D\
            and self.e1 == other_nut.e1

    def outer_contour(self):
        return vmw.Circle2D(vm.O2D, 0.5*self.A)

    def inner_contour(self):
        return [vmw.Circle2D(vm.O2D, 0.5*self.D)]

    def volmdlr_primitives(self, center=vm.O3D, x=vm.X3D, y=vm.Y3D, z=vm.Z3D):
        extrusion = volmdlr.primitives3d.ExtrudedProfile(center, x, y,
                                                         self.outer_contour(),
                                                         [self.inner_contour()],
                                                         z*self.e1,
                                                         name=self.name)
        return [extrusion]


class Bolt(DessiaObject):
    def __init__(self, screw:HexagonScrew, nut:HexagonNut,
                 nut_position:vm.Point3D, washer:FlatWasher=None, name:str=''):
        self.screw = screw
        self.nut = nut
        self.nut_position = nut_position
        self.washer = washer

        DessiaObject.__init__(self, name=name)

class ScrewAssembly(DessiaObject):
    def __init__(self, screws:List[HexagonScrew], positions:List[vm.Point3D],
                 axis:int, name:str=''):
        self.screws = screws
        self.positions = positions
        self.axis = axis
        DessiaObject.__init__(self, name=name)


class BoltAssembly(DessiaObject):
    def __init__(self, bolts:List[Bolt], positions:List[vm.Point3D],
                 axis:int, name:str=''):
        self.bolts = bolts
        self.positions = positions
        self.axis = axis
        DessiaObject.__init__(self, name=name)

class ScrewCatalog(DessiaObject):
    _standalone_in_db = True

    def __init__(self, screws:List[HexagonScrew], name:str=''):
        self.screws = screws
        self.name = name

    def screws_by_diameter(self):
        d = {}
        for s in self.screws:
            if s.d in d:
                d[s.d].append(s)
            else:
                d[s.d] = [s]
        return d

    def volmdlr_primitives(self):
        center = vm.O3D
        primitives = []

        screws_by_diameter = self.screws_by_diameter()

        for d in sorted(screws_by_diameter.keys()):
            center.x = 0
            for screw in screws_by_diameter[d]:
                primitives.extend(screw.volmdlr_primitives(center=center))
                center.x += 1.5*screw.t
            center.y += 2.5*screw.t
        return primitives