#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import math
import volmdlr as vm
import volmdlr.primitives3D
from typing import List

from dessia_common.core import DessiaObject

class HexagonNut(DessiaObject):
    _standalone_in_db = True
    _editable_variables = ['d', 't', 'h', 'e', 'name']
    _ordered_variables = ['d', 'name', 'h']
    _titled_variables = {'d' : 'd', 'name' : 'Name'}

    def __init__(self, d:float, t:float, h:float, name:str=''):
        self.d = d
        self.t = t
        self.e = 2*self.t / math.sqrt(3)
        self.h = h
        DessiaObject.__init__(self, name=name)

    def to_dict(self):
        dict_ = DessiaObject.to_dict(self)
        dict_.update({'d' : self.d, 't' : self.t, 'h' : self.h})
        return dict_

    @classmethod
    def dict_to_object(cls, dict_):
        return cls(d=dict_['d'], t=dict_['t'], h=dict_['h'], name=dict_['name'])

    def __hash__(self):
        return round(1000*self.d + 2123*self.t + 782*self.h)

    def __eq__(self, other_nut):
        return self.d == other_nut.d and self.t == other_nut.t\
            and self.h == other_nut.h

    def outer_contour(self):
        p1 = vm.Point2D((0., -0.5*self.e))
        p2 = vm.Point2D((0.5*self.t, -0.25*self.e))
        l1 = vm.LineSegment2D(p1, p2)
        p3 = vm.Point2D((0.5*self.t, 0.25*self.e))
        l2 = vm.LineSegment2D(p2, p3)
        p4 = vm.Point2D((0., 0.5*self.e))
        l3 = vm.LineSegment2D(p3, p4)
        p5 = vm.Point2D((-0.5*self.t, 0.25*self.e))
        l4 = vm.LineSegment2D(p4, p5)
        p6 = vm.Point2D((-0.5*self.t, -0.25*self.e))
        l5 = vm.LineSegment2D(p5, p6)
        l6 = vm.LineSegment2D(p6, p1)
        return vm.Contour2D([l1, l2, l3, l4, l5, l6])

    def inner_contour(self):
        return vm.Contour2D([vm.Circle2D(vm.o2D, 0.5*self.d)])

    def volmdlr_model(self, center=vm.o3D, x=vm.x3D, y=vm.y3D, z=vm.z3D):
        extrusion = volmdlr.primitives3D.ExtrudedProfile(center, x, y,
                                                         self.outer_contour(),
                                                         [self.inner_contour()],
                                                         z*self.h,
                                                         name=self.name)
        groups = [(self.name, [extrusion])]
        model = vm.VolumeModel(groups)
        return model

    def cad_export(self, fcstd_filepath='An unamed hexagon nut', python_path='python',
                   freecad_lib_path='/usr/lib/freecad/lib', export_types=['fcstd']):
        model = self.volmdlr_model()
        model.FreeCADExport(fcstd_filepath, python_path=python_path,
                            freecad_lib_path=freecad_lib_path, export_types=export_types)

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

    def to_dict(self):
        dict_ = {}
        dict_.update({'d' : self.d, 'L' : self.L, 'a' : self.a,
                      's' : self.s, 't' : self.t, 'name' : self.name})
        return dict_

    @classmethod
    def dict_to_object(cls, dict_):
        return cls(dict_['d'], dict_['L'], dict_['a'],
                   dict_['s'], dict_['t'], dict_['name'])


    def head_outer_contour(self):
        p1 = vm.Point2D((0., -0.5*self.e))
        p2 = vm.Point2D((0.5*self.t, -0.25*self.e))
        l1 = vm.LineSegment2D(p1, p2)
        p3 = vm.Point2D((0.5*self.t, 0.25*self.e))
        l2 = vm.LineSegment2D(p2, p3)
        p4 = vm.Point2D((0., 0.5*self.e))
        l3 = vm.LineSegment2D(p3, p4)
        p5 = vm.Point2D((-0.5*self.t, 0.25*self.e))
        l4 = vm.LineSegment2D(p4, p5)
        p6 = vm.Point2D((-0.5*self.t, -0.25*self.e))
        l5 = vm.LineSegment2D(p5, p6)
        l6 = vm.LineSegment2D(p6, p1)
        return vm.Contour2D([l1, l2, l3, l4, l5, l6])

    def body_outer_contour(self):
        return vm.Contour2D([vm.Circle2D(vm.o2D, 0.5*self.d)])

    def volmdlr_model(self, center=vm.o3D, x=vm.x3D, y=vm.y3D, z=vm.z3D):
        head = volmdlr.primitives3D.ExtrudedProfile(center, x, y,
                                                    self.head_outer_contour(),
                                                    [],
                                                    z*self.a,
                                                    name='head '+self.name)
        body_without_thead = volmdlr.primitives3D.ExtrudedProfile(center+z*self.a, x, y,
                                                                  self.body_outer_contour(),
                                                                  [],
                                                                  z*self.s,
                                                                  name='body '+self.name)
        body_with_thread = volmdlr.primitives3D.ExtrudedProfile(center+z*(self.a+self.s), x, y,
                                                                self.body_outer_contour(),
                                                                [],
                                                                z*(self.L-self.a),
                                                                name='thread '+self.name)
        groups = [(self.name, [head, body_without_thead, body_with_thread])]
        model = vm.VolumeModel(groups)
        return model

    def cad_export(self, fcstd_filepath='An unamed hexagon screw', python_path='python',
                   freecad_lib_path='/usr/lib/freecad/lib', export_types=['fcstd']):
        model = self.volmdlr_model()
        model.FreeCADExport(fcstd_filepath, python_path=python_path,
                            freecad_lib_path=freecad_lib_path, export_types=export_types)

class FlatWasher(DessiaObject):
    _standalone_in_db = True

    def __init__(self, D:float, A:float, e1:float, name=''):
        self.D = D
        self.A = A
        self.e1 = e1
        DessiaObject.__init__(self, name=name)

    def to_dict(self):
        dict_ = DessiaObject.to_dict(self)
        dict_.update({'D' : self.D, 'A' : self.A, 'e1' : self.e1})
        return dict_

    @classmethod
    def dict_to_object(cls, dict_):
        return cls(dict_['D'], dict_['A'], dict_['e1'], dict_['name'])


    def __hash__(self):
        return round(1000*self.A + 2123*self.D + 782*self.e1)

    def __eq__(self, other_nut):
        return self.A == other_nut.A and self.D == other_nut.D\
            and self.e1 == other_nut.e1

    def outer_contour(self):
        return vm.Contour2D([vm.Circle2D(vm.o2D, 0.5*self.A)])

    def inner_contour(self):
        return vm.Contour2D([vm.Circle2D(vm.o2D, 0.5*self.D)])

    def volmdlr_model(self, center=vm.o3D, x=vm.x3D, y=vm.y3D, z=vm.z3D):
        extrusion = volmdlr.primitives3D.ExtrudedProfile(center, x, y,
                                                         self.outer_contour(),
                                                         [self.inner_contour()],
                                                         z*self.e1,
                                                         name=self.name)
        groups = [(self.name, [extrusion])]
        model = vm.VolumeModel(groups)
        return model

    def cad_export(self, fcstd_filepath='An unamed flat washer', python_path='python',
                   freecad_lib_path='/usr/lib/freecad/lib', export_types=['fcstd']):
        model = self.volmdlr_model()
        model.FreeCADExport(fcstd_filepath, python_path=python_path,
                            freecad_lib_path=freecad_lib_path, export_types=export_types)

class Bolt(DessiaObject):
    def __init__(self, screw:HexagonScrew, nut:HexagonNut,
                 nut_position:float, washer:FlatWasher=None, name:str=''):
        self.screw = screw
        self.nut = nut
        self.nut_position = nut_position
        self.washer = washer

        DessiaObject.__init__(self, name=name)

    def to_dict(self):
        dict_ = DessiaObject.to_dict(self)
        dict_.update({'screw' : self.screw.to_dict(),
                      'nut' : self.nut.to_dict(),
                      'nut_position' : self.nut_position})
        if self.washer is not None:
            dict_['washer'] = self.washer.to_dict()
        return dict_

    @classmethod
    def dict_to_object(cls, dict_):
        screw = HexagonScrew.dict_to_object(dict_['screw'])
        nut = HexagonNut.dict_to_object(dict_['nut'])
        if "washer" in dict_:
            if dict_['washer'] is not None:
                washer = FlatWasher.dict_to_object(dict_['washer'])
                return cls(screw=screw, nut=nut, nut_position=dict_['nut_position'],
                           washer=washer, name=dict_['name'])
        return cls(screw=screw, nut=nut, nut_position=dict_['nut_position'], name=dict_['name'])


class ScrewAssembly(DessiaObject):
    def __init__(self, screws:List[HexagonScrew], positions:List[float],
                 axis:int, name:str=''):
        self.screws = screws
        self.positions = positions
        self.axis = axis
        DessiaObject.__init__(self, name=name)

    def to_dict(self):
        dict_ = DessiaObject.to_dict(self)
        dict_.update({'screws' : [s.to_dict() for s in self.screws],
                      'positions' : self.positions,
                      'axis' : self.axis})
        return dict_

    @classmethod
    def dict_to_object(cls, dict_):
        screws = [HexagonScrew.dict_to_object(s) for s in dict_['screws']]
        return cls(screws=screws, positions=dict_['positions'],
                   axis=dict_['axis'], name=dict_['name'])

class BoltAssembly(DessiaObject):
    def __init__(self, bolts:List[Bolt], positions:List[float],
                 axis:int, name:str=''):
        self.bolts = bolts
        self.positions = positions
        self.axis = axis
        DessiaObject.__init__(self, name=name)


    def to_dict(self):
        dict_ = DessiaObject.to_dict(self)
        dict_.update({'bolts' : [bolt.to_dict() for bolt in self.bolts],
                      'positions' : self.positions,
                      'axis' : self.axis})
        return dict_

    @classmethod
    def dict_to_object(cls, dict_):
        bolts = [Bolt.dict_to_object(s) for s in dict_['bolts']]
        return cls(bolts=bolts, positions=dict_['positions'],
                   axis=dict_['axis'], name=dict_['name'])
