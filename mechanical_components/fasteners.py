#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import math
import volmdlr as vm
import volmdlr.primitives3D

class HexagonNut:
    def __init__(self, d, t, h, name=''):
        self.d = d
        self.t = t
        self.e = 2*self.t / math.sqrt(3)
        self.h = h
        self.name = name
        
#        self.face_length = self.t / math.sqrt(3)
        
    def to_dict(self):
        d = {}
        d['d'] = self.d
        d['t'] = self.t
#        d['e'] = self.e
        d['h'] = self.h
        return d
    
    @classmethod
    def dict_to_object(cls, d):
        return cls(d['d'], d['t'], d['h'], d['name'])
    
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
        model=vm.VolumeModel(groups)        
        return model   
    
    def cad_export(self, fcstd_filepath='An unamed hexagon nut', python_path='python', 
                      freecad_lib_path='/usr/lib/freecad/lib', export_types=['fcstd']):
        model = self.volmdlr_model()
        model.FreeCADExport(fcstd_filepath, python_path=python_path,
                            freecad_lib_path=freecad_lib_path, export_types=export_types)

class Screw:
    def __init__(self, d, L, a, s, name=''):
        self.d = d
        self.L = L
        self.a = a
        self.s = s
        self.name = name
        
    def to_dict(self):
        d = {}
        d['d'] = self.d
        d['L'] = self.L
        d['a'] = self.a
        d['s'] = self.s
        d['name'] = self.name
        return d
    
    @classmethod
    def dict_to_object(cls, d):
        return cls(d['d'], d['L'], d['a'], d['s'], d['name'])

class FlatWasher:
    def __init__(self, D, A, e1, name=''):
        self.D = D
        self.A = A
        self.e1 = e1
        self.name = name
        
        
    def to_dict(self):
        d = {}
        d['D'] = self.D
        d['A'] = self.A
        d['e1'] = self.e1
        d['name'] = self.name
        return d
    
    @classmethod
    def dict_to_object(cls, d):
        return cls(d['D'], d['A'], d['e1'], d['name'])

class Bolt:
    
    def __init__(self, screw, nut, nut_position, washer=None):
        self.screw = screw
        self.nut = nut
        self.nut_position = nut_position
        self.washer = washer
        
    def to_dict(self):
        d = {}
        d['screw'] = self.screw.to_dict()
        d['nut'] = self.nut.to_dict()
        d['nut_position'] = self.nut_position
        if self.washer is not None:
            d['washer'] = self.washer.to_dict()
        return d
    
    @classmethod
    def dict_to_object(cls, d):
        screw = Screw.dict_to_object(d['screw'])
        nut = HexagonNut.dict_to_object(d['nut'])
        nut_position = d['nut_position']
        if "washer" in d:
            if d['washer'] is not None:
                washer = FlatWasher.dict_to_object(d['washer'])
                return cls(screw, nut, nut_position, washer)
        return cls(screw, nut, nut_position)
        
    
class ScrewAssembly:

    def __init__(self, screws, positions, axis, name=''):
        self.screws = screws
        self.positions = positions
        self.axis = axis
        self.name = name
    
    def to_dict(self):
        d = {}
        d['screws'] = [s.to_dict() for s in self.screws]
        d['positions'] = self.positions
        d['axis'] = self.axis
        d['name'] = self.name
        return d
    
    @classmethod
    def dict_to_object(cls, d):
        screws = [Screw.dict_to_object(s) for s in d['screws']]
        positions = d['positions']
        axis = d['axis']
        name = d['name']
        return cls(screws, positions, axis, name)
    
class BoltAssembly:
    
    def __init__(self, bolts, positions, axis, name=''):
        self.bolts = bolts
        self.positions = positions
        self.axis = axis
        self.name = name
        
        
    def to_dict(self):
        d = {}
        d['bolts'] = [bolt.to_dict() for bolt in self.bolts]
        d['positions'] = self.positions
        d['axis'] = self.axis
        d['name'] = self.name
        return d
    
    @classmethod
    def dict_to_object(cls, d):
        bolts = [Bolt.dict_to_object(s) for s in d['bolts']]
        positions = d['positions']
        axis = d['axis']
        name = d['name']
        return cls(bolts, positions, axis, name)
    
iso_nuts = [HexagonNut(1.6, 3.2, 1.3)]