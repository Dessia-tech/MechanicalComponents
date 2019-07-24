#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import math
import volmdlr as vm
import volmdlr.primitives3D

class HexagonNut:
    
    _standalone_in_db = True
    
    _jsonschema = {
        "definitions": {},
        "$schema": "http://json-schema.org/draft-07/schema#",
        "type": "object",
        "title": "Hexagon nut",
        "required": [
            "d",
            "t",
            "h",
            'name'
          ],
          "properties": {
            "d": {"type": "number",  "examples": [0.008],
                  "physical_quantity": 'distance', "unit": "m",
                  "editable": True, "description": "Diameter"},
            "t": {"type": "number",  "examples": [0.012],
                  "physical_quantity": 'distance', "unit": "m",
                  "editable": True, "description": "Tool width"},
            "h": {"type": "number",  "examples": [0.008],
                  "physical_quantity": 'distance', "unit": "m",
                  "editable": True,
                  "description": "Height"},
            "e": {"type": "number", "examples": [0.003],
                  "physical_quantity": 'distance', "unit": "m",
                  "editable": False},
            "name": {'type': 'string', 'editable': True, "examples": ['M4'],
                     }
            }
          }
    
    def __init__(self, d, t, h, name=''):
        self.d = d
        self.t = t
        self.e = 2*self.t / math.sqrt(3)
        self.h = h
        self.name = name
        
    def to_dict(self):
        d = {}
        d['d'] = self.d
        d['t'] = self.t
        d['h'] = self.h
        d['name'] = self.name
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

class HexagonScrew:
    _standalone_in_db = True
    
    _jsonschema = {
        "definitions": {},
        "$schema": "http://json-schema.org/draft-07/schema#",
        "type": "object",
        "title": "Hexagon screw",
        "required": [
            "d",
            "L",
            "a",
            's',
            't',
            'name'
          ],
          "properties": {
            "d": {"type": "number",  "examples": [0.008],
                  "physical_quantity": 'distance', "unit": "m",
                  "editable": True, "description": "Diameter"},
            "L": {"type": "number",  "examples": [0.012],
                  "physical_quantity": 'distance', "unit": "m",
                  "editable": True, "description": "Body length"},
            "a": {"type": "number",  "examples": [0.008],
                  "physical_quantity": 'distance', "unit": "m",
                  "editable": True,
                  "description": "Head length"},
            "s": {"type": "number",  "examples": [0.004],
                  "physical_quantity": 'distance', "unit": "m",
                  "editable": True,
                  "description": "Body length without thread"},
            "t": {"type": "number",  "examples": [0.012],
                  "physical_quantity": 'distance', "unit": "m",
                  "editable": True,
                  "description": "Tool width"},
            "e": {"type": "number", "examples": [0.003],
                  "physical_quantity": 'distance', "unit": "m",
                  "editable": False},
            "name": {'type': 'string', 'editable': True, "examples": ['M4'],
                     }
            }
          }
    
    
    def __init__(self, d, L, a, s, t, name=''):
        self.d = d
        self.L = L
        self.a = a
        self.s = s
        self.t = t
        self.e = 2*self.t / math.sqrt(3)
        self.name = name
        
    def __hash__(self):
        return round(1000*self.d + 2123*self.t + 782*self.L + 2839*self.s + 3829*self.a)

    def __eq__(self, other_screw):
        return self.d == other_screw.d and self.L == other_screw.L\
            and self.a == other_screw.a and self.s == other_screw.s\
            and self.t == other_screw.t
        
    def to_dict(self):
        d = {}
        d['d'] = self.d
        d['L'] = self.L
        d['a'] = self.a
        d['s'] = self.s
        d['t'] = self.t
        d['name'] = self.name
        return d
    
    @classmethod
    def dict_to_object(cls, d):
        return cls(d['d'], d['L'], d['a'], d['s'], d['t'], d['name'])


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
        model=vm.VolumeModel(groups)        
        return model   
    
    def cad_export(self, fcstd_filepath='An unamed hexagon nut', python_path='python', 
                      freecad_lib_path='/usr/lib/freecad/lib', export_types=['fcstd']):
        model = self.volmdlr_model()
        model.FreeCADExport(fcstd_filepath, python_path=python_path,
                            freecad_lib_path=freecad_lib_path, export_types=export_types)

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
    
iso_nuts = [HexagonNut(0.0016, 0.0032, 0.0013)]
iso_hexscrews = [HexagonScrew(0.008, 0.025, 0.008, 0.005, 0.012)]