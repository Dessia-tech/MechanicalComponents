#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 10:30:12 2018

@author: steven
"""

class Nut:
    def __init__(self, d, t, e, h, name=''):
        self.d = d
        self.t = t
        self.e = e
        self.h = h
        self.name = name
        
    def Dict(self):
        d = self.__dict__.copy()
        return d
    
    @classmethod
    def DictToObject(cls, d):
        return cls(d['d'], d['t'], d['e'], d['h'], d['name'])

class Screw:
    def __init__(self, d, L, a, s, name=''):
        self.d = d
        self.L = L
        self.a = a
        self.s = s
        self.name = name
        
    def Dict(self):
        d = self.__dict__.copy()
        return d
    
    @classmethod
    def DictToObject(cls, d):
        return cls(d['d'], d['L'], d['a'], d['s'], d['name'])

class FlatWasher:
    def __init__(self, D, A, e1, name=''):
        self.D = D
        self.A = A
        self.e1 = e1
        self.name = name
        
        
    def Dict(self):
        d = self.__dict__.copy()
        return d
    
    @classmethod
    def DictToObject(cls, d):
        return cls(d['D'], d['A'], d['e1'], d['name'])

class Bolt:
    
    dessia_db_attributes = [{'name':'screw',
                     'class':'mechanical_components.fasteners.Screw',
                     'type':'object'}]
    
    dessia_db_attributes = [{'name':'nut',
                 'class':'mechanical_components.fasteners.Nut',
                 'type':'object'}]
    
    dessia_db_attributes = [{'name':'washer',
                     'class':'mechanical_components.fasteners.Washer',
                     'type':'object'}]
    
    def __init__(self, screw, nut, nut_position, washer=None):
        self.screw = screw
        self.nut = nut
        self.nut_position = nut_position
        self.washer = washer
        
    def Dict(self):
        d = {}
        d['screw'] = self.screw.Dict()
        d['nut'] = self.nut.Dict()
        d['nut_position'] = self.nut_position
        if self.washer is not None:
            d['washer'] = self.washer.Dict()
        return d
    
    @classmethod
    def DictToObject(cls, d):
        screw = Screw.DictToObject(d['screw'])
        nut = Nut.DictToObject(d['screw'])
        nut_position = d['nut_position']
        if d['washer'] is not None:
            washer = FlatWasher.DictToObject(d['washer'])
        return cls(screw, nut, nut_position, washer)
        
    
class ScrewAssembly:
    
    dessia_db_attributes = [{'name':'screws',
                         'class':'mechanical_components.fasteners.Screw',
                         'type':'list'}]

    def __init__(self, screws, positions, axis, name=''):
        self.screws = screws
        self.positions = positions
        self.axis = axis
        self.name = name
    
    def Dict(self):
        d = {}
        d['screws'] = [s.Dict() for s in self.screws]
        d['positions'] = self.positions
        d['axis'] = self.axis
        d['name'] = self.name
        return d
    
    @classmethod
    def DictToObject(cls, d):
        screws = [Screw.DictToObject(s) for s in d['screws']]
        positions = d['positions']
        axis = d['axis']
        name = d['name']
        return cls(screws, positions, axis, name)
    
class BoltAssembly:
    
    dessia_db_attributes = [{'name':'bolts',
                     'class':'mechanical_components.fasteners.Bolt',
                     'type':'list'}]
    
    def __init__(self, bolts, positions, axis, name=''):
        self.bolts = bolts
        self.positions = positions
        self.axis = axis
        self.name = name
        
        
    def Dict(self):
        d = {}
        d['bolts'] = [s.Dict() for s in self.bolts]
        d['positions'] = self.positions
        d['axis'] = self.axis
        d['name'] = self.name
        return d
    
    @classmethod
    def DictToObject(cls, d):
        bolts = [Screw.DictToObject(s) for s in d['bolts']]
        positions = d['positions']
        axis = d['axis']
        name = d['name']
        return cls(bolts, positions, axis, name)
    
    
iso_nuts = [Nut(1.6, 3.2, 3.4, 1.3),
        Nut(2, 4, 4.4, 1.6),
        Nut(2.5, 5, 5.4, 2),
        Nut(3, 5.5, 6, 2.4),
        Nut(4, 7, 7.6, 3.2)]


screws = [Screw(1.6, 16, 0.35, 1.5),
          Screw(2, 18, 0.35, 1.5),
          Screw(2.5, 20, 0.35, 1.5)]

washers = [FlatWasher(5, 1.7, 0.5),
           FlatWasher(6, 2.2, 0.5),
           FlatWasher(7, 2.7, 0.5)]


bolts =[Bolt(screws[0], iso_nuts[0], 1.3),
        Bolt(screws[1], iso_nuts[1], 1.3, washers[0])]

bolts_assembly = BoltAssembly(bolts, [(0,0,0), (0.1, 0, 0)], (1,0,0))
screw_assembly = ScrewAssembly(screws, [(0,0,0), (0.1, 0, 0), (-0.05, 0.05, 0)], (1,0,0))