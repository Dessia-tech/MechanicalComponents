#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 02:13:01 2018

@author: Pierre-Emmanuel Dumouchel
"""

import numpy as npy
npy.seterr(all='raise', under='ignore')
from scipy import interpolate

import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D

import math
from scipy.linalg import norm
from scipy.optimize import fsolve, minimize
import networkx as nx

from dessia_common import DessiaObject
import mechanical_components.tools as tools
import json
import copy

#data_coeff_YB_Iso
evol_coeff_yb_iso={'data':[[0.0,1.0029325508401201],
                           [4.701492563229561,0.9310850480431024],
                           [9.104477651442416,0.8782991233021732],
                           [14.5522388104227,0.8240469255458759],
                           [19.328358165913905,0.784457481990179],
                           [23.955224059269884,0.7609970656504502],
                           [28.507462609617956,0.7521994155347822],
                           [34.029850499665855,0.7507331425194141],
                           [40.0,0.7492668574805859]
                          ], 'x':'Linear','y':'Linear'}

#data_wholer_curve
wholer_hardened_alloy_steel={'data':[[4.296196199237153,1.9797762011105589],
                                     [4.824840106199563,1.9413306094362142],
                                     [5.3344338175705674,1.908892154601565],
                                     [6.115493253679078,1.8632380197445122],
                                     [6.596511629990596,1.8560294765618042],
                                     [7.144205815889171,1.8536266428508523],
                                     [7.691899918442984,1.8524252154829133],
                                     [8.010991340520903,1.8524252154829133]
                                    ], 'x':'Log','y':'Log'}

wholer_nitrided_alloy_steel={'data':[[4.104865629699472,1.9252942042661974],
                                     [4.568697315952783,1.8521640228225367],
                                     [4.887581626297173,1.8046294185503593],
                                     [5.438381821440599,1.7900033864666123],
                                     [6.402282079596832,1.7918316299646175],
                                     [7.264719174616821,1.7918316299646175],
                                     [7.989456220850952,1.793659894487549]
                                    ], 'x':'Log','y':'Log'}

wholer_through_hardened_steel={'data':[[4.172369719531124,1.895676495604088],
                                       [4.677200861168087,1.7983611100752137],
                                       [4.9677168648417585,1.741894170956562],
                                       [5.329671247836526,1.6842258044699714],
                                       [5.439210101685194,1.672211551815507],
                                       [6.091680488353632,1.6734129791834462],
                                       [7.139443246155129,1.671010124447568],
                                       [8.00146620105282,1.6758158339193243]
                                      ],'x':'Log','y':'Log'}

wholer_surface_hardened_steel={'data':[[4.281908490035029,1.7611169667937343],
                                       [4.701013626493532,1.6998443182033265],
                                       [5.015342395492649,1.6553916107142128],
                                       [5.358246582896013,1.6109389032250994],
                                       [5.620187251510196,1.5857089915731581],
                                       [6.020242109032534,1.5748961873115592],
                                       [6.567936294931109,1.5748961873115592],
                                       [7.263269725861159,1.5772990210225108],
                                       [7.996703631318779,1.5772990210225108]
                                      ], 'x':'Log','y':'Log'}

wholer_carbon_steel={'data':[[4.307791955971963,1.6419147590563592],
                             [5.242702822291173,1.535876005424268],
                             [5.938450393343521,1.4700588400224806],
                             [6.518240063668731,1.431665495290182],
                             [7.221234961844144,1.4334937598131132],
                             [7.989456220850952,1.4353220033111185]
                            ], 'x':'Log','y':'Log'}

wholer_cast_iron={'data':[[4.307791955971963,1.6419147590563592],
                          [5.242702822291173,1.535876005424268],
                          [5.938450393343521,1.4700588400224806],
                          [6.518240063668731,1.431665495290182],
                          [7.221234961844144,1.4334937598131132],
                          [7.989456220850952,1.4353220033111185]
                         ], 'x':'Log','y':'Log'}

wholer_bronze={'data':[[4.307791955971963,1.6419147590563592],
                       [5.242702822291173,1.535876005424268],
                       [5.938450393343521,1.4700588400224806],
                       [6.518240063668731,1.431665495290182],
                       [7.221234961844144,1.4334937598131132],
                       [7.989456220850952,1.4353220033111185]
                      ], 'x':'Log','y':'Log'}

wholer_grey_iron={'data':[[4.307791955971963,1.6419147590563592],
                          [5.242702822291173,1.535876005424268],
                          [5.938450393343521,1.4700588400224806],
                          [6.518240063668731,1.431665495290182],
                          [7.221234961844144,1.4334937598131132],
                          [7.989456220850952,1.4353220033111185]
                         ], 'x':'Log','y':'Log'}
#data_gear_material

sigma_hardened_alloy_steel={'data':[[1.8422104370714443,1.4645831828946267],
                                    [1.948612010770208,1.5219116983411152],
                                    [2.0605171321606295,1.5810895335609718]
                                    ,[2.141235568740199,1.6254729099758645]
                                   ], 'x':'Log','y':'Log'}

sigma_nitrided_alloy_steel={'data':[[1.8458794622934307,1.4349942652846983],
                                    [1.943108482795906,1.488624180937243],
                                    [2.0201578941534892,1.5274596179084272],
                                    [2.128393990321924,1.5866374531282839]
                                   ],'x':'Log','y':'Log'}

sigma_through_hardened_steel={'data':[[1.7798371068844516,1.292597616678765],
                                      [1.921094370898698,1.3850629693024938],
                                      [2.032999472571764,1.4627338829976548],
                                      [2.1650841833897223,1.5533499158480155]
                                     ],'x':'Log','y':'Log'}

sigma_surface_hardened_steel={'data':[[1.8312033811228403,1.115064130895591],
                                      [1.932101426847302,1.200132264055036],
                                      [2.038503000546066,1.2852003773380847]
                                     ], 'x':'Log','y':'Log'}

sigma_carbon_steel={'data':[[1.677104538690319,1.1002696720906269],
                            [1.7633265032441903,1.1723926463420797],
                            [1.8385414118494579,1.2389677010262203],
                            [1.8844041581135444,1.2796524577707729]
                           ], 'x':'Log','y':'Log'}

sigma_cast_iron={'data':[[1.4734739247717241,0.922736186307453],
                         [1.5468543306246763,0.9837633214242817],
                         [1.6073931580593532,1.0336946174064863],
                         [1.6404143456225206,1.0688314545837265]
                        ], 'x':'Log','y':'Log'}

sigma_bronze={'data':[[1.313871566195314,0.7858874572688317],
                      [1.3890864826875238,0.8487638922826322],
                      [1.4294457009773085,0.8802021097895326],
                      [1.4551288380965028,0.9097910273994609]
                     ], 'x':'Log','y':'Log'}

sigma_grey_iron={'data':[[1.354230792372041,0.7100658633470387],
                         [1.4276111785076375,0.7766409180311793],
                         [1.4936535339166166,0.84691459238566],
                         [1.5431853054026896,0.8986951882648367],
                         [1.5725374677438706,0.933832025442077]
                        ], 'x':'Log','y':'Log'}

class Material(DessiaObject):
    """
    Gear material

    :param volumic_mass: A float to define the gear volumic mass
    :param data_coeff_YB_Iso: a dictionary to define the YB parameter of the ISO description
    :param data_wholer_curve: a dictionary to define the wholer slope of the ISO description
    :param data_gear_material: a dictionary to define the maximum gear stress

    :data_coeff_YB_Iso: - **'data'** matrix define points of the YB curve in the plane (YB, helix_angle)
        - **'x'** string define the x axis evolution ('Log' or 'Linear')
        - **'y'** string define the y axis evolution ('Log' or 'Linear')

    :data_wholer_curve: - **'data'** matrix define points of the wholer slope in the plane (wholer slope, number of cycle)
        - **'x'** string define the x axis evolution ('Log' or 'Linear')
        - **'y'** string define the y axis evolution ('Log' or 'Linear')

    :data_gear_material: - **'data'** matrix define points of the maximum gear stress (maximum gear stress, wholer slope)
        - **'x'** string define the x axis evolution ('Log' or 'Linear')
        - **'y'** string define the y axis evolution ('Log' or 'Linear')

    >>> volumic_mass=7800
    >>> data_coeff_YB_Iso={'data':[[0.0,1.0029325508401201],
                           [4.701492563229561,0.9310850480431024],
                           [23.955224059269884,0.7609970656504502],
                           [40.0,0.7492668574805859]
                          ], 'x':'Linear','y':'Linear'}
    >>> data_wholer_curve={'data':[[4.307791955971963,1.6419147590563592],
                       [6.518240063668731,1.431665495290182],
                       [7.989456220850952,1.4353220033111185]
                      ], 'x':'Log','y':'Log'}
    >>> data_gear_material={'data':[[1.313871566195314,0.7858874572688317],
                      [1.4294457009773085,0.8802021097895326],
                      [1.4551288380965028,0.9097910273994609]
                     ], 'x':'Log','y':'Log'}
    >>> material1=Material(volumic_mass, data_coeff_YB_Iso,
                           data_wholer_curve, data_gear_material)
    """
    _standalone_in_db = False

    def __init__(self, volumic_mass, data_coeff_YB_Iso, data_wholer_curve,
                 data_gear_material, name=''):
        self.volumic_mass = volumic_mass
        self.data_coeff_YB_Iso = data_coeff_YB_Iso
        self.data_wholer_curve = data_wholer_curve
        self.data_gear_material = data_gear_material

        DessiaObject.__init__(self, name=name)

    def __eq__(self, other_eb):
        equal = (self.volumic_mass == other_eb.volumic_mass
                 and self.data_coeff_YB_Iso == other_eb.data_coeff_YB_Iso
                 and self.data_wholer_curve == other_eb.data_wholer_curve
                 and self.data_gear_material == other_eb.data_gear_material)
        return equal

    def __hash__(self):
        material_hash = hash(self.volumic_mass)
        return material_hash

    def FunCoeff(self,x,data,type_x='Linear',type_y='Linear'):
        """ Interpolation of material data

        :param x: value of the interpolation
        :param data: dictionary of the input data
        :param type_x: type of the x axis of the data matrix ('Log' or 'Linear')
        :param type_y: type of the y axis of the data matrix ('Log' or 'Linear')

        :returns:  interpolation value

        >>> interp1=material1.FunCoeff(x = 5.2,data = data_wholer_curve,
                                       type_x = 'Log',type_y = 'Log')
        """
        if type_x == 'Log':
            x = math.log10(x)
        f = interpolate.interp1d(list(data[:,0]),list(data[:,1]),
                                 fill_value='extrapolate')
        sol=float(f(x))
        if type_y=='Log':
            sol=10**sol
        return sol

    def Dict(self):

        d = {'name' : self.name} # TODO Change this to DessiaObject.__init__
        d['volumic_mass'] = self.volumic_mass
        d['data_coeff_YB_Iso'] = self.data_coeff_YB_Iso
        d['data_wholer_curve'] = self.data_wholer_curve
        d['data_gear_material'] = self.data_gear_material
        return d

    @classmethod
    def DictToObject(cls, d):
        material = cls(volumic_mass = d['volumic_mass'],
                   data_coeff_YB_Iso = d['data_coeff_YB_Iso'],
                   data_wholer_curve = d['data_wholer_curve'],
                   data_gear_material = d['data_gear_material'],
                   name=d['name'])
        return material

hardened_alloy_steel=Material(7850, evol_coeff_yb_iso,
                              wholer_hardened_alloy_steel,
                              sigma_hardened_alloy_steel,
                              name='Hardened alloy steel')

nitrided_alloy_steel=Material(7850, evol_coeff_yb_iso, wholer_nitrided_alloy_steel,
                              sigma_nitrided_alloy_steel,
                              name='Nitrided alloy steel')

through_hardened_steel=Material(7850, evol_coeff_yb_iso,
                                wholer_through_hardened_steel,
                                sigma_through_hardened_steel,
                                name='Through hardened steel')

surface_hardened_steel=Material(7850, evol_coeff_yb_iso,
                                wholer_surface_hardened_steel,
                                sigma_surface_hardened_steel,
                                name='Surface hardened steel')

carbon_steel=Material(7850, evol_coeff_yb_iso, wholer_carbon_steel, sigma_carbon_steel,
                      name='Carbon steel')

cast_iron=Material(7200, evol_coeff_yb_iso, wholer_cast_iron, sigma_cast_iron,
                   name='Cast iron')

bronze=Material(8200, evol_coeff_yb_iso, wholer_bronze, sigma_bronze,
                name='Bronze')

grey_iron=Material(7200, evol_coeff_yb_iso, wholer_grey_iron, sigma_grey_iron,
                   name='Grey iron')

class Rack(DessiaObject):
    """
    Gear rack definition

    :param transverse_pressure_angle: definition of the transverse pressure angle of the rack
    :type transverse_pressure_angle: radian

    >>> Rack1=Rack(20/180.*math.pi) #definition of an ISO rack
    """
    _standalone_in_db = True
    _jsonschema = {
        "definitions": {},
        "$schema": "http://json-schema.org/draft-07/schema#",
        "type": "object",
        "title": "mechanical_components.meshes.Rack Base Schema",
        "required": ['transverse_pressure_angle'],
        "properties": {
            'transverse_pressure_angle': {
                "type": "number",
                "step" : 0.01,
                "examples": [0.34],
                "editable": True,
                "description": "Transverse pressure angle",
                "order" : 1,
                "unit": 'rad',
                },
            'module': {
                "type": "number",
                "step" : 1,
                "examples": [2],
                "editable": True,
                "description": "Module",
                "order" : 2
                },
            'coeff_gear_addendum': {
                "type": "number",
                "step" : 0.01,
                "examples": [1],
                "editable": True,
                "description": "Coeff gear addendum",
                "order" : 3
                },
            'coeff_gear_dedendum': {
                "type": "number",
                "step" : 0.01,
                "examples": [1.25],
                "editable": True,
                "description": "Coeff gear dedendum",
                "order" : 4
                },
            'coeff_root_radius': {
                "type": "number",
                "step" : 0.01,
                "examples": [0.38],
                "editable": True,
                "description": "Coeff root radius",
                "order" : 5
                },
            'coeff_circular_tooth_thickness': {
                "type": "number",
                "step" : 0.01,
                "examples": [0.5],
                "editable": True,
                "description": "Coeff circular tooth thickness",
                "order" : 6
                },
             }
         }

    _display_angular = []

    def __init__(self, transverse_pressure_angle, module=None,
                 coeff_gear_addendum=1, coeff_gear_dedendum=1.25,
                 coeff_root_radius=0.38, coeff_circular_tooth_thickness=0.5, name=''):
        self.transverse_pressure_angle = transverse_pressure_angle
        self.module = module
        self.coeff_gear_addendum = coeff_gear_addendum
        self.coeff_gear_dedendum = coeff_gear_dedendum
        self.coeff_root_radius = coeff_root_radius
        self.coeff_circular_tooth_thickness = coeff_circular_tooth_thickness


        DessiaObject.__init__(self, name=name)

        if module is not None:
            self.Update(module, transverse_pressure_angle, coeff_gear_addendum,
                        coeff_gear_dedendum, coeff_root_radius,
                        coeff_circular_tooth_thickness)

    def __eq__(self, other_eb):
        equal = (self.transverse_pressure_angle == other_eb.transverse_pressure_angle
                 and self.coeff_gear_addendum == other_eb.coeff_gear_addendum
                 and self.coeff_gear_dedendum == other_eb.coeff_gear_dedendum
                 and self.coeff_root_radius == other_eb.coeff_root_radius
                 and self.coeff_circular_tooth_thickness == other_eb.coeff_circular_tooth_thickness)
        return equal

    def __hash__(self):
        rack_hash = hash(self.transverse_pressure_angle) + hash(self.coeff_gear_addendum)\
                 + hash(self.coeff_gear_dedendum) + hash(self.coeff_root_radius)\
                 + hash(self.coeff_circular_tooth_thickness)
        return rack_hash

    def RackParam(self, transverse_pressure_angle, coeff_gear_addendum,
                  coeff_gear_dedendum, coeff_root_radius, coeff_circular_tooth_thickness):

        self.transverse_pressure_angle = transverse_pressure_angle
        self.transverse_radial_pitch = self.module*math.pi
        self.gear_addendum = coeff_gear_addendum*self.module
        self.gear_dedendum = coeff_gear_dedendum*self.module
        self.root_radius = coeff_root_radius*self.module
        self.circular_tooth_thickness = coeff_circular_tooth_thickness*self.transverse_radial_pitch

        self.tooth_space = self.transverse_radial_pitch-self.circular_tooth_thickness
        self.whole_depth = self.gear_addendum+self.gear_dedendum
        self.clearance = self.root_radius-self.root_radius*math.sin(self.transverse_pressure_angle)


        # trochoide parameter
        self.a = (self.tooth_space/2.
                 - self.gear_dedendum * math.tan(self.transverse_pressure_angle)
                 - self.root_radius * math.tan(0.5*math.atan(math.cos(self.transverse_pressure_angle)
                                                 /(math.sin(self.transverse_pressure_angle)))))
        self.b = self.gear_dedendum - self.root_radius

    def Update(self, module, transverse_pressure_angle=None, coeff_gear_addendum=None,
               coeff_gear_dedendum=None, coeff_root_radius=None,
               coeff_circular_tooth_thickness=None):
        """
        Update of the gear rack

        :param module: update of the module of the rack define on the pitch factory diameter
        :type module: m
        :param transverse_pressure_angle: update of the transverse pressure angle of the rack
        :type transverse_pressure_angle: radian
        :param coeff_gear_addendum: update of the gear addendum coefficient (gear_addendum = coeff_gear_addendum*module)
        :param coeff_gear_dedendum: update of the gear dedendum coefficient (gear_dedendum = coeff_gear_dedendum*module) (top of the rack)
        :param coeff_root_radius: update of the root radius coefficient (root_radius = coeff_root_radius*module)
        :param coeff_circular_tooth_thickness: update of the circular tooth thickness coefficient (circular_tooth_thickness = coeff_circular_tooth_thickness*transverse_radial_pitch)

        >>> input={'module':2*1e-3,'transverse_pressure_angle':21/180.*math.pi}
        >>> Rack1.Update(**input) # Update of the rack definition
        """
        if transverse_pressure_angle == None:
            transverse_pressure_angle = self.transverse_pressure_angle
        if coeff_gear_addendum == None:
            coeff_gear_addendum = self.coeff_gear_addendum
        if coeff_gear_dedendum == None:
            coeff_gear_dedendum = self.coeff_gear_dedendum
        if coeff_root_radius == None:
            coeff_root_radius = self.coeff_root_radius
        if coeff_circular_tooth_thickness == None:
            coeff_circular_tooth_thickness = self.coeff_circular_tooth_thickness
        self.module = module

        self.RackParam(transverse_pressure_angle, coeff_gear_addendum,
                       coeff_gear_dedendum, coeff_root_radius, coeff_circular_tooth_thickness)

    ### Optimization Method

    def CheckRackViable(self):
        """ Check the viability of the rack toward the top and the root

        :results: boolean variable, and a list of element to be positive for the optimizer
        """
        list_ineq=[]
        list_ineq.append(self.transverse_radial_pitch-self.circular_tooth_thickness
                         -2*self.gear_dedendum*math.tan(self.transverse_pressure_angle)
                         -2*(self.root_radius*math.cos(self.transverse_pressure_angle)-math.tan(self.transverse_pressure_angle)
                         *self.root_radius*(1-math.sin(self.transverse_pressure_angle))))
        list_ineq.append(self.circular_tooth_thickness-2*(self.gear_addendum*math.tan(self.transverse_pressure_angle)))
        check=False
        if min(list_ineq)>0:
            check=True
        return check,list_ineq

    def ListeIneq(self):
        """ Compilation method for inequality list used by the optimizer

        :results: vector of data that should be positive
        """
        check,ineq=self.CheckRackViable
        return ineq

    def Contour(self,number_pattern):
        """ Construction of the volmdr 2D rack profile

        :param number_pattern: number of rack pattern to define
        """
        p1=vm.Point2D((0,0))
        p2=p1.Translation((self.gear_addendum*math.tan(self.transverse_pressure_angle),self.gear_addendum))
        p4=p1.Translation((self.circular_tooth_thickness,0))
        p3=p4.Translation((-self.gear_addendum*math.tan(self.transverse_pressure_angle),self.gear_addendum))
        p5=p4.Translation((self.gear_dedendum*math.tan(self.transverse_pressure_angle),-self.gear_dedendum))
        p7=p4.Translation((self.tooth_space,0))
        p6=p7.Translation((-self.gear_dedendum*math.tan(self.transverse_pressure_angle),-self.gear_dedendum))
        L=primitives2D.OpenedRoundedLineSegments2D([p1,p2,p3,p4,p5,p6,p7],{4:self.root_radius,5:self.root_radius},False)

        Rack_Elem=[]
        for i in range(number_pattern):
            Rack_Elem.append(L.Translation(((i)*(p7.vector-p1.vector))))
        p10=Rack_Elem[0].points[0]
        p15=Rack_Elem[-1].points[-1]
        p11=p10.Translation((-self.circular_tooth_thickness,0))
        p12=p11.Translation((0,2*self.whole_depth))
        p14=p15.Translation((self.circular_tooth_thickness,0))
        p13=p14.Translation((0,2*self.whole_depth))
        Rack_Elem.append(primitives2D.OpenedRoundedLineSegments2D([p10,p11,p12,p13,p14,p15],{},False))

        return Rack_Elem

    def Plot(self,number_pattern):
        """ Plot function of the rack

        :param number_pattern: number of rack pattern to draw
        """
        Rack_Elem=self.Contour(number_pattern)
        RackElem=vm.Contour2D(Rack_Elem)
        RackElem.MPLPlot()

    def CSVExport(self):
        """
        Export CSV format

         :returns:  list of all element in dict() function
        """
        d=self.__dict__.copy()
        return list(d.keys()),list(d.values())

    def Dict(self):
        d = {'name' : self.name} # TODO Change this to DessiaObject.__init__
        d['transverse_pressure_angle'] = self.transverse_pressure_angle
        d['module'] = self.module
        d['coeff_gear_addendum'] = self.coeff_gear_addendum
        d['coeff_gear_dedendum'] = self.coeff_gear_dedendum
        d['coeff_root_radius'] = self.coeff_root_radius
        d['coeff_circular_tooth_thickness'] = self.coeff_circular_tooth_thickness
        return d

    @classmethod
    def DictToObject(cls, d):
        rack = cls(transverse_pressure_angle = d['transverse_pressure_angle'],
                   module = d['module'],
                   coeff_gear_addendum = d['coeff_gear_addendum'],
                   coeff_gear_dedendum = d['coeff_gear_dedendum'],
                   coeff_root_radius = d['coeff_root_radius'],
                   coeff_circular_tooth_thickness = d['coeff_circular_tooth_thickness'],
                   name=d['name'])
        return rack

class Mesh(DessiaObject):
    """
    Gear mesh definition

    :param z: number of tooth
    :param db: base diameter
    :type db: m
    :param cp: coefficient profile shift of the rack
    :param transverse_pressure_angle_rack: transverse pressure angle of the rack
    :type transverse_pressure_angle_rack: radian
    :param coeff_gear_addendum: update of the gear addendum coefficient (gear_addendum = coeff_gear_addendum*module)
    :param coeff_gear_dedendum: update of the gear dedendum coefficient (gear_dedendum = coeff_gear_dedendum*module)
    :param coeff_root_radius: update of the root radius coefficient (root_radius = coeff_root_radius*module)
    :param coeff_circular_tooth_thickness: update of the circular tooth thickness coefficient (circular_tooth_thickness = coeff_circular_tooth_thickness*transverse_radial_pitch)
    :param material: class material define the gear mesh material
    :param gear_width: gear mesh width

    >>> input={z:13, db:40*1e-3, cp:0.3, transverse_pressure_angle_rack:20/180.*math.pi,
                 coeff_gear_addendum:1, coeff_gear_dedendum:1, coeff_root_radius:1,
                 coeff_circular_tooth_thickness:1}
    >>> mesh1=Mesh(**input) # generation of one gear mesh
    """
    _standalone_in_db = True
    _jsonschema = {
        "definitions": {},
        "$schema": "http://json-schema.org/draft-07/schema#",
        "type": "object",
        "title": "mechanical_components.meshes.Mesh Base Schema",
        "required": ['z', 'db', 'coefficient_profile_shift', 'rack'],
        "properties": {
            'z': {
                "type": "number",
                "step" : 1,
                "examples": [15],
                "editable": True,
                "description": "Number tooth",
                "order" : 1,
                },
            'db': {
                "type": "number",
                "step" : 0.001,
                "examples": [0.06],
                "editable": True,
                "description": "Base diameter",
                "order" : 2,
                "unit" : "m",
                },
            'coefficient_profile_shift': {
                "type": "number",
                "step" : 0.01,
                "examples": [0.1],
                "editable": True,
                "description": "Coefficient profile shift",
                "order" : 3
                },
            "rack" : {
                "type" : "object",
                "classes" : ["mechanical_components.meshes.Rack"],
                "description" : "Rack definition",
                "editable" : True,
                "order" : 4
                },
            'gear_width': {
                "type": "number",
                "step" : 0.001,
                "examples": [0.02],
                "editable": True,
                "description": "Gear width",
                "order" : 5,
                "unit" : "m",
                },
             }
         }

    _display_angular = []

    # TODO: rename variables in init accordingly to their attribute name
    # Ex: cp -> coefficient_profile_shift
    def __init__(self, z, db, coefficient_profile_shift, rack, material=None,
                 gear_width=1, name=''):

        self.rack = rack
        self.GearParam(z, db, coefficient_profile_shift)

        # Definition of default parameters
        if material is None:
            self.material = hardened_alloy_steel
        self.material = material

        self.gear_width = gear_width

        DessiaObject.__init__(self, name=name)

    def __eq__(self, other_eb):
        equal = (self.z == other_eb.z
                 and self.db == other_eb.db
                 and self.coefficient_profile_shift == other_eb.coefficient_profile_shift
                 and self.rack == other_eb.rack
                 and self.material == other_eb.material)
        return equal

    def __hash__(self):
        rack_hash = hash(self.z) + hash(int(self.db*1e3))\
                 + hash(int(self.coefficient_profile_shift*1e3)) + hash(self.rack)
        return rack_hash

    def Update(self, z, db, coefficient_profile_shift, transverse_pressure_angle_rack,
               coeff_gear_addendum, coeff_gear_dedendum, coeff_root_radius,
               coeff_circular_tooth_thickness, material, gear_width=1):
        """ Update of the gear mesh

        :param all: same parameters of this class initialisation

        >>> input={z:14, db:42*1e-3, cp:0.5}
        >>> mesh1.Update(**input)
        """
        self.GearParam(z, db, coefficient_profile_shift)
        self.rack.Update(self.rack.module, transverse_pressure_angle_rack, coeff_gear_addendum,
                         coeff_gear_dedendum, coeff_root_radius,
                         coeff_circular_tooth_thickness)

        self.gear_width = gear_width

    ### geometry definition

    def GearParam(self,z, db, coefficient_profile_shift):

        self.z = z
        self.db = db
        self.dff = self.db/math.cos(self.rack.transverse_pressure_angle)
        module_rack = self.dff/self.z
        self.rack.Update(module_rack)
        self.coefficient_profile_shift = coefficient_profile_shift

        self.outside_diameter=(self.dff
                               +2*(self.rack.gear_addendum
                                   +self.rack.module*self.coefficient_profile_shift))
        self.alpha_outside_diameter = math.acos(self.db/self.outside_diameter)
        self.root_diameter=(self.dff
                            - 2*(self.rack.gear_dedendum
                                - self.rack.module*self.coefficient_profile_shift))
        self.root_diameter_active,self.phi_trochoide=self._RootDiameterActive()
        self.alpha_root_diameter_active=math.acos(self.db/self.root_diameter_active)

        self.alpha_pitch_diameter=math.acos(self.db/self.dff)
        self.circular_tooth_thickness = (self.rack.circular_tooth_thickness
                                       +self.rack.module*self.coefficient_profile_shift
                                           *math.tan(self.rack.transverse_pressure_angle)
                                       +self.rack.module*self.coefficient_profile_shift
                                           *math.tan(self.rack.transverse_pressure_angle))
        self.tooth_space=self.rack.transverse_radial_pitch-self.circular_tooth_thickness
        self.outside_active_angle = (2*self.circular_tooth_thickness/self.dff-2
                                     *(math.tan(self.alpha_outside_diameter)
                                        -self.alpha_outside_diameter
                                        -math.tan(self.alpha_pitch_diameter)
                                        +self.alpha_pitch_diameter))
        self.base_circular_tooth_thickness = (self.db/2
                                              *(2*self.circular_tooth_thickness/self.dff
                                                +2*(math.tan(self.alpha_pitch_diameter)
                                                -self.alpha_pitch_diameter)))

        self.root_angle=self.tooth_space/(self.dff/2)-2*(math.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter)
        self.root_gear_angle=self.circular_tooth_thickness/(self.dff/2)+2*(math.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter)

    def GearSection(self,diameter):
        """ Definition of the gear section

        :param diameter: diameter of the gear section calculation
        :type diameter: m

        :results: gear section in m

        >>> gs=mesh1.GearSection(44*1e-3)
        """
        alpha_diameter=math.acos(self.db/diameter)
        theta1=(math.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter)-(math.tan(alpha_diameter)-alpha_diameter)
        return diameter/2*(2*theta1+self.outside_active_angle)

    def _RootDiameterActive(self):
        a=self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.dff/2
        phi=-(a+b*math.tan(math.pi/2-self.rack.transverse_pressure_angle))/r
        root_diameter_active=2*norm(self._Trochoide(phi))
        return root_diameter_active,phi

    ### Optimization Method

    def ListeIneq(self):
        """ Compilation method for inequality list used by the optimizer

        :results: vector of data that should be positive
        """
        check,ineq=self.rack.CheckRackViable()
        return ineq

    ### Trace method

    def Contour(self,discret=10,list_number=[None]):
        """ Definition of the gear contour for volmdlr

        :param discret: number of discretization points on the gear mesh involute
        :param list_number: list of gear tooth to include on the graph

        :results: volmdlr profile

        >>> C1=mesh1.Contour(10)
        >>> G1=vm.Contour2D(C1)
        >>> G1.MPLPlot() # generate a plot with matplotlib
        """
        # Analytical tooth profil
        if list_number==[None]:
            list_number=npy.arange(int(self.z))
        L=[self._OutsideTrace(0)]
        L.append(self._InvoluteTrace(discret,0,'T'))
        L.append(self._TrochoideTrace(2*discret,0,'T'))
        L.append(self._RootCircleTrace(0))
        L.append(self._TrochoideTrace(2*discret,0,'R'))
        L.append(self._InvoluteTrace(discret,1,'R'))
        for i in list_number[1::]:
            L.append(self._OutsideTrace(i))
            L.append(self._InvoluteTrace(discret,i,'T'))
            L.append(self._TrochoideTrace(2*discret,i,'T'))
            L.append(self._RootCircleTrace(i))
            L.append(self._TrochoideTrace(2*discret,i,'R'))
            L.append(self._InvoluteTrace(discret,i+1,'R'))
        return L

    def _InvoluteTrace(self,discret,number,ind='T'):

        if ind=='T':
            drap=1
            theta=npy.linspace(math.tan(self.alpha_outside_diameter),
                               math.tan(self.alpha_root_diameter_active),discret)
        else:
            drap=-1
            theta=npy.linspace(math.tan(self.alpha_root_diameter_active),
                               math.tan(self.alpha_outside_diameter),discret)

        sol=self._Involute(drap*theta)
        x=sol[0]
        y=sol[1]
        p=[vm.Point2D((x[0],y[0]))]
        for i in range(1,discret):
            p.append(vm.Point2D((x[i],y[i])))
        ref=primitives2D.OpenedRoundedLineSegments2D(p,{},False)

        if ind=='T':
            L=ref.Rotation(vm.Point2D((0,0)),-number*2*math.pi/self.z)
            self.rac=L.points[-1]
        else:
            L=ref.Rotation(vm.Point2D((0,0)),
                           self.base_circular_tooth_thickness*2/self.db)
            L=L.Rotation(vm.Point2D((0,0)),-number*2*math.pi/self.z)
            L.points[0]=self.rac
        return L

    def _TrochoideTrace(self, discret, number, type_flank='T'):
        # Function evolution of the trochoide
        if type_flank=='T':
            indice_flank=1
        else:
            indice_flank=-1

        a=indice_flank*self.rack.a # indice a in the ISO definition of the rack
        phi0=a/(self.dff/2)

        list_2D=[]
        if type_flank=='R':
            theta=npy.linspace(phi0,indice_flank*self.phi_trochoide,discret)
        else:
            theta=npy.linspace(indice_flank*self.phi_trochoide,phi0,discret)
        for t in theta:
            list_2D.append(vm.Point2D((self._Trochoide(t,type_flank))))
        list_2D=primitives2D.OpenedRoundedLineSegments2D(list_2D,{},False)

        list_2D=list_2D.Rotation(vm.Point2D((0,0)),-self.root_angle/2)

        if type_flank=='T':
            export_2D=list_2D.Rotation(vm.Point2D((0,0)),-number*2*math.pi/self.z)
            export_2D.points[0]=self.rac
        else:
            export_2D=list_2D.Rotation(vm.Point2D((0,0)),-number*2*math.pi/self.z)
            self.rac=export_2D.points[-1]
        return export_2D

    def _RootCircleTrace(self,number):
        # 2D trace of the connection between the two trochoide

        # on the drive flank
        indice_flank=1
        a=indice_flank*self.rack.a
        phi0=a/(self.dff/2)
        p1=vm.Point2D((self._Trochoide(phi0,'T')))
        p1=p1.Rotation(vm.Point2D((0,0)),-self.root_angle/2)

        # on the coast flank
        indice_flank=-1
        a=indice_flank*self.rack.a
        phi0=a/(self.dff/2)
        p2=vm.Point2D((self._Trochoide(phi0,'R')))
        p2=p2.Rotation(vm.Point2D((0,0)),-self.root_angle/2)

        list_2D=primitives2D.OpenedRoundedLineSegments2D([p1,p2],{},False)

        export_2D=list_2D.Rotation(vm.Point2D((0,0)),-number*2*math.pi/self.z)
        return export_2D

    def _OutsideTrace(self,number):
        # Trace of the top of the gear mesh
        theta4=math.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter
        p1=vm.Point2D((self.outside_diameter/2*math.cos(theta4),self.outside_diameter/2*math.sin(theta4)))
        p2=p1.Rotation(vm.Point2D((0,0)),self.outside_active_angle/2)
        p3=p2.Rotation(vm.Point2D((0,0)),self.outside_active_angle/2)
        list_2D=primitives2D.OpenedRoundedLineSegments2D([p3,p2,p1],{},False)

        export_2D=list_2D.Rotation(vm.Point2D((0,0)),-number*2*math.pi/self.z)
        return export_2D

    def _Involute(self,tan_alpha):
        """ Involute function estimation

        :param tan_alpha: tan of the pressure angle
        :results: tuple of the involute point (x,y) with the origin of the involute position at the point (x=base_diameter/2, y=0) and pressure angle is positive in the direction counter clockwise
        """
        x, y = [], []
        for ta in tan_alpha:
            x.append(self.db/2*math.cos(ta)+self.db/2*ta*math.sin(ta))
            y.append(self.db/2*math.sin(ta)-self.db/2*ta*math.cos(ta))
        return (x, y)

    def _Trochoide(self,phi,type_flank='T'):
        # function generation of trochoide point
        if type_flank=='T':
            indice_flank=1
        else:
            indice_flank=-1
        a=indice_flank*self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.dff/2
        rho=self.rack.root_radius
        x2=rho*math.sin(math.atan((a-r*phi)/b)-phi)+a*math.cos(phi)-b*math.sin(phi)+r*(math.sin(phi)-phi*math.cos(phi))
        y2=-rho*math.cos(math.atan((a-r*phi)/b)-phi)-a*math.sin(phi)-b*math.cos(phi)+r*(math.cos(phi)+phi*math.sin(phi))
        export_point=(y2,x2)
        return export_point

    ### Method for ISO stress

    def _ISO_YS(self,s_thickness_iso):
        # stress concentration factor for ISO approach
        rho_f=self.rack.root_radius+self.rack.b**2/(self.dff/2+self.rack.b)
        coeff_ys_iso=1+0.15*s_thickness_iso/rho_f
        return coeff_ys_iso

    def GearISOSection(self,angle):
        """ Calculation of the ISO section

        :param angle: pressure angle of the ISO section calculation
        :type angle: radian

        :results: ISO section and ISO height
        """
        a=self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.dff/2
        theta0 = fsolve((lambda theta:a + b*math.tan(theta) + r*(-angle - self.root_angle/2 - theta + math.pi/2)) ,0)[0]
        phi0=(a-b*math.tan(theta0))/r
        pt_iso = self._Trochoide(phi0)
        angle0 = math.atan(pt_iso[1]/pt_iso[0])-self.root_angle/2
        angle_iso = self.root_gear_angle-2*angle0
        diameter_iso = 2*norm(pt_iso)
        s_thickness_iso = diameter_iso*math.sin(angle_iso/2)
        h_height_iso = (s_thickness_iso/2)/math.tan(angle)
        return s_thickness_iso, h_height_iso

    ### Export and graph method

    def PlotData(self, x, heights, y, z, labels = True):
        transversal_plot_data = []
        axial_plot_data = []

#        points = []
#        for p in self.Contour(2):
#            for point in p.points:
#                points.extend(point.vector)
#        transversal_plot_data.append({'type' : 'line',
#                              'data' : points,
#                              'color' : [0, 0, 0],
#                              'dash' : 'none',
#                              'stroke_width' : 2,
#                              'marker' : '.',
#                              'size' : 1})
        # Outer diameter
        transversal_plot_data.append({'type' : 'circle',
                          'cx' : y,
                          'cy' : z,
                          'r' : 0.5 * self.outside_diameter,
                          'color' : [0, 0, 0],
                          'size' : 1,
                          'group' : 3,
                          'dash' : 'none',})

        return transversal_plot_data, axial_plot_data

    def Dict(self):
        d = {'name' : self.name} # TODO Change this to DessiaObject.__init__
        d['z'] = self.z
        d['db'] = self.db
        d['coefficient_profile_shift'] = self.coefficient_profile_shift
        d['gear_width'] = self.gear_width
        d['rack'] = self.rack.Dict()
        d['material'] = self.material.Dict()
        return d

    @classmethod
    def DictToObject(cls, d):
        material = Material.DictToObject(d['material'])
        rack = Rack.DictToObject(d['rack'])
        mesh = cls(z = d['z'],
                   db = d['db'],
                   coefficient_profile_shift = d['coefficient_profile_shift'],
                   rack = rack,
                   gear_width = d['gear_width'],
                   material = material,
                   name=d['name'])
        return mesh

class MeshCombination(DessiaObject):
    def __init__(self,  center_distance, connections, meshes,
                 torque=None, cycle=None, safety_factor=1, name=''):

        self.center_distance = center_distance
        self.connections = connections
        self.meshes = meshes
        self.torque = torque
        self.cycle = cycle
        self.safety_factor = safety_factor

        self.minimum_gear_width = 10e-3
        self.helix_angle = 0

        # NetworkX graph construction
        list_gear = []
        for gs in self.connections:
            for g in gs:
                if g not in list_gear:
                    list_gear.append(g)
        gear_graph = nx.Graph()
        gear_graph.add_nodes_from(list_gear)
        gear_graph.add_edges_from(self.connections)
        self.gear_graph = gear_graph
        self.list_gear = list_gear

        self.transverse_pressure_angle = []
        for num_gear, (num1, num2) in enumerate(self.connections):
            mesh_first = self.meshes[num1]
            mesh_second = self.meshes[num2]
            df_first = 2*self.center_distance[num_gear]*mesh_first.z/mesh_second.z/(1+mesh_first.z/mesh_second.z)
            self.transverse_pressure_angle.append(math.acos(mesh_first.db/df_first))
        transverse_pressure_angle_0 = self.transverse_pressure_angle[0]

        self.Z = {i: mesh.z for i, mesh in meshes.items()}
        self.material = {i: mesh.material for i, mesh in meshes.items()}
        self.gear_width = {}
        self.DB = {}
        for num_mesh, mesh in self.meshes.items():
            self.gear_width[num_mesh] = mesh.gear_width
            self.DB[num_mesh] = mesh.db

        self.DF, DB_new, self.connections_dfs, transverse_pressure_angle_new\
            = MeshCombination.GearGeometryParameter(self.Z, transverse_pressure_angle_0, center_distance,
                                                    connections, gear_graph)
        if len(cycle.keys())<len(list_gear): # the gear mesh optimizer calculate this dictionary
            self.cycle = MeshCombination.CycleParameter(cycle, self.Z, list_gear)
        self.torque, self.normal_load, self.tangential_load, self.radial_load = MeshCombination.GearTorque(self.Z, torque, self.DB,
                    gear_graph, list_gear, connections, self.DF, self.transverse_pressure_angle)
        self.linear_backlash, self.radial_contact_ratio = \
            MeshCombination.GearContactRatioParameter(self.Z, self.DF, self.transverse_pressure_angle,
                                                      center_distance,
                                                      meshes, self.connections_dfs, connections)

        gear_width_new, self.sigma_iso, self.sigma_lim = MeshCombination.GearWidthDefinition(self.safety_factor,
                                            self.minimum_gear_width,
                                            list_gear, self.tangential_load, meshes,
                                            connections,
                                            self.material, self.cycle, self.radial_contact_ratio, self.helix_angle,
                                            self.transverse_pressure_angle)
        self.check()

        DessiaObject.__init__(self, name=name)

    def __eq__(self, other_eb):
        equal = (self.center_distance == other_eb.center_distance
                 and self.connections == other_eb.connections
                 and self.meshes == other_eb.meshes
                 and self.torque == other_eb.torque
                 and self.cycle == other_eb.cycle
                 and self.safety_factor == other_eb.safety_factor)
        return equal

    def __hash__(self):
        mc_hash = 0
        for num_mesh, mesh in self.meshes.items():
            mc_hash += hash(mesh)
        for center_distance in self.center_distance:
            mc_hash += hash(int(center_distance*1e3))
        return mc_hash

    def check(self):
        valid = True
        gear_width, _, _ = MeshCombination.GearWidthDefinition(self.safety_factor,
                                            self.minimum_gear_width,
                                            self.list_gear, self.tangential_load, self.meshes,
                                            self.connections,
                                            self.material, self.cycle, self.radial_contact_ratio,
                                            self.helix_angle,
                                            self.transverse_pressure_angle)
        for num_mesh, mesh in self.meshes.items():
            if abs(gear_width[num_mesh] - mesh.gear_width) > 1e-6:
                valid = False
        self.DF, DB_new, self.connections_dfs, transverse_pressure_angle_new\
            = MeshCombination.GearGeometryParameter(self.Z, self.transverse_pressure_angle[0], self.center_distance,
                                                    self.connections, self.gear_graph)
        for num_mesh, mesh in self.meshes.items():
            if abs(DB_new[num_mesh] - mesh.db) > 1e-6:
                valid = False
        for num_gear, connection in enumerate(self.connections):
            if abs(transverse_pressure_angle_new[num_gear] - self.transverse_pressure_angle[num_gear]) > 1e-6:
                valid = False
        return valid

    @classmethod
    def create(cls, Z, center_distance, connections, transverse_pressure_angle_0,
               coefficient_profile_shift, transverse_pressure_angle_rack,
               coeff_gear_addendum, coeff_gear_dedendum, coeff_root_radius,
               coeff_circular_tooth_thickness, material=None, torque=None, cycle=None,
               safety_factor=1):

        # NetworkX graph construction
        list_gear = []
        for gs in connections:
            for g in gs:
                if g not in list_gear:
                    list_gear.append(g)
        gear_graph = nx.Graph()
        gear_graph.add_nodes_from(list_gear)
        gear_graph.add_edges_from(connections)

        # Definition of default parameters
        minimum_gear_width = 10e-3
        helix_angle = 0

        if material == None:
            material = {list_gear[0]:hardened_alloy_steel}
        for ne in list_gear:
            if ne not in material.keys():
                material[ne] = hardened_alloy_steel

        if torque == None:
            torque = [{list_gear[0]:100,list_gear[1]:'output'}]

        if cycle == None:
            cycle = {list_gear[0]:1e6}

        DF, DB, connections_dfs, transverse_pressure_angle\
            = cls.GearGeometryParameter(Z, transverse_pressure_angle_0, center_distance,
                                         connections, gear_graph)

        if len(cycle.keys())<len(list_gear): # the gear mesh optimizer calculate this dictionary
            cycle = cls.CycleParameter(cycle, Z, list_gear)

        torque, normal_load, tangential_load, radial_load = cls.GearTorque(Z, torque, DB,
                    gear_graph, list_gear, connections, DF, transverse_pressure_angle)

        meshes={}
        for num_engr in list_gear:
            z = Z[num_engr]
            db = DB[num_engr]
            cp = coefficient_profile_shift[num_engr]
#            ngp=self.list_gear.index(num_engr)
            tpa = transverse_pressure_angle_rack[num_engr]
            cga = coeff_gear_addendum[num_engr]
            cgd = coeff_gear_dedendum[num_engr]
            crr = coeff_root_radius[num_engr]
            cct = coeff_circular_tooth_thickness[num_engr]
            mat = material[num_engr]
            rack = Rack(transverse_pressure_angle = tpa,
                        coeff_gear_addendum = cga, coeff_gear_dedendum = cgd,
                        coeff_root_radius = crr, coeff_circular_tooth_thickness = cct)
            meshes[num_engr] = Mesh(z, db, cp, rack, mat)

        linear_backlash, radial_contact_ratio = \
            cls.GearContactRatioParameter(Z, DF, transverse_pressure_angle,
                                          center_distance,
                                          meshes, connections_dfs, connections)

        gear_width, sigma_iso, sigma_lim = cls.GearWidthDefinition(safety_factor,
                                            minimum_gear_width,
                                            list_gear, tangential_load, meshes,
                                            connections,
                                            material, cycle, radial_contact_ratio, helix_angle,
                                            transverse_pressure_angle)

        for num_gear in list_gear:
            meshes[num_gear].gear_width = gear_width[num_gear]
        mesh_combination = cls(center_distance, connections, meshes, torque, cycle)
        return mesh_combination

    def Update(self, Z, center_distance, connections, transverse_pressure_angle_0,
               coefficient_profile_shift,
               transverse_pressure_angle_rack, coeff_gear_addendum,
               coeff_gear_dedendum, coeff_root_radius, coeff_circular_tooth_thickness,
               material, torque, cycle, safety_factor):
        """ Update of the gear mesh assembly

        :param all: same parameters of this class initialisation

        >>> Z={1:13,2:46,4:38}
        >>> center_distance=[0.118,0.125]
        >>> mesh_assembly1.Update(Z=Z,center_distance=center_distance)
        """
        self.center_distance = center_distance
        self.transverse_pressure_angle_0 = transverse_pressure_angle_0
        self.DF, self.DB, self.connections_dfs, self.transverse_pressure_angle\
            = MeshCombination.GearGeometryParameter(Z, transverse_pressure_angle_0, center_distance,
                                         connections, self.gear_graph)
        for num_engr in self.list_gear:
            z = Z[num_engr]
            db = self.DB[num_engr]
            cp = coefficient_profile_shift[num_engr]
            tpa = transverse_pressure_angle_rack[num_engr]
            cga = coeff_gear_addendum[num_engr]
            cgd = coeff_gear_dedendum[num_engr]
            crr = coeff_root_radius[num_engr]
            cct = coeff_circular_tooth_thickness[num_engr]
            mat = self.material[num_engr]
            self.meshes[num_engr].Update(z, db, cp, tpa, cga, cgd,
                                         crr, cct, mat)
        self.linear_backlash, self.radial_contact_ratio = \
            MeshCombination.GearContactRatioParameter(Z, self.DF, self.transverse_pressure_angle,
                                                      center_distance,
                                                      self.meshes, self.connections_dfs, connections)

    ### Optimization Method
    def CheckMinimumBacklash(self,backlash_min=2*1e-4):
        """ Define constraint and functional for the optimizer on backlash

        :param backlash_min: maximum backlash available
        :results:
            * check is a boolean (True if 0<backlash<backlash_min)
            * list_ineq a list of element that should be positive for the optimizer
            * obj is a functional on the backlash used for the optimizer
        """
        list_ineq=[] # liste of value to evaluate backlash
        obj=0
        for lb in self.linear_backlash:
            list_ineq.append(lb) # backlash > 0
            list_ineq.append(backlash_min-lb) # backlash < backlash_min so (backlash_min-backlash)>0
            obj+=10*(lb-backlash_min)**2
        check=False
        if min(list_ineq)>0:
            check=True
        return check,list_ineq,obj

    def CheckRadialContactRatio(self,radial_contact_ratio_min=1):
        """ Define constraint and functional for the optimizer on radial contact ratio

        :param radial_contact_ratio_min: minimum radial contact ratio available
        :results:
            * check is a boolean (True if radial_contact_ratio_min<radial_contact_ratio)
            * list_ineq a list of element that should be positive for the optimizer
            * obj is a functional on the backlash used for the optimizer
        """
        list_ineq=[]
        obj=0
        for num_mesh,(eng1,eng2) in enumerate(self.connections):
            rca=self.radial_contact_ratio[num_mesh]
            list_ineq.append(rca-radial_contact_ratio_min)
            if rca>radial_contact_ratio_min:
                obj+=0.001*(rca-radial_contact_ratio_min)
            else:
                obj+=1000*(radial_contact_ratio_min-rca)
        check=False
        if min(list_ineq)>0:
            check=True
        return check,list_ineq,obj

    def ListeIneq(self):
        """ Compilation method for inequality list used by the optimizer

        :results: vector of data that should be positive
        """
        _,ineq,_=self.CheckMinimumBacklash(4*1e-4)
        _,list_ineq,_=self.CheckRadialContactRatio(1)
        ineq.extend(list_ineq)

        for num_gear,mesh in self.meshes.items():
            list_ineq=mesh.ListeIneq()
            ineq.extend(list_ineq)

        return ineq

    def Functional(self):
        """ Compilation method for a part of the functional used by the optimizer

        :results: scalar add to the global functional of the optimizer
        """
        check1,ineq1,obj1 = self.CheckMinimumBacklash(4*1e-4)
        check2,ineq2,obj2 = self.CheckRadialContactRatio(1)
        obj = obj1 + obj2
        return obj

    ### Method gear mesh calculation
    @classmethod
    def GearContactRatioParameter(cls, Z, DF, transverse_pressure_angle,
                           center_distance,
                           meshes, connections_dfs, connections):


        linear_backlash = []
        radial_contact_ratio = []
        for engr1,engr2 in connections_dfs:
            if (engr1,engr2) in connections:
                num_mesh = connections.index((engr1,engr2))
            elif (engr2,engr1) in connections:
                num_mesh = connections.index((engr2,engr1))
            else:
                raise RuntimeError
            circular_tooth_thickness1 = meshes[engr1].GearSection(DF[num_mesh][engr1])
            circular_tooth_thickness2 = meshes[engr2].GearSection(DF[num_mesh][engr2])
            transverse_radial_pitch1=math.pi*DF[num_mesh][engr1]/meshes[engr1].z
            space_width1=transverse_radial_pitch1-circular_tooth_thickness1
            space_width2=transverse_radial_pitch1-circular_tooth_thickness2
            linear_backlash.append(min(space_width1-circular_tooth_thickness2,space_width2-circular_tooth_thickness1))
            transverse_pressure_angle1 = transverse_pressure_angle[num_mesh]
            center_distance1 = center_distance[num_mesh]
            radial_contact_ratio.append((1/2.*(math.sqrt(meshes[engr1].outside_diameter**2
                                                       - meshes[engr1].db**2)
                                              + math.sqrt(meshes[engr2].outside_diameter**2
                                                         - meshes[engr2].db**2)
                                        - 2*center_distance1*math.sin(transverse_pressure_angle1))
                                        /(transverse_radial_pitch1*math.cos(transverse_pressure_angle1))))
        return linear_backlash,radial_contact_ratio

    @classmethod
    def GearGeometryParameter(cls, Z, transverse_pressure_angle_0, center_distance, connections, gear_graph):
        # Construction of pitch and base diameter
        DF = {}
        db = {}
        dict_transverse_pressure_angle = {0: transverse_pressure_angle_0}
        connections_dfs = list(nx.edge_dfs(gear_graph,
                            [connections[0][0], connections[0][1]]))
        for num_dfs,((engr1,engr2),cd) in enumerate(zip(connections_dfs, center_distance)):
            if (engr1,engr2) in connections:
                num_mesh = connections.index((engr1,engr2))
            else:
                num_mesh = connections.index((engr2,engr1))
            Z1 = Z[engr1]
            Z2 = Z[engr2]
            DF1 = 2*cd*Z1/Z2/(1+Z1/Z2)
            DF2 = 2*cd-DF1
            DF[num_mesh] = {}
            DF[num_mesh][engr1] = DF1
            DF[num_mesh][engr2] = DF2
            if num_mesh == 0:
                db1 = float(DF1*math.cos(transverse_pressure_angle_0))
                db2 = float(DF2*math.cos(transverse_pressure_angle_0))
            else:
                db1 = db[engr1]
                try:
                    dict_transverse_pressure_angle[num_mesh] = math.acos(db1/DF1)
                except:
                    print('Error Diameter DB {}, DF {}, Z1 {}, Z2 {}, pa {}'.format(db1, DF1, Z1, Z2, transverse_pressure_angle_0))
                    raise ValidGearDiameterError()
                db2 = DF2*math.cos(dict_transverse_pressure_angle[num_mesh])
            db[engr1] = db1
            db[engr2] = db2
        transverse_pressure_angle = []
        for num_mesh in sorted(dict_transverse_pressure_angle.keys()):
            tpa = dict_transverse_pressure_angle[num_mesh]
            transverse_pressure_angle.append(tpa)

        return DF, db, connections_dfs, transverse_pressure_angle

    @classmethod
    def GearTorque(cls, Z, torque, db, gear_graph, list_gear, connections, DF, transverse_pressure_angle):
        """ Calculation of the gear mesh torque

        :param Z: dictionary define the number of teeth {node1:Z1, node2:Z2, mesh3:Z3 ...}
        :param torque: dictionary defining all input torque, one node where the torque is not specified is define as the 'output' {node1:torque1, node2:torque2, node3:'output'}
        :param db: dictionary define the base diameter {mesh1: {node1:db1_a, node2:db2_a}, mesh2: {node2:db2_b, node3:db3_b}}
        :type db: m

        :results:
            * **torque1** - dictionary of the applied torque on gear mesh (torque applied by node_x on node_y) {node1:tq1, node3:tq3 ...}
            * **torque2** - dictionary of the counter drive torque on gear mesh (torque applied by node_x on node_y) {node2:-tq1, node4:-tq3 ...}
            * **normal_load** - dictionary define the normal load for each gear mesh (applied torque define the direction) {mesh1 : [fn_x1,fn_y1,fn_z1],mesh2 : [fn_x2,fn_y2,fn_z2] ...}
            * **tangential_load** - dictionary define the tangential load for each gear mesh (applied torque define the direction) {mesh1 : [ft_x1,ft_y1,ft_z1],mesh2 : [ft_x2,ft_y2,ft_z2] ...}
            * **radial_load** - dictionary define the radial load for each gear mesh (applied torque define the direction) {mesh1 : [fr_x1,fr_y1,fr_z1],mesh2 : [fr_x2,fr_y2,fr_z2] ...}

        be careful, due to the parameters of the gear mesh assembly (define one pressure angle for each mesh) the diameter db2_a is different to db2_b (you have to define correctly transverse_pressure_angle to have db2_a=db2_b)
        """
        if 'output' in torque.values():
            for num_gear,tq in torque.items():
                if tq=='output':
                    node_output=num_gear
            torque_graph_dfs=list(nx.dfs_edges(gear_graph,node_output))
            order_torque_calculation=[(eng2,eng1) for (eng1,eng2) in torque_graph_dfs[::-1]]
            # calculation torque distribution
            temp_torque={}
            for eng1 in list_gear:
                temp_torque[eng1]=0
            for num_mesh_tq,(eng1,eng2) in enumerate(order_torque_calculation):
                if eng1 in torque.keys():
                    temp_torque[eng1]+=torque[eng1]
                temp_torque[eng2]+=-temp_torque[eng1]*Z[eng2]/float(Z[eng1])
            dic_torque={}
            for num_mesh_tq,(eng1,eng2) in enumerate(order_torque_calculation):
                dic_torque[(eng1,eng2)]=temp_torque[eng1]

        normal_load={}
        tangential_load={}
        radial_load={}

        for num_mesh,(eng1,eng2) in enumerate(connections):
            if 'output' not in torque.values():
                dic_torque=torque
            try:
                tq=dic_torque[(eng1,eng2)]
            except:
                tq=dic_torque[(eng2,eng1)]
            normal_load[num_mesh]=abs(tq)*2/(db[eng1])
            tangential_load[num_mesh]=abs(tq)*2/(DF[num_mesh][eng1])
            radial_load[num_mesh]=math.tan(transverse_pressure_angle[num_mesh])*tangential_load[num_mesh]
        return dic_torque, normal_load, tangential_load, radial_load

    @classmethod
    def CycleParameter(cls, cycle, Z, list_gear):
        """ Calculation of the gear mesh cycle

        :param Z: dictionary define the number of teeth {node1:Z1, node2:Z2, node3:Z3 ...}
        :param cycle: Dictionary defining the number of cycle for one node {node3: number_cycle3}

        :results: dictionary define the number of cycle for each gear mesh {node1:cycle1, node2:cycle2, node3:cycle3 ...}
        """
        eng_init=list(cycle.keys())[0]
        for eng in list_gear:
            if eng not in cycle.keys():
                cycle[eng]=cycle[eng_init]*Z[eng_init]/float(Z[eng])
        return cycle

    @classmethod
    def GearWidthDefinition(cls, safety_factor, minimum_gear_width,
                            list_gear, tangential_load, meshes, connections,
                            material, cycle, radial_contact_ratio, helix_angle,
                            transverse_pressure_angle):
        """ Calculation of the gear width

        :param safety_factor: Safety factor used for the ISO design

        :results:
            * **gear_width** - dictionary define the gear mesh width {node1 : gw1, node2 : gw2, node3 : gw3 ...}
            * **sigma_iso** - dictionary define the ISO stress {mesh1 : {node1 sig_iso1: , node2 : sig_iso2_1}, mesh2 : {node2 : sig_iso2_2, node3 : sig_iso3} ...}
            * **sigma_lim** - dictionary define the limit material stress {mesh1 : {node1 sig_lim1: , node2 : sig_lim2}, mesh2 : {node2 : sig_lim2, node3 : sig_lim3} ...}

        in this function, we define the gear width for each gear mesh to respect sig_lim = sig_iso for each gear mesh
        """
        coeff_yf_iso = cls._CoeffYFIso(connections, meshes, transverse_pressure_angle)
        coeff_ye_iso = cls._CoeffYEIso(connections, radial_contact_ratio)
        coeff_yb_iso = cls._CoeffYBIso(connections, material, helix_angle)

        sigma_lim = cls.SigmaMaterialISO(safety_factor, connections,
                                          material, cycle, meshes)
        gear_width = {}
        for eng in list_gear:
            gear_width[eng] = minimum_gear_width

        for num_mesh,(eng1,eng2) in enumerate(connections):
            gear_width1 = abs(tangential_load[num_mesh]
                        / (sigma_lim[num_mesh][eng1]
                        * meshes[eng1].rack.module)
                        *coeff_yf_iso[num_mesh][eng1]
                        *coeff_ye_iso[num_mesh]
                        *coeff_yb_iso[num_mesh][eng1])

            gear_width2=abs(tangential_load[num_mesh]
                            /(sigma_lim[num_mesh][eng2]
                                *meshes[eng2].rack.module)
                            *coeff_yf_iso[num_mesh][eng2]
                            *coeff_ye_iso[num_mesh]
                            *coeff_yb_iso[num_mesh][eng2])

            gear_width_set=max(gear_width1,gear_width2)
            gear_width[eng1]=max(gear_width[eng1],gear_width_set)
            gear_width[eng2]=max(gear_width[eng2],gear_width_set)
        sigma_iso=sigma_lim
        return gear_width, sigma_iso, sigma_lim

    @classmethod
    def SigmaMaterialISO(cls, safety_factor, connections, material, cycle,
                         meshes):
        """ Calculation of the material limit stress

        :param safety_factor: Safety factor used for the ISO design

        :results:
            * **sigma_lim** - dictionary define the limit material stress {mesh1 : {node1 sig_lim1: , node2 : sig_lim2}, mesh2 : {node2 : sig_lim2, node3 : sig_lim3} ...}

        in this function, we use the FunCoeff function of the Material class to interpolate the material parameters
        """
        angle=30./180.*math.pi
        sigma_lim={}
        for num_mesh,(eng1,eng2) in enumerate(connections):
            sigma_lim[num_mesh] = {}

            matrice_wholer = material[eng1].data_wholer_curve
            matrice_material = material[eng1].data_gear_material
            sgla = material[eng1].FunCoeff(cycle[eng1],npy.array(matrice_wholer['data']),matrice_wholer['x'],matrice_wholer['y'])
            sgl1 = material[eng1].FunCoeff(sgla,npy.array(matrice_material['data']),matrice_material['x'],matrice_material['y'])
            s_thickness_iso_1,h_height_iso_1 = meshes[eng1].GearISOSection(angle)
            coeff_ys_iso = meshes[eng1]._ISO_YS(s_thickness_iso_1)
            sigma_lim[num_mesh][eng1] = float((sgl1/(safety_factor*coeff_ys_iso))*10**7)

            matrice_wholer = material[eng2].data_wholer_curve
            matrice_material = material[eng2].data_gear_material
            sglb = material[eng2].FunCoeff(cycle[eng2], npy.array(matrice_wholer['data']), matrice_wholer['x'],matrice_wholer['y'])
            sgl2 = material[eng2].FunCoeff(sglb, npy.array(matrice_material['data']),matrice_material['x'],matrice_material['y'])
            s_thickness_iso_2,h_height_iso_2 = meshes[eng2].GearISOSection(angle)
            coeff_ys_iso = meshes[eng2]._ISO_YS(s_thickness_iso_2)
            sigma_lim[num_mesh][eng2] = float((sgl2/(safety_factor*coeff_ys_iso))*10**7)
        return sigma_lim

    @classmethod
    def _CoeffYFIso(cls, connections, meshes, transverse_pressure_angle):
        # shape factor for ISO stress calculation
        angle=30./180.*math.pi
        coeff_yf_iso={}
        for num_mesh,(eng1,eng2) in enumerate(connections):
            coeff_yf_iso[num_mesh]={}
            s_thickness_iso_1,h_height_iso_1 = meshes[eng1].GearISOSection(angle)
            s_thickness_iso_2,h_height_iso_2 = meshes[eng2].GearISOSection(angle)
            coeff_yf_iso[num_mesh][eng1]=((6*(h_height_iso_1/meshes[eng1].rack.module)*math.cos(transverse_pressure_angle[num_mesh]))
                                    /((s_thickness_iso_1/meshes[eng1].rack.module)**2
                                       *math.cos(meshes[eng1].rack.transverse_pressure_angle)))
            coeff_yf_iso[num_mesh][eng2]=((6*(h_height_iso_2/meshes[eng2].rack.module)*math.cos(transverse_pressure_angle[num_mesh]))
                                    /((s_thickness_iso_2/meshes[eng2].rack.module)**2
                                       *math.cos(meshes[eng2].rack.transverse_pressure_angle)))
        return coeff_yf_iso

    @classmethod
    def _CoeffYEIso(cls, connections, radial_contact_ratio):
        #  radial contact ratio factor for ISO stress calculation
        coeff_ye_iso=[]
        for ne,eng in enumerate(connections):
            coeff_ye_iso.append(1/radial_contact_ratio[ne])
        return coeff_ye_iso

    @classmethod
    def _CoeffYBIso(cls, connections, material, helix_angle):
        # gear widht factor impact for ISO stress calculation
        coeff_yb_iso={}
        for num_mesh,(eng1,eng2) in enumerate(connections):
            coeff_yb_iso[num_mesh] = {}
            matrice_YB = material[eng1].data_coeff_YB_Iso
            coeff_yb_iso[num_mesh][eng1] = material[eng1].FunCoeff(helix_angle,npy.array(matrice_YB['data']),matrice_YB['x'],matrice_YB['y'])
            matrice_YB = material[eng2].data_coeff_YB_Iso
            coeff_yb_iso[num_mesh][eng2] = material[eng2].FunCoeff(helix_angle,npy.array(matrice_YB['data']),matrice_YB['x'],matrice_YB['y'])
        return coeff_yb_iso

    ### Function graph and export

    def GearRotate(self, list_gear, list_center, list_rot):
        """ Displacement of the volmdlr gear profile (rotation and translation)

        :param list_gear: list of volmdlr contour [meshes.Contour, meshes.Contour ...], each contour is centered on the origin
        :param list_center: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
        :param list_rot: list of rotation for each gear mesh [node1 : rot1, node2 : rot2 ...]

        :results: list of volmdlr component
        """
        export=[]
        for (i,center,k) in zip(list_gear,list_center,list_rot):
            model_export=[]
            for m in i:
                center = vm.Point2D(center)
                model_trans = m.Translation(center)
                model_trans_rot = model_trans.Rotation(center, k)
                model_export.append(model_trans_rot)
            export.append(model_export)
        return export

    def InitialPosition(self,set_pos,liste_eng=()):
        """ Calculation of the rotation for two gear mesh to initiate the contact

        :param list_gear: list of volmdlr contour [meshes.Contour, meshes.Contour ...], each contour is centered on the origin
        :param list_center: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
        :param list_rot: list of rotation for each gear mesh [node1 : rot1, node2 : rot2 ...]

        :results: list of volmdlr component
        """
        Angle1=math.acos(self.meshes[liste_eng[0]].db/self.DF[set_pos][liste_eng[0]])
        Angle2=math.acos(self.meshes[liste_eng[1]].db/self.DF[set_pos][liste_eng[1]])
        Gear1Angle=-(math.tan(Angle1)-Angle1)
        Gear2Angle=-(math.tan(Angle2)-Angle2)+math.pi
        return [Gear1Angle,Gear2Angle]

    # TODO: use volmdlr Vector and points
    def VolumeModel(self, centers = {}, axis = (1,0,0), name = ''):
        """ Generation of the 3D volume for all the gear mesh

        :param center: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
        :param axis: direction of gear mesh rotation

        :results: list of 3D volmdlr component
        """
        x = vm.Vector3D(axis)
        y = x.RandomUnitNormalVector()
        z = vm.Vector3D(npy.cross(x.vector, y.vector))
        if len(centers)==0:
            centers = {}
            center_var = self.PosAxis({self.list_gear[0]:[0,0]})
            for engr_num in center_var.keys():
                centers[engr_num]=tuple(center_var[engr_num][0]*y.vector+center_var[engr_num][1]*z.vector)
        else:
            center_var={}
            for engr_num in centers.keys():
                center_var[engr_num]=(npy.dot(centers[engr_num],x.vector),npy.dot(centers[engr_num],y.vector),npy.dot(centers[engr_num],z.vector))
            centers=center_var

        Gears3D={}
        Struct=[]
        Rotation={}
        primitives=[]

        for set_pos_dfs,(eng1,eng2) in enumerate(self.connections_dfs):
            position1 = centers[eng1]
            position2 = centers[eng2]

            if (eng1,eng2) in self.connections:
                set_pos=self.connections.index((eng1,eng2))
                list_rot=self.InitialPosition(set_pos,(eng1,eng2))
            elif (eng2,eng1) in self.connections:
                set_pos=self.connections.index((eng2,eng1))
                list_rot=self.InitialPosition(set_pos,(eng2,eng1))
            Rotation[set_pos]={}
            if set_pos_dfs==0:
                Gears3D[eng1]=self.meshes[eng1].Contour(3)
            Struct.append(vm.Circle2D(vm.Point2D(position1),self.DF[set_pos][eng1]/2.))
            Gears3D[eng2]=self.meshes[eng2].Contour(3)
            Struct.append(vm.Circle2D(vm.Point2D(position2),self.DF[set_pos][eng2]/2.))

            if position2[1]==position1[1]:
                if position2[2]-position1[2]>0:
                    angle0=math.pi/2.
                else:
                    angle0=-math.pi/2.

            else:
                angle0=-math.atan((position2[2]-position1[2])/(position2[1]-position1[1]))
                if (position2[2]-position1[2])<0:
                    angle0=angle0+math.pi
            if set_pos_dfs==0:
                Rotation[set_pos][eng1]=list_rot[0]+angle0
                Rotation[set_pos][eng2]=list_rot[1]+angle0
            else:
                for k1,rot in Rotation.items():
                    if eng1 in rot.keys():
                        Rotation[set_pos][eng1]=rot[eng1]
                        delta_rot=Rotation[set_pos][eng1]-(list_rot[0]-angle0)
                Rotation[set_pos][eng2]=list_rot[1]-angle0-delta_rot*((self.meshes[eng1].z)/(self.meshes[eng2].z))
            Gears3D_Rotate=self.GearRotate([Gears3D[eng1],Gears3D[eng2]],[(position1[1::]),(position2[1::])],
                                       list_rot=[Rotation[set_pos][eng1],Rotation[set_pos][eng2]])

            C1=vm.Contour2D(Gears3D_Rotate[0])
            # print(Gears3D_Rotate[0])
            C2=vm.Contour2D(Gears3D_Rotate[1])

            extrusion_vector1 = (self.gear_width[eng1]*x)
            extrusion_vector2 = (self.gear_width[eng2]*x)

            if set_pos_dfs==0:
                vect_x = -0.5*self.gear_width[eng1]*x + vm.Vector3D((x.Dot(vm.Vector3D(centers[eng1])), 0,0))
                t1=primitives3D.ExtrudedProfile(vm.Vector3D(vect_x), y, z, C1, [], vm.Vector3D(extrusion_vector1))
                primitives.append(t1)
            vect_x = -0.5*self.gear_width[eng2]*x + vm.Vector3D((x.Dot(vm.Vector3D(centers[eng2])), 0,0))
            t2=primitives3D.ExtrudedProfile(vm.Vector3D(vect_x),y,z, C2, [], vm.Vector3D(extrusion_vector2))
            primitives.append(t2)

        model = vm.VolumeModel(primitives, name)
        return model

    def Mass(self):
        """
        Estimation of gear mesh mass

        :results: mass of all gear mesh
        """
        DF = {}
        for i,(ic1, ic2) in enumerate(self.connections):
            DF[ic1] = self.DF[i][ic1]
            DF[ic2] = self.DF[i][ic2]

        mass = 0.
        for i,df in DF.items():
            mass +=  self.gear_width[i] * self.material[i].volumic_mass* math.pi * (0.5*DF[i])**2
        return mass

    # Waiting for meshes to know how to plot themselves
#    def PlotData(self, x, heights, ys, zs, labels = True):
#        transversal_plot_data = []
#        axial_plot_data = []
#        # TODO remove when meshes would be a list
#        imesh = []
#        meshes = []
#        for ic, connec in enumerate(self.connections):
#            imesh.extend(connec)
#        imesh = list(set(imesh))
#        for ic, (ic1, ic2) in enumerate(self.connections):
#            if ic1 in imesh:
#                meshes.append(self.meshes[ic][ic1])
#                imesh.remove(ic1)
#            if ic2 in imesh:
#                meshes.append(self.meshes[ic][ic2])
#                imesh.remove(ic2)
#
#        for imesh, mesh in enumerate(meshes):
#            t, a = mesh.PlotData(x, heights, ys[imesh], zs[imesh], labels)
#            transversal_plot_data.extend(t)
#            axial_plot_data.extend(a)
#
#        # Ploting axial because mesh doesn't know its width
#
#        return axial_plot_data, transversal_plot_data

    def FreeCADExport(self, fcstd_filepath, centers = {}, axis = (1,0,0),
                      python_path='python', path_lib_freecad='/usr/lib/freecad/lib',
                      export_types=['fcstd']):
        """ Export 3D volume to FreeCAD

        :param file_path: file path for the freecad file
        :param center: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
        :param axis: direction of gear mesh rotation

        :results: export of a FreeCAD file
        """
        model = self.VolumeModel(centers, axis)
        model.FreeCADExport(fcstd_filepath,python_path,path_lib_freecad,export_types)

    # TODO change this function to make it like PlotData in PWT: output is dict of geometrical shapes
    def SVGExport(self,name,position):
        """ Export SVG graph of all gear mesh

        :param name: name of the svg file
        :param position: dictionary define some center position {node2 : [0,0], node4 : [0.12,0] ..}

        :results: SVG graph

        in the position dictionary, you have to be coherent with the center position

            * for exemple, if the center-distance of the mesh1 (node1, node2) is 0.117 m you can define position such as:

                * {node1 : [0,0], node2 : [0.117,0]}
                * {node1 : [0,0]}
        """
        x_opt=position
        TG={}
        L1=[]
        Struct=[]
        Rot={}
        for num,en in enumerate(self.connections_dfs):
            position1=(x_opt[en[0]][0],x_opt[en[0]][1])
            position2=(x_opt[en[1]][0],x_opt[en[1]][1])
            #tuple1 et 2 correspondent a la position des centres
            ne=self.connections.index(en)
            Rot[ne]={}
            if num==0:
                TG[en[0]]=self.meshes[en[0]].Contour(5)
            Struct.append(vm.Circle2D(vm.Point2D(position1),self.DF[ne][en[0]]/2.))
            TG[en[1]]=self.meshes[en[1]].Contour(5)
            Struct.append(vm.Circle2D(vm.Point2D(position2),self.DF[ne][en[1]]/2.))
            #Definition de la position angulaire initiale
            list_rot=self.InitialPosition(ne,en)
            if position2[0]==position1[0]:
                if position2[1]-position1[1]>0:
                    angle=math.pi/2.
                else:
                    angle=-math.pi/2.
            else:
                angle=-math.atan((position2[1]-position1[1])/(position2[0]-position1[0]))
            if num==0:
                Rot[ne][en[0]]=list_rot[0]-angle
                Rot[ne][en[1]]=list_rot[1]-angle
            else:
                for k1,v1 in Rot.items():
                    if en[0] in v1.keys():
                        Rot[ne][en[0]]=v1[en[0]]
                        delta_rot=Rot[ne][en[0]]-(list_rot[0]-angle)
                Rot[ne][en[1]]=list_rot[1]-angle-delta_rot*((self.meshes[en[0]].z)/(self.meshes[en[1]].z))
            sol=self.GearRotate([TG[en[0]],TG[en[1]]],[position1,position2],list_rot=[Rot[ne][en[0]],Rot[ne][en[1]]])
            if num==0:
                L1.extend(sol[0])
            L1.extend(sol[1])
        L1.extend(Struct)
#        G1=vm.Contour2D(L1)
#        G1.MPLPlot()
        return L1

    def JSON(self, filepath, indent = 2):
        with open(filepath+'.json', 'w') as j:
            json.dump(tools.StringifyDictKeys(self.Dict()), j, indent = indent)

    def Dict(self):
        d = {'name' : self.name} # TODO Change this to DessiaObject.__init__
        d['center_distance'] = self.center_distance
        d['connections'] = self.connections
        d['meshes'] = {}
        for num_mesh, mesh in self.meshes.items():
            d['meshes'][num_mesh] = mesh.Dict()
        d['keys_torque'] = []
        d['torque'] = []
        for keys, value in self.torque.items():
            d['keys_torque'].append(keys)
            d['torque'].append(value)
        d['cycle'] = self.cycle
        d['safety_factor'] = self.safety_factor
        return d

    @classmethod
    def DictToObject(cls, d):
        meshes = {}
        for num_mesh, mesh in d['meshes'].items():
            meshes[int(num_mesh)] = Mesh.DictToObject(mesh)
        torques = {}
        for keys, value in zip(d['keys_torque'], d['torque']):
            ki = ()
            for k in keys:
                if k.__class__ == str:
                    ki = ki + (int(k),)
                else:
                    ki = ki + (k,)
            torques[ki] = value
        cycle = {}
        for num_mesh, cy in d['cycle'].items():
            if num_mesh.__class__ == str:
                cycle[int(num_mesh)] = cy
            else:
                cycle[num_mesh] = cy
        connections = []
        for connection in d['connections']:
            connections.append(tuple(connection))
        mesh_combination = cls(center_distance = d['center_distance'],
                   connections = connections,
                   meshes = meshes,
                   torque = torques,
                   cycle = cycle,
                   safety_factor = d['safety_factor'],
                   name=d['name'])
        return mesh_combination


class MeshAssembly(DessiaObject):
    def __init__(self, connections, mesh_combinations, torque=None,
                 cycle=None, strong_links=None, safety_factor=1, name=''):

        self.connections = connections
        self.mesh_combinations = mesh_combinations
        self.torque = torque
        self.cycle = cycle
        self.strong_links = strong_links
        self.safety_factor = safety_factor

        self.center_distance = []
        for num_cd, list_connection in enumerate(self.connections):
            for num_mesh_iter, gs in enumerate(list_connection):
                valid = False
                for mesh_combination in mesh_combinations:
                    for num_mesh_local, gs_local in enumerate(mesh_combination.connections):
                        if set(gs) == set(gs_local):
                            self.center_distance.append(mesh_combination.center_distance[num_mesh_local])
                            valid = True
                        if valid:
                            break
                    if valid:
                        break
                if valid:
                    break

        transverse_pressure_angle = []
        for num_cd, list_connection in enumerate(self.connections):
            for num_mesh_iter, gs in enumerate(list_connection):
                valid = False
                for mesh_combination in mesh_combinations:
                    for num_mesh_local, gs_local in enumerate(mesh_combination.connections):
                        if set(gs) == set(gs_local):
                            transverse_pressure_angle.append(mesh_combination.transverse_pressure_angle[num_mesh_local])
                            valid = True
                        if valid:
                            break
                    if valid:
                        break

        list_gear = {}
        for mesh_combination in self.mesh_combinations:
            for num_mesh, mesh in mesh_combination.meshes.items():
                list_gear[num_mesh] = mesh
        coefficient_profile_shift = {}
        for num_mesh, mesh in list_gear.items():
            coefficient_profile_shift[num_mesh] = mesh.coefficient_profile_shift
        Z = {}
        for num_mesh, mesh in list_gear.items():
            Z[num_mesh] = mesh.z
        material = {}
        for num_mesh, mesh in list_gear.items():
            material[num_mesh] = mesh.material
        transverse_pressure_angle_rack = {}
        for num_mesh, mesh in list_gear.items():
            transverse_pressure_angle_rack[num_mesh] = mesh.rack.transverse_pressure_angle
        coeff_gear_addendum = {}
        for num_mesh, mesh in list_gear.items():
            coeff_gear_addendum[num_mesh] = mesh.rack.coeff_gear_addendum
        coeff_gear_dedendum = {}
        for num_mesh, mesh in list_gear.items():
            coeff_gear_dedendum[num_mesh] = mesh.rack.coeff_gear_dedendum
        coeff_root_radius = {}
        for num_mesh, mesh in list_gear.items():
            coeff_root_radius[num_mesh] = mesh.rack.coeff_root_radius
        coeff_circular_tooth_thickness = {}
        for num_mesh, mesh in list_gear.items():
            coeff_circular_tooth_thickness[num_mesh] = mesh.rack.coeff_circular_tooth_thickness

        self.general_data = []
        for num_graph,list_sub_graph in enumerate(self.sub_graph_dfs):
            num_mesh=0
            general_data={'Z': {}, 'connections': [],
                 'material':{},'torque':{},'cycle':{},
                 'safety_factor':safety_factor}
            input_data={'center_distance':[],'transverse_pressure_angle_0':0,
                 'coefficient_profile_shift':{},'transverse_pressure_angle_rack':{},
                 'coeff_gear_addendum':{},'coeff_gear_dedendum':{},
                 'coeff_root_radius':{},'coeff_circular_tooth_thickness':{}}
            li_connection=[]
            for num_cd, list_connection in enumerate(connections):
                for num_mesh_iter,gs in enumerate(list_connection):
                    if (gs in list_sub_graph) or (gs[::-1] in list_sub_graph):
                        li_connection.append(gs)
                        for num_gear in gs:
                            if num_gear in coefficient_profile_shift.keys():
                                input_data['coefficient_profile_shift'][num_gear] = coefficient_profile_shift[num_gear]
                            if num_gear in transverse_pressure_angle_rack.keys():
                                input_data['transverse_pressure_angle_rack'][num_gear] = transverse_pressure_angle_rack[num_gear]
                            if num_gear in coeff_gear_addendum.keys():
                                input_data['coeff_gear_addendum'][num_gear] = coeff_gear_addendum[num_gear]
                            if num_gear in coeff_gear_dedendum.keys():
                                input_data['coeff_gear_dedendum'][num_gear] = coeff_gear_dedendum[num_gear]
                            if num_gear in coeff_root_radius.keys():
                                input_data['coeff_root_radius'][num_gear] = coeff_root_radius[num_gear]
                            if num_gear in coeff_circular_tooth_thickness.keys():
                                input_data['coeff_circular_tooth_thickness'][num_gear]=coeff_circular_tooth_thickness[num_gear]
                            if num_gear in Z.keys():
                                general_data['Z'][num_gear]=Z[num_gear]
                            if num_gear in material.keys():
                                general_data['material'][num_gear]=material[num_gear]
                        if num_mesh == 0:
                            input_data['transverse_pressure_angle_0'] = transverse_pressure_angle[num_mesh]
                    num_mesh+=1
                input_data['center_distance'].append(self.center_distance[num_cd])
            general_data['connections']=li_connection
            for (eng1,eng2) in list_sub_graph:
                if (eng1,eng2) in torque.keys():
                    general_data['torque'][(eng1,eng2)]=torque[(eng1,eng2)]
                if (eng2,eng1) in torque.keys():
                    general_data['torque'][(eng2,eng1)]=torque[(eng2,eng1)]
                if eng1 not in general_data['cycle'].keys():
                    general_data['cycle'][eng1]=cycle[eng1]
                if eng2 not in general_data['cycle'].keys():
                    general_data['cycle'][eng2]=cycle[eng2]
            self.general_data.append(general_data)

        DessiaObject.__init__(self, name=name)

    def __eq__(self, other_eb):
        equal = (self.connections == other_eb.connections
                 and self.mesh_combinations == other_eb.mesh_combinations
                 and self.torque == other_eb.torque
                 and self.cycle == other_eb.cycle
                 and self.strong_links == other_eb.strong_links
                 and self.safety_factor == other_eb.safety_factor)
        return equal

    def __hash__(self):
        ma_hash = 0
        for mesh_combination in self.mesh_combinations:
            ma_hash += hash(mesh_combination)
        if self.strong_links is not None:
            for strong_link in self.strong_links:
                for sl in strong_link:
                    ma_hash += hash(strong_link)
        return ma_hash

    def check(self):
        valid = True
        for mesh_combination in self.mesh_combinations:
            valid = (valid and mesh_combination.check())
        return valid

    @classmethod
    def create(cls, center_distance, connections, transverse_pressure_angle,
                 coefficient_profile_shift, transverse_pressure_angle_rack,
                 coeff_gear_addendum, coeff_gear_dedendum, coeff_root_radius,
                 coeff_circular_tooth_thickness, Z, strong_links=None, material=None,
                 torque=None, cycle=None,
                 safety_factor=1):


#        self.connections = connections
#        self.center_distance = center_distance
#        self.mesh_combinations = []
#        self.general_data = []

        mesh_combinations = []
        output_data = []

        graph_dfs,_ = gear_graph_simple(connections)
        for num_graph,list_sub_graph in enumerate(graph_dfs):
            num_mesh=0
            general_data={'Z': {}, 'connections': [],
                 'material':{},'torque':{},'cycle':{},
                 'safety_factor':safety_factor}
            input_data={'center_distance':[],'transverse_pressure_angle_0':0,
                 'coefficient_profile_shift':{},'transverse_pressure_angle_rack':{},
                 'coeff_gear_addendum':{},'coeff_gear_dedendum':{},
                 'coeff_root_radius':{},'coeff_circular_tooth_thickness':{}}
            li_connection=[]
            for num_cd, list_connection in enumerate(connections):
                for num_mesh_iter,gs in enumerate(list_connection):
                    if (gs in list_sub_graph) or (gs[::-1] in list_sub_graph):
                        li_connection.append(gs)
                        for num_gear in gs:
                            if num_gear in coefficient_profile_shift.keys():
                                input_data['coefficient_profile_shift'][num_gear] = coefficient_profile_shift[num_gear]
                            if num_gear in transverse_pressure_angle_rack.keys():
                                input_data['transverse_pressure_angle_rack'][num_gear] = transverse_pressure_angle_rack[num_gear]
                            if num_gear in coeff_gear_addendum.keys():
                                input_data['coeff_gear_addendum'][num_gear] = coeff_gear_addendum[num_gear]
                            if num_gear in coeff_gear_dedendum.keys():
                                input_data['coeff_gear_dedendum'][num_gear] = coeff_gear_dedendum[num_gear]
                            if num_gear in coeff_root_radius.keys():
                                input_data['coeff_root_radius'][num_gear] = coeff_root_radius[num_gear]
                            if num_gear in coeff_circular_tooth_thickness.keys():
                                input_data['coeff_circular_tooth_thickness'][num_gear]=coeff_circular_tooth_thickness[num_gear]
                            if num_gear in Z.keys():
                                general_data['Z'][num_gear]=Z[num_gear]
                            if num_gear in material.keys():
                                general_data['material'][num_gear]=material[num_gear]
                        if num_mesh == 0:
                            input_data['transverse_pressure_angle_0'] = transverse_pressure_angle[num_mesh]
                    num_mesh+=1
                input_data['center_distance'].append(center_distance[num_cd])
            general_data['connections']=li_connection
            for (eng1,eng2) in list_sub_graph:
                if (eng1,eng2) in torque.keys():
                    general_data['torque'][(eng1,eng2)]=torque[(eng1,eng2)]
                if (eng2,eng1) in torque.keys():
                    general_data['torque'][(eng2,eng1)]=torque[(eng2,eng1)]
                if eng1 not in general_data['cycle'].keys():
                    general_data['cycle'][eng1]=cycle[eng1]
                if eng2 not in general_data['cycle'].keys():
                    general_data['cycle'][eng2]=cycle[eng2]

            output_data.append(general_data)
            xt = dict(list(input_data.items()) + list(general_data.items()))
            mesh_combinations.append(MeshCombination.create(**xt))
        mesh_assembly = cls(connections, mesh_combinations, torque, cycle,
                            strong_links, safety_factor)
        return mesh_assembly

    def _get_graph_dfs(self):
        _graph_dfs,_ = gear_graph_simple(self.connections)
        return _graph_dfs
    sub_graph_dfs = property(_get_graph_dfs)

    def _get_list_gear(self):
        _,_list_gear = gear_graph_simple(self.connections)
        return _list_gear
    list_gear = property(_get_list_gear)

    def SVGExport(self,name,position):
        centers=self.PosAxis(position)
        L=[]
        for mesh_assembly_iter in self.mesh_combinations:
            position_svg={}
            for num_gear,pos in centers.items():
                if num_gear in mesh_assembly_iter.Z.keys():
                    position_svg[num_gear]=pos
            L.extend(mesh_assembly_iter.SVGExport('gear',position_svg))
        G1=vm.Contour2D(L)
        G1.MPLPlot()

    def FreeCADExport(self, fcstd_filepath, centers = {}, axis = (1,0,0), export_types=['fcstd'], python_path = 'python',
                      path_lib_freecad = '/usr/lib/freecad/lib'):
        """ Export 3D volume to FreeCAD

        :param file_path: file path for the freecad file
        :param center: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
        :param axis: direction of gear mesh rotation

        :results: export of a FreeCAD file
        """
        for ma in self.mesh_combinations:
            ma.FreeCADExport(fcstd_filepath, centers, axis, python_path, path_lib_freecad, export_types)

    def Update(self, optimizer_data):
        output_x=[]
        for num_graph,list_sub_graph in enumerate(self.sub_graph_dfs):
            num_mesh = 0
            input_data={'center_distance':[],'transverse_pressure_angle_0':[],
                 'coefficient_profile_shift':{},'transverse_pressure_angle_rack':{},
                 'coeff_gear_addendum':{},'coeff_gear_dedendum':{},
                 'coeff_root_radius':{},'coeff_circular_tooth_thickness':{}}
            li_connection=[]
            for num_cd,list_connection in enumerate(self.connections):
                for num_mesh_iter,(eng1,eng2) in enumerate(list_connection):
                    if ((eng1,eng2) in list_sub_graph) or ((eng2,eng1) in list_sub_graph):
                        li_connection.append((eng1,eng2))
                        for key,list_value in optimizer_data.items():
                            if key in ['coefficient_profile_shift',
                                       'transverse_pressure_angle_rack',
                                       'coeff_gear_addendum','coeff_gear_dedendum',
                                       'coeff_root_radius','coeff_circular_tooth_thickness']:
                                input_data[key][eng1]=optimizer_data[key][eng1]
                                input_data[key][eng2]=optimizer_data[key][eng2]
                            elif key in ['center_distance']:
                                input_data[key].append(optimizer_data[key][num_cd])
                            elif key in ['transverse_pressure_angle']:
                                input_data['transverse_pressure_angle_0'].append(optimizer_data[key][num_mesh])
                    num_mesh += 1
            input_data['transverse_pressure_angle_0'] = input_data['transverse_pressure_angle_0'][0]
            xt = dict(list(input_data.items())+list(self.general_data[num_graph].items()))
            output_x.append(xt)

#            if self.save!=optimizer_data:
            self.mesh_combinations[num_graph].Update(**xt)
        return output_x

    def PosAxis(self,position):
        # Definition of the initial center for all gear (when not given by the user)

        gear_graph=nx.Graph()
        gear_graph.add_nodes_from(self.list_gear)
        for num_cd,list_connections in enumerate(self.connections):
            (eng1_m,eng2_m)=list_connections[0]
            if len(list_connections)>1:
                for (eng1,eng2) in list_connections[1:]:
                    gear_graph.add_edges_from([(eng1_m,eng1),(eng2_m,eng2)])
                    eng1_m=eng1
                    eng2_m=eng2
#        list_line=list(nx.connected_component_subgraphs(gear_graph))
        list_line = [gear_graph.subgraph(c).copy() for c in nx.connected_components(gear_graph)]
        dict_line={}
        for num_line,list_num_eng in enumerate(list_line):
            for num_eng in list_num_eng:
                dict_line[num_eng]=num_line
        def fun(x):
            obj=0
            for num_cd,list_connections in enumerate(self.connections):
                eng1=dict_line[list_connections[0][0]]
                eng2=dict_line[list_connections[0][1]]
                obj+=(((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5-self.center_distance[num_cd])**2
            return obj
        def eg(x):
            ine=[]
            for k,val in position.items():
                key=dict_line[k]
                ine.append(x[2*int(key)]-val[0])
                ine.append(x[2*int(key)+1]-val[1])
            return ine
        def ineg(x):
            ine=[]
            for num_cd,list_connections in enumerate(self.connections):
                eng1=dict_line[list_connections[0][0]]
                eng2=dict_line[list_connections[0][1]]
                ine.append(((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5-0.999*self.center_distance[num_cd])
                ine.append(1.001*self.center_distance[num_cd]-((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5)
            return ine
        cons = ({'type': 'eq','fun' : eg},{'type': 'ineq','fun' : ineg})
        drap=1
        while drap==1:
            x0=tuple(npy.random.random(2*len(list_line))*1)
            Bound=[[0,1]]*(len(list_line)*2)
            res = minimize(fun,x0, method='SLSQP', bounds=Bound,constraints=cons)
            if (min(ineg(res.x))>0) and (max(eg(res.x))<1e-7):
                drap=0
        x_opt=res.x
        centers={}
        for num_pos,num_eng in enumerate(self.list_gear):
            opt_pos=dict_line[num_eng]
            centers[num_eng]=[x_opt[2*opt_pos],x_opt[2*opt_pos+1]]
        return centers

    def Dict(self):
        d = {'name' : self.name} # TODO Change this to DessiaObject.__init__
        d['connections'] = self.connections
        d['mesh_combinations'] = []
        for mesh_combination in self.mesh_combinations:
            d['mesh_combinations'].append(mesh_combination.Dict())
        d['keys_torque'] = []
        d['torque'] = []
        for keys, value in self.torque.items():
            d['keys_torque'].append(keys)
            d['torque'].append(value)
#        d['torque'] = self.torque
        d['cycle'] = self.cycle
        d['strong_links'] = self.strong_links
        d['safety_factor'] = self.safety_factor
        return d

    @classmethod
    def DictToObject(cls, d):
        mesh_combinations = []
        for mesh_combination in d['mesh_combinations']:
            mesh_combinations.append(MeshCombination.DictToObject(mesh_combination))
        torques = {}
        for keys, value in zip(d['keys_torque'], d['torque']):
            ki = ()
            for k in keys:
                if k.__class__ == str:
                    ki = ki + (int(k),)
                else:
                    ki = ki + (k,)
            torques[ki] = value
        cycle = {}
        for num_mesh, cy in d['cycle'].items():
            if num_mesh.__class__ == str:
                cycle[int(num_mesh)] = cy
            else:
                cycle[num_mesh] = cy
        mesh_assembly = cls(connections = d['connections'],
                   mesh_combinations = mesh_combinations,
                   torque = torques,
                   cycle = cycle,
                   strong_links = d['strong_links'],
                   safety_factor = d['safety_factor'],
                   name=d['name'])
        return mesh_assembly

def gear_graph_simple(connections):
    # NetworkX graph construction
    list_gear=[] # list of all gears
    compt_mesh=0 # number of gear mesh
    for gs in connections:
        for (eng1,eng2) in gs:
            compt_mesh+=1
            if eng1 not in list_gear:
                list_gear.append(eng1)
            if eng2 not in list_gear:
                list_gear.append(eng2)
    # Construction of one graph include all different connection type (gear_mesh, same_speed, same_shaft)
    gear_graph=nx.Graph()
    gear_graph.add_nodes_from(list_gear)
    for list_edge in connections:
        gear_graph.add_edges_from(list_edge)
#    sub_graph=list(nx.connected_component_subgraphs(gear_graph))
    sub_graph = [gear_graph.subgraph(c).copy() for c in nx.connected_components(gear_graph)]
    sub_graph_dfs = []
    for s_graph in sub_graph:
        node_init=list(s_graph.nodes())[0]
        sub_graph_dfs.append(list(nx.dfs_edges(s_graph,node_init)))
    return sub_graph_dfs,list_gear

def gear_graph_complex(connections,strong_link):
    # Construction of one graph include all different connection type (gear_mesh, same_speed, same_shaft)
    _,list_gear=gear_graph_simple(connections)
    gear_graph=nx.Graph()
    gear_graph.add_nodes_from(list_gear)
    for list_edge in connections:
        gear_graph.add_edges_from(list_edge,typ='gear_mesh')
        li_shaft1=[]
        li_shaft2=[]
        for eng1,eng2 in list_edge:
            li_shaft1.append(eng1)
            li_shaft2.append(eng2)
        if len(li_shaft1)>1:
            for pos_gear,num_gear in enumerate(li_shaft1[1:]):
                valid_strong_ling=False
                for list_strong_link in strong_link:
                    if (num_gear in list_strong_link) and (li_shaft1[pos_gear] in list_strong_link):
                        valid_strong_ling=True
                if valid_strong_ling:
                    gear_graph.add_edges_from([(num_gear,li_shaft1[pos_gear])],typ='same_speed')
                else:
                    gear_graph.add_edges_from([(num_gear,li_shaft1[pos_gear])],typ='same_shaft')
        if len(li_shaft2)>1:
            for pos_gear,num_gear in enumerate(li_shaft2[1:]):
                valid_strong_ling=False
                for list_strong_link in strong_link:
                    if (num_gear in list_strong_link) and (li_shaft2[pos_gear] in list_strong_link):
                        valid_strong_ling=True
                if valid_strong_ling:
                    gear_graph.add_edges_from([(num_gear,li_shaft2[pos_gear])],typ='same_speed')
                else:
                    gear_graph.add_edges_from([(num_gear,li_shaft2[pos_gear])],typ='same_shaft')
    connections_dfs=list(nx.dfs_edges(gear_graph,list_gear[0]))
    # construction of a graph without same_shaft attribute
    gear_graph_kinematic=copy.deepcopy(gear_graph)
    for edge,typ in nx.get_edge_attributes(gear_graph_kinematic,'typ').items():
        if typ=='same_shaft':
            gear_graph_kinematic.remove_edges_from([edge])
    connections_kinematic_dfs=list(nx.dfs_edges(gear_graph_kinematic,list_gear[0]))
    return connections_dfs,connections_kinematic_dfs,gear_graph

class ValidGearDiameterError(Exception):
    def __init__(self):
        super().__init__('Fail base diameter is greater than pitch diameter')