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
import volmdlr.primitives3d as primitives3D
import volmdlr.primitives2d as primitives2D
import plot_data as vmp
import math
from scipy.linalg import norm
from scipy.optimize import fsolve, minimize
import networkx as nx
import matplotlib.pyplot as plt
from dessia_common import DessiaObject
import mechanical_components.tools as tools
import json
import copy
from typing import  List, Tuple
from scipy.optimize import fsolve


class Data(DessiaObject):
    _standalone_in_db = False

    def __init__(self, data: List[Tuple[float, float]], x: str, y: str):
        self.data = data
        self.x = x
        self.y = y
        DessiaObject.__init__(self, name='data')




#data_coeff_YB_Iso
evol_coeff_yb_iso = Data(data=[[0.0, 1.0029325508401201],
                               [4.701492563229561, 0.9310850480431024],
                               [9.104477651442416, 0.8782991233021732],
                               [14.5522388104227, 0.8240469255458759],
                               [19.328358165913905, 0.784457481990179],
                               [23.955224059269884, 0.7609970656504502],
                               [28.507462609617956, 0.7521994155347822],
                               [34.029850499665855, 0.7507331425194141],
                               [40.0, 0.7492668574805859]
                               ], x='Linear', y='Linear')

#data_wohler_curve
wohler_hardened_alloy_steel = Data(data=[[4.296196199237153, 2.9797762011105589],
                                         [4.824840106199563, 2.9413306094362142],
                                         [5.3344338175705674, 2.908892154601565],
                                         [6.115493253679078, 2.8632380197445122],
                                         [6.596511629990596, 2.8560294765618042],
                                         [7.144205815889171, 2.8536266428508523],
                                         [7.691899918442984, 2.8524252154829133],
                                         [8.010991340520903, 2.8524252154829133]
                                         ], x='Log', y='Log')

wohler_nitrided_alloy_steel = Data(data=[[4.104865629699472, 2.9252942042661974],
                                         [4.568697315952783, 2.8521640228225367],
                                         [4.887581626297173, 2.8046294185503593],
                                         [5.438381821440599, 2.7900033864666123],
                                         [6.402282079596832, 2.7918316299646175],
                                         [7.264719174616821, 2.7918316299646175],
                                         [7.989456220850952, 2.793659894487549]
                                         ], x='Log', y='Log')

wohler_through_hardened_steel = Data(data=[[4.172369719531124, 2.895676495604088],
                                           [4.677200861168087, 2.7983611100752137],
                                           [4.9677168648417585, 2.741894170956562],
                                           [5.329671247836526, 2.6842258044699714],
                                           [5.439210101685194, 2.672211551815507],
                                           [6.091680488353632, 2.6734129791834462],
                                           [7.139443246155129, 2.671010124447568],
                                           [8.00146620105282, 2.6758158339193243]
                                           ], x='Log', y='Log')

wohler_surface_hardened_steel = Data(data=[[4.281908490035029, 2.7611169667937343],
                                           [4.701013626493532, 2.6998443182033265],
                                           [5.015342395492649, 2.6553916107142128],
                                           [5.358246582896013, 2.6109389032250994],
                                           [5.620187251510196, 2.5857089915731581],
                                           [6.020242109032534, 2.5748961873115592],
                                           [6.567936294931109, 2.5748961873115592],
                                           [7.263269725861159, 2.5772990210225108],
                                           [7.996703631318779, 2.5772990210225108]
                                           ], x='Log', y='Log')

wohler_carbon_steel = Data(data=[[4.307791955971963, 2.6419147590563592],
                                 [5.242702822291173, 2.535876005424268],
                                 [5.938450393343521, 2.4700588400224806],
                                 [6.518240063668731, 2.431665495290182],
                                 [7.221234961844144, 2.4334937598131132],
                                 [7.989456220850952, 2.4353220033111185]
                                 ], x='Log', y='Log')

wohler_cast_iron = Data(data=[[4.307791955971963, 2.6419147590563592],
                              [5.242702822291173, 2.535876005424268],
                              [5.938450393343521, 2.4700588400224806],
                              [6.518240063668731, 2.431665495290182],
                              [7.221234961844144, 2.4334937598131132],
                              [7.989456220850952, 2.4353220033111185]
                              ], x='Log', y='Log')

wohler_bronze = Data(data=[[4.307791955971963, 2.6419147590563592],
                           [5.242702822291173, 2.535876005424268],
                           [5.938450393343521, 2.4700588400224806],
                           [6.518240063668731, 2.431665495290182],
                           [7.221234961844144, 2.4334937598131132],
                           [7.989456220850952, 2.4353220033111185]
                           ], x='Log', y='Log')

wohler_grey_iron = Data(data=[[4.307791955971963, 2.6419147590563592],
                              [5.242702822291173, 2.535876005424268],
                              [5.938450393343521, 2.4700588400224806],
                              [6.518240063668731, 2.431665495290182],
                              [7.221234961844144, 2.4334937598131132],
                              [7.989456220850952, 2.4353220033111185]
                              ], x='Log', y='Log')
#data_gear_material

sigma_hardened_alloy_steel = Data(data=[[1.8422104370714443, 1.4645831828946267],
                                        [1.948612010770208, 1.5219116983411152],
                                        [2.0605171321606295, 1.5810895335609718],
                                        [2.141235568740199, 1.6254729099758645]
                                        ], x='Log', y='Log')

sigma_nitrided_alloy_steel = Data(data=[[1.8458794622934307, 1.4349942652846983],
                                        [1.943108482795906, 1.488624180937243],
                                        [2.0201578941534892, 1.5274596179084272],
                                        [2.128393990321924, 1.5866374531282839]
                                        ], x='Log', y='Log')

sigma_through_hardened_steel = Data(data=[[1.7798371068844516, 1.292597616678765],
                                          [1.921094370898698, 1.3850629693024938],
                                          [2.032999472571764, 1.4627338829976548],
                                          [2.1650841833897223, 1.5533499158480155]
                                          ], x='Log', y='Log')

sigma_surface_hardened_steel = Data(data=[[1.8312033811228403, 1.115064130895591],
                                          [1.932101426847302, 1.200132264055036],
                                          [2.038503000546066, 1.2852003773380847]
                                          ], x='Log', y='Log')

sigma_carbon_steel = Data(data=[[1.677104538690319, 1.1002696720906269],
                                [1.7633265032441903, 1.1723926463420797],
                                [1.8385414118494579, 1.2389677010262203],
                                [1.8844041581135444, 1.2796524577707729]
                                ], x='Log', y='Log')

sigma_cast_iron = Data(data=[[1.4734739247717241, 0.922736186307453],
                             [1.5468543306246763, 0.9837633214242817],
                             [1.6073931580593532, 1.0336946174064863],
                             [1.6404143456225206, 1.0688314545837265]
                             ], x='Log', y='Log')

sigma_bronze = Data(data=[[1.313871566195314, 0.7858874572688317],
                          [1.3890864826875238, 0.8487638922826322],
                          [1.4294457009773085, 0.8802021097895326],
                          [1.4551288380965028, 0.9097910273994609]
                          ], x='Log', y='Log')

sigma_grey_iron = Data(data=[[1.354230792372041, 0.7100658633470387],
                             [1.4276111785076375, 0.7766409180311793],
                             [1.4936535339166166, 0.84691459238566],
                             [1.5431853054026896, 0.8986951882648367],
                             [1.5725374677438706, 0.933832025442077]
                             ], x='Log', y='Log')

class Material(DessiaObject):
    """
    Gear material

    :param volumic_mass: A float to define the gear volumic mass
    :param data_coeff_YB_Iso: a dictionary to define the YB parameter of the ISO description
    :param data_wohler_curve: a dictionary to define the wohler slope of the ISO description
    :param data_gear_material: a dictionary to define the maximum gear stress

    :data_coeff_YB_Iso: - **'data'** matrix define points of the YB curve in the plane (YB, helix_angle)
        - **'x'** string define the x axis evolution ('Log' or 'Linear')
        - **'y'** string define the y axis evolution ('Log' or 'Linear')

    :data_wohler_curve: - **'data'** matrix define points of the wohler slope in the plane (wohler slope, number of cycle)
        - **'x'** string define the x axis evolution ('Log' or 'Linear')
        - **'y'** string define the y axis evolution ('Log' or 'Linear')

    :data_gear_material: - **'data'** matrix define points of the maximum gear stress (maximum gear stress, wohler slope)
        - **'x'** string define the x axis evolution ('Log' or 'Linear')
        - **'y'** string define the y axis evolution ('Log' or 'Linear')

    >>> volumic_mass=7800
    >>> data_coeff_YB_Iso={'data':[[0.0,1.0029325508401201],
                           [4.701492563229561,0.9310850480431024],
                           [23.955224059269884,0.7609970656504502],
                           [40.0,0.7492668574805859]
                          ], 'x':'Linear','y':'Linear'}
    >>> data_wohler_curve={'data':[[4.307791955971963,1.6419147590563592],
                       [6.518240063668731,1.431665495290182],
                       [7.989456220850952,1.4353220033111185]
                      ], 'x':'Log','y':'Log'}
    >>> data_gear_material={'data':[[1.313871566195314,0.7858874572688317],
                      [1.4294457009773085,0.8802021097895326],
                      [1.4551288380965028,0.9097910273994609]
                     ], 'x':'Log','y':'Log'}
    >>> material1=Material(volumic_mass, data_coeff_YB_Iso,
                           data_wohler_curve, data_gear_material)
    """
    _standalone_in_db = False


    def __init__(self, volumic_mass: float = 7850, data_coeff_YB_Iso: Data = evol_coeff_yb_iso, data_wohler_curve: Data = wohler_hardened_alloy_steel,
                 data_gear_material: Data = sigma_hardened_alloy_steel, name: str = ''):
        self.volumic_mass = volumic_mass
        self.data_coeff_YB_Iso = data_coeff_YB_Iso
        self.data_wohler_curve = data_wohler_curve
        self.data_gear_material = data_gear_material

        DessiaObject.__init__(self, name=name)

    # def __eq__(self, other_eb):
    #     equal = (self.volumic_mass == other_eb.volumic_mass
    #              and self.data_coeff_YB_Iso == other_eb.data_coeff_YB_Iso
    #              and self.data_wohler_curve == other_eb.data_wohler_curve
    #              and self.data_gear_material == other_eb.data_gear_material)
    #     return equal

    def __hash__(self):
        material_hash = hash(self.volumic_mass)
        return material_hash

    def FunCoeff(self, x, data, type_x='Linear', type_y='Linear'):
        """ Interpolation of material data

        :param x: value of the interpolation
        :param data: dictionary of the input data
        :param type_x: type of the x axis of the data matrix ('Log' or 'Linear')
        :param type_y: type of the y axis of the data matrix ('Log' or 'Linear')

        :returns:  interpolation value

        >>> interp1=material1.FunCoeff(x = 5.2,data = data_wohler_curve,
                                       type_x = 'Log',type_y = 'Log')
        """
        if type_x == 'Log':

            x = math.log10(abs(x)) #TODO

        f = interpolate.interp1d(list(data[:, 0]), list(data[:, 1]),
                                 fill_value='extrapolate')
        sol = float(f(x))
        if type_y == 'Log':
            sol = 10**sol
        return sol

    # def Dict(self):

    #     d = {'name' : self.name} # TODO Change this to DessiaObject.__init__
    #     d['volumic_mass'] = self.volumic_mass
    #     d['data_coeff_YB_Iso'] = self.data_coeff_YB_Iso
    #     d['data_wohler_curve'] = self.data_wohler_curve
    #     d['data_gear_material'] = self.data_gear_material
    #     return d

    # @classmethod
    # def dict_to_object(cls, d):
    #     material = cls(volumic_mass=d['volumic_mass'],
    #                    data_coeff_YB_Iso=d['data_coeff_YB_Iso'],
    #                    data_wohler_curve=d['data_wohler_curve'],
    #                    data_gear_material=d['data_gear_material'],
    #                    name=d['name'])
    #     return material

hardened_alloy_steel = Material(7850, evol_coeff_yb_iso,
                                wohler_hardened_alloy_steel,
                                sigma_hardened_alloy_steel,
                                name='Hardened alloy steel')

nitrided_alloy_steel = Material(7850, evol_coeff_yb_iso, wohler_nitrided_alloy_steel,
                                sigma_nitrided_alloy_steel,
                                name='Nitrided alloy steel')

through_hardened_steel = Material(7850, evol_coeff_yb_iso,
                                  wohler_through_hardened_steel,
                                  sigma_through_hardened_steel,
                                  name='Through hardened steel')

surface_hardened_steel = Material(7850, evol_coeff_yb_iso,
                                  wohler_surface_hardened_steel,
                                  sigma_surface_hardened_steel,
                                  name='Surface hardened steel')

carbon_steel = Material(7850, evol_coeff_yb_iso, wohler_carbon_steel, sigma_carbon_steel,
                        name='Carbon steel')

cast_iron = Material(7200, evol_coeff_yb_iso, wohler_cast_iron, sigma_cast_iron,
                     name='Cast iron')

bronze = Material(8200, evol_coeff_yb_iso, wohler_bronze, sigma_bronze,
                  name='Bronze')

grey_iron = Material(7200, evol_coeff_yb_iso, wohler_grey_iron, sigma_grey_iron,
                     name='Grey iron')






class Rack(DessiaObject):
    """
    Gear rack definition

    :param transverse_pressure_angle_0: definition of the transverse pressure angle of the rack
    :type transverse_pressure_angle: radian
    :param coeff_gear_addendum: update of the gear addendum coefficient (gear_addendum = coeff_gear_addendum*module)
    :param coeff_gear_dedendum: update of the gear dedendum coefficient (gear_dedendum = coeff_gear_dedendum*module)
    :param coeff_root_radius: update of the root radius coefficient (root_radius = coeff_root_radius*module)
    :param coeff_circular_tooth_thickness: update of the circular tooth thickness coefficient (circular_tooth_thickness = coeff_circular_tooth_thickness*transverse_radial_pitch)
    :param helix_angle: float  define the helix_angle of the rack
    >>> Rack1=Rack(20/180.*math.pi) #definition of an ISO rack
    """
    _standalone_in_db = True
    _eq_is_data_eq = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']

    def __init__(self, transverse_pressure_angle_0: float, module: float = None,
                 coeff_gear_addendum: float = 1, coeff_gear_dedendum: float = 1.25,
                 coeff_root_radius: float = 0.38, coeff_circular_tooth_thickness: float = 0.5,
                 helix_angle: float = 0.0, name: str = ''):
        self.transverse_pressure_angle_0 = transverse_pressure_angle_0
        self.module = module
        self.coeff_gear_addendum = coeff_gear_addendum
        self.coeff_gear_dedendum = coeff_gear_dedendum
        self.coeff_root_radius = coeff_root_radius
        self.coeff_circular_tooth_thickness = coeff_circular_tooth_thickness
        self.helix_angle = helix_angle

        if module is not None:
            self.update(module, transverse_pressure_angle_0, coeff_gear_addendum,
                        coeff_gear_dedendum, coeff_root_radius,
                        coeff_circular_tooth_thickness)

        DessiaObject.__init__(self, name=name)

    def rack_param(self, transverse_pressure_angle_0, coeff_gear_addendum,
                   coeff_gear_dedendum, coeff_root_radius, coeff_circular_tooth_thickness):

        self.transverse_pressure_angle_0 = transverse_pressure_angle_0
        self.transverse_radial_pitch = self.module*math.pi
        self.gear_addendum = coeff_gear_addendum*self.module
        self.gear_dedendum = coeff_gear_dedendum*self.module


        self.root_radius = coeff_root_radius*self.module
        self.circular_tooth_thickness = coeff_circular_tooth_thickness*self.transverse_radial_pitch

        self.tooth_space = self.transverse_radial_pitch-self.circular_tooth_thickness
        self.whole_depth = self.gear_addendum+self.gear_dedendum
        self.clearance = self.root_radius-self.root_radius*math.sin(self.transverse_pressure_angle_0)


        # trochoide parameter
        self.a = (self.tooth_space/2.
                  - self.gear_dedendum * math.tan(self.transverse_pressure_angle_0)
                  - self.root_radius * math.tan(0.5*math.atan(math.cos(self.transverse_pressure_angle_0)
                                                              /(math.sin(self.transverse_pressure_angle_0)))))
        self.b = self.gear_dedendum - self.root_radius

    def update(self, module, transverse_pressure_angle_0=None, coeff_gear_addendum=None,
               coeff_gear_dedendum=None, coeff_root_radius=None,
               coeff_circular_tooth_thickness=None):
        """
        update of the gear rack

        :param module: update of the module of the rack define on the pitch factory diameter
        :type module: m
        :param transverse_pressure_angle: update of the transverse pressure angle of the rack
        :type transverse_pressure_angle: radian
        :param coeff_gear_addendum: update of the gear addendum coefficient (gear_addendum = coeff_gear_addendum*module)
        :param coeff_gear_dedendum: update of the gear dedendum coefficient (gear_dedendum = coeff_gear_dedendum*module) (top of the rack)
        :param coeff_root_radius: update of the root radius coefficient (root_radius = coeff_root_radius*module)
        :param coeff_circular_tooth_thickness: update of the circular tooth thickness coefficient (circular_tooth_thickness = coeff_circular_tooth_thickness*transverse_radial_pitch)

        >>> input={'module':2*1e-3,'transverse_pressure_angle':21/180.*math.pi}
        >>> Rack1.update(**input) # update of the rack definition
        """
        if transverse_pressure_angle_0 == None:
            transverse_pressure_angle_0 = self.transverse_pressure_angle_0
        if coeff_gear_addendum == None:
            coeff_gear_addendum = self.coeff_gear_addendum
        if coeff_gear_dedendum == None:
            coeff_gear_dedendum = self.coeff_gear_dedendum
        if coeff_root_radius == None:
            coeff_root_radius = self.coeff_root_radius
        if coeff_circular_tooth_thickness == None:
            coeff_circular_tooth_thickness = self.coeff_circular_tooth_thickness
        self.module = module
        self.coeff_gear_addendum = coeff_gear_addendum
        self.coeff_gear_dedendum = coeff_gear_dedendum
        self.coeff_root_radius = coeff_root_radius
        self.coeff_circular_tooth_thickness = coeff_circular_tooth_thickness

        self.rack_param(transverse_pressure_angle_0, coeff_gear_addendum,
                        coeff_gear_dedendum, coeff_root_radius, coeff_circular_tooth_thickness)


    def update_helix_angle(self, helix_angle):
        """
        update of the gear rack

        :param module: update of the module of the rack define on the pitch factory diameter
        :type module: m
        :param transverse_pressure_angle: update of the transverse pressure angle of the rack
        :type transverse_pressure_angle: radian
        :param coeff_gear_addendum: update of the gear addendum coefficient (gear_addendum = coeff_gear_addendum*module)
        :param coeff_gear_dedendum: update of the gear dedendum coefficient (gear_dedendum = coeff_gear_dedendum*module) (top of the rack)
        :param coeff_root_radius: update of the root radius coefficient (root_radius = coeff_root_radius*module)
        :param coeff_circular_tooth_thickness: update of the circular tooth thickness coefficient (circular_tooth_thickness = coeff_circular_tooth_thickness*transverse_radial_pitch)

        >>> input={'module':2*1e-3,'transverse_pressure_angle':21/180.*math.pi}
        >>> Rack1.update(**input) # update of the rack definition
        """
        self.helix_angle = helix_angle

    ### Optimization Method

    def check_rack_viable(self):
        """ Check the viability of the rack toward the top and the root

        :results: boolean variable, and a list of element to be positive for the optimizer
        """
        list_ineq = []

        list_ineq.append(abs(self.transverse_radial_pitch)-abs(self.circular_tooth_thickness)
                         -2*abs(self.gear_dedendum)*math.tan(self.transverse_pressure_angle_0)
                         -2*(abs(self.root_radius)*math.cos(self.transverse_pressure_angle_0)-math.tan(self.transverse_pressure_angle_0)
                             *abs(self.root_radius)*(1-math.sin(self.transverse_pressure_angle_0))))
        list_ineq.append(abs(self.circular_tooth_thickness)-2*(abs(self.gear_addendum)*math.tan(self.transverse_pressure_angle_0)))
        check = False
        if min(list_ineq) > 0:
            check = True
        return check, list_ineq

    def liste_ineq(self):
        """ Compilation method for inequality list used by the optimizer

        :results: vector of data that should be positive
        """
        check, ineq = self.CheckRackViable
        return ineq

    def contour(self, number_pattern):
        """ Construction of the volmdr 2D rack profile

        :param number_pattern: number of rack pattern to define
        """
        p1 = vm.Point2D((0, 0))
        p2 = p1.Translation((self.gear_addendum*math.tan(self.transverse_pressure_angle_0), self.gear_addendum))
        p4 = p1.Translation((self.circular_tooth_thickness, 0))
        p3 = p4.Translation((-self.gear_addendum*math.tan(self.transverse_pressure_angle_0), self.gear_addendum))
        p5 = p4.Translation((self.gear_dedendum*math.tan(self.transverse_pressure_angle_0), -self.gear_dedendum))
        p7 = p4.Translation((self.tooth_space, 0))
        p6 = p7.Translation((-self.gear_dedendum*math.tan(self.transverse_pressure_angle_0), -self.gear_dedendum))
        L = primitives2D.OpenedRoundedLineSegments2D([p1, p2, p3, p4, p5, p6, p7], {4:self.root_radius, 5:self.root_radius}, False)

        Rack_Elem = []
        for i in range(number_pattern):
            Rack_Elem.append(L.Translation(((i)*(p7.vector-p1.vector))))
        p10 = Rack_Elem[0].points[0]
        p15 = Rack_Elem[-1].points[-1]
        p11 = p10.Translation((-self.circular_tooth_thickness, 0))
        p12 = p11.Translation((0, 2*self.whole_depth))
        p14 = p15.Translation((self.circular_tooth_thickness, 0))
        p13 = p14.Translation((0, 2*self.whole_depth))
        Rack_Elem.append(primitives2D.OpenedRoundedLineSegments2D([p10, p11, p12, p13, p14, p15], {}, False))

        return Rack_Elem

    def plot(self, number_pattern):
        """ Plot function of the rack

        :param number_pattern: number of rack pattern to draw
        """
        Rack_Elem = self.Contour(number_pattern)
        RackElem = vm.Contour2D(Rack_Elem)
        RackElem.MPLPlot()

    def CSV_export(self):
        """
        Export CSV format

         :returns:  list of all element in dict() function
        """
        d = self.__dict__.copy()
        return list(d.keys()), list(d.values())

# class Axe(DessiaObject):

#     def __init__(self, speed:float,torque:float):
#         self.torque=torque
#         self.speed=speed


class Mesh(DessiaObject):
    """
    Gear mesh definition

    :param z: number of tooth
    :type z: int
    :param db: base diameter
    :type db: float
    :param coefficient_profile_shift: coefficient profile shift of the rack
    :type coefficient_profile_shift: float
    :param rack: class rack define the rack of the mesh
    :type rack: meshe.Rack

    :param material: class material define the gear mesh material
    :type material: meshe.Material
    :param gear_width: gear mesh width
    :type gear_width: float
    :param external_torque: the torque
    :type external_torque: float



    >>> input={'z':13, 'db':40*1e-3, 'coefficient_profile_shift':0.3, 'rack':Rack1
                 coeff_gear_addendum:1, coeff_gear_dedendum:1, coeff_root_radius:1,
                 coeff_circular_tooth_thickness:1}
    >>> mesh1=Mesh(**input) # generation of one gear mesh
    """
    _standalone_in_db = True
    _eq_is_data_eq = True
    _non_serializable_attributes = ['rac', 'reference_point_trochoide', 'reference_point_outside']
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']

    def __init__(self, z: int, db: float, coefficient_profile_shift: float, rack: Rack,
                 material: Material = None,
                 gear_width: float = 1, external_torque: float = None, cycle: float = None,
                 name: str = ''):

        self.rack = rack
        self.gear_param(z, db, coefficient_profile_shift)

        # Definition of default parameters
        self.material = material
        if material is None:
            self.material = hardened_alloy_steel

        self.gear_width = gear_width
        self.external_torque = external_torque
        self.cycle = cycle
        self.reference_point_trochoide = vm.Point2D(0, 0)
        self.reference_point_outside = vm.Point2D(0, 0)
        DessiaObject.__init__(self, name=name)

    def update(self, z, db, coefficient_profile_shift, transverse_pressure_angle_rack,
               coeff_gear_addendum, coeff_gear_dedendum, coeff_root_radius,
               coeff_circular_tooth_thickness, material, gear_width=1):
        """ update of the gear mesh

        :param all: same parameters of this class initialisation

        >>> input={z:14, db:42*1e-3, cp:0.5}
        >>> mesh1.update(**input)
        """
        self.rack.update(self.rack.module, transverse_pressure_angle_rack, coeff_gear_addendum,
                         coeff_gear_dedendum, coeff_root_radius,
                         coeff_circular_tooth_thickness)
        self.gear_param(z, db, coefficient_profile_shift)


        self.gear_width = gear_width

    def update_helix_angle(self, helix_angle, gear_width):
        """ update of the gear mesh

        :param all: same parameters of this class initialisation

        >>> input={z:14, db:42*1e-3, cp:0.5}
        >>> mesh1.update(**input)
        """
        self.rack.update_helix_angle(helix_angle)

        self.gear_width = gear_width

    ### geometry definition

    def gear_param(self, z, db, coefficient_profile_shift):

        self.z = z

        self.db = abs(db)

        self.dff = abs(self.db/math.cos(self.rack.transverse_pressure_angle_0))

        module_rack = abs(self.dff/self.z)
        self.rack.update(module_rack)
        self.coefficient_profile_shift = coefficient_profile_shift

        self.outside_diameter = abs((self.dff
                                     +2*(self.rack.gear_addendum
                                         +self.rack.module*self.coefficient_profile_shift)*(self.z/abs(self.z))))


        self.alpha_outside_diameter = math.acos(self.db/self.outside_diameter)


        self.root_diameter = (self.dff
                              - 2*(self.rack.gear_dedendum
                                   - self.rack.module*self.coefficient_profile_shift)*(self.z/abs(self.z)))


        self.root_diameter_active, self.phi_trochoide = self._root_diameter_active()


        self.alpha_root_diameter_active = math.acos(self.db/self.root_diameter_active)

        self.alpha_pitch_diameter = math.acos(self.db/self.dff)
        self.circular_tooth_thickness = (self.rack.circular_tooth_thickness
                                         +(self.rack.module*self.coefficient_profile_shift
                                           *math.tan(self.rack.transverse_pressure_angle_0)
                                           +self.rack.module*self.coefficient_profile_shift
                                           *math.tan(self.rack.transverse_pressure_angle_0)*self.z/abs(self.z)))
        self.tooth_space = self.rack.transverse_radial_pitch-self.circular_tooth_thickness

        self.outside_active_angle = (2*self.circular_tooth_thickness/self.dff-2
                                     *abs(math.tan(self.alpha_outside_diameter)
                                          -self.alpha_outside_diameter
                                          -math.tan(self.alpha_pitch_diameter)
                                          +self.alpha_pitch_diameter))


        self.base_circular_tooth_thickness = (self.db/2
                                              *(2*self.circular_tooth_thickness/self.dff
                                                +2*(math.tan(self.alpha_pitch_diameter)
                                                    -self.alpha_pitch_diameter)))

        self.root_angle = self.tooth_space/(self.dff/2)-2*(math.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter)
        self.root_gear_angle = self.circular_tooth_thickness/(self.dff/2)+2*(math.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter)

    def gear_section(self, diameter):
        """ Definition of the gear section

        :param diameter: diameter of the gear section calculation
        :type diameter: m

        :results: gear section in m

        >>> gs=mesh1.gear_section(44*1e-3)
        """

        alpha_diameter = math.acos(self.db/diameter)

        theta1 = (math.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter)-(math.tan(alpha_diameter)-alpha_diameter)

        return diameter/2*(2*theta1+abs(self.outside_active_angle)) #TODO

    def _root_diameter_active(self):
        a = self.rack.a*self.z/abs(self.z)
        b = (self.rack.b-(self.rack.module*self.coefficient_profile_shift))*self.z/abs(self.z)

        r = self.dff/2
        # if self.z<0:
        #     r=-r
        phi = -(self.z/abs(self.z))*(a+b*math.tan(math.pi/2-self.rack.transverse_pressure_angle_0))/r

        root_diameter_active = 2*norm(self._trochoide(phi))

        return root_diameter_active, phi

    ### Optimization Method

    def liste_ineq(self):
        """ Compilation method for inequality list used by the optimizer

        :results: vector of data that should be positive
        """
        check, ineq = self.rack.check_rack_viable()

        return ineq

    ### Trace method

    def update_reference_point(self, discret=3):


        self.reference_point_outside = copy.copy(self._outside_trace(0).points[int(len(self._outside_trace(0).points)/2)])
        self._involute_trace(discret, 0, 'T')
        if self.z > 0:

            last_point = self._trochoide_trace(4*discret, 0, 'T').points[-1]
            first_point = self._trochoide_trace(4*discret, 0, 'R').points[0]
            self.reference_point_trochoide = vm.Point2D((first_point[0]-last_point[0])/2+last_point[0], (first_point[1]-last_point[1])/2+last_point[1])




    def contour(self, discret=1, list_number=None):
        """ Definition of the gear contour for volmdlr

        :param discret: number of discretization points on the gear mesh involute
        :param list_number: list of gear tooth to include on the graph

        :results: volmdlr profile

        >>> C1=mesh1.contour(10)
        >>> G1=vm.Contour2D(C1)
        >>> G1.MPLPlot() # generate a plot with matplotlib
        """
        # Analytical tooth profil
        list_number_origin = 1
        if not list_number:
            list_number_origin = 0
            list_number = npy.arange(int(abs(self.z)))
        L = [self._outside_trace(list_number[0])]


        L.append(self._involute_trace(discret, list_number[0], 'T'))
        if self.z > 0:
            L.append(self._trochoide_root_circle_trace(2*discret, list_number[0]))


        L.append(self._involute_trace(discret, list_number[0]+1, 'R'))
        for i in list_number[1::]:
            L.append(self._outside_trace(i))
            L.append(self._involute_trace(discret, i, 'T'))

            if self.z > 0:
                L.append(self._trochoide_root_circle_trace(2*discret, i))

            L.append(self._involute_trace(discret, i+1, 'R'))
        L2 = []
        primitives = []
        for element in L:
                    for point in element.points:
                        if not point in L2:
                                L2.append(point)
                        # else:
                        #     index=len(L2)-1
        # a=vm.edges.LineSegment2D(start=L2[0], end=L2[1]).plot()
        for point_1, point_2 in zip(L2[:-1], L2[1:]):
            line_segment = vm.edges.LineSegment2D(start=point_1, end=point_2)

            primitives.append(line_segment)

        if not list_number_origin:
            primitives.append(vm.edges.LineSegment2D(start=L2[-1], end=L2[0]))



        return primitives


    def contour_circle(self, list_number=None):
        """ Definition of the gear contour for volmdlr

        :param discret: number of discretization points on the gear mesh involute
        :param list_number: list of gear tooth to include on the graph

        :results: volmdlr profile

        >>> C1=mesh1.contour(10)
        >>> G1=vm.Contour2D(C1)
        >>> G1.MPLPlot() # generate a plot with matplotlib
        """
        # Analytical tooth profil
        if not list_number:
            list_number = npy.arange(int(abs(self.z)))
        L = [self._outside_trace(list_number[0])]


        for i in list_number[1::]:
            L.append(self._outside_trace(i))
        arc2d = vm.edges.Arc2D(L[0].points[0], L[int(len(L)/2)].points[int(len(L[int(len(L)/2)].points)/2)], L[-1].points[-1])

        return [arc2d]

    def _involute_trace(self, discret, number, ind='T'):

        if ind == 'T':
            drap = 1
            theta = npy.linspace(math.tan(self.alpha_outside_diameter),
                                 math.tan(self.alpha_root_diameter_active), discret)
        else:
            drap = -1
            theta = npy.linspace(math.tan(self.alpha_root_diameter_active),
                                 math.tan(self.alpha_outside_diameter), discret)



        sol = self._involute(drap*theta)
        x = sol[0]
        y = sol[1]
        p = [vm.Point2D(x[0], y[0])]

        for i in range(1, discret):
            p.append(vm.Point2D(x[i], y[i]))

        ref = primitives2D.OpenedRoundedLineSegments2D(p, {}, False)

        if ind == 'T':
            L = ref.rotation(vm.Point2D(0, 0), -number*2*math.pi/self.z)
            if self.z > 0:
                self.rac = L.points[-1]
        else:
            L = ref.rotation(vm.Point2D(0, 0),
                             self.base_circular_tooth_thickness*2/self.db)
            L = L.rotation(vm.Point2D(0, 0), -number*2*math.pi/self.z)
            if self.z > 0:
                L.points[0] = self.rac
        # if self.z<0:
        #     x=[]
        #     y=[]
        #     for point in L.points:
        #         x.append(point.vector[0])
        #         y.append(point.vector[1])

        # plt.plot(x,y)
        return L

    def _trochoide_trace(self, discret, number, type_flank='T'):
        # Function evolution of the trochoide
        if type_flank == 'T':
            indice_flank = 1
        else:
            indice_flank = -1

        a = indice_flank*self.rack.a  # indice a in the ISO definition of the rack

        phi0 = a/(self.dff/2)

        list_2D = []

        if type_flank == 'R':
            theta = npy.linspace(phi0, indice_flank*self.phi_trochoide, discret)
        else:
            theta = npy.linspace(indice_flank*self.phi_trochoide, phi0, discret)
        for t in theta:
            point = self._trochoide(t, type_flank)
            list_2D.append(vm.Point2D(point[0], point[1]))
        list_2D = primitives2D.OpenedRoundedLineSegments2D(list_2D, {}, False)

        list_2D = list_2D.rotation(vm.Point2D(0, 0), -self.root_angle/2)

        if type_flank == 'T':
            export_2D = list_2D.rotation(vm.Point2D(0, 0), -number*2*math.pi/self.z)
            export_2D.points[0] = self.rac
        else:
            export_2D = list_2D.rotation(vm.Point2D(0, 0), -number*2*math.pi/self.z)
            self.rac = export_2D.points[-1]
        # if self.z<0:
        #     x=[]
        #     y=[]
        #     for point in export_2D.points:
        #         x.append(point.vector[0])
        #         y.append(point.vector[1])

        #     plt.plot(x,y)
        return export_2D


    def _trochoide_root_circle_trace(self, discret, number):
        # Function evolution of the trochoide

        list_2D = []
        a_t = self.rack.a
        phi0_t = a_t*(self.z/abs(self.z))/(self.dff/2)
        trochoide_start = self._trochoide(phi0_t, 'T')

        p1 = vm.Point2D(trochoide_start[0], trochoide_start[1])
        p1 = p1.rotation(vm.Point2D(0, 0), -self.root_angle/2)

        a_r = -1*self.rack.a
        phi0_r = a_r*(self.z/abs(self.z))/(self.dff/2)
        trochoide_end = (self._trochoide(phi0_r, 'R'))

        p2 = vm.Point2D(trochoide_end[0], trochoide_end[1])
        p2 = p1.rotation(vm.Point2D(0, 0), -self.root_angle/2)


        space = p2[1]-p1[1]

        if space > 0:
            theta_t = npy.linspace(1*self.phi_trochoide, phi0_t, discret)
            for t in theta_t:
                point = self._trochoide(t, 'T')
                list_2D.append(vm.Point2D(point[0], point[1]))
            theta_r = npy.linspace(phi0_r, -1*self.phi_trochoide, discret)
            for t in theta_r:
                point = self._trochoide(t, 'R')
                list_2D.append(vm.Point2D(point[0], point[1]))
        else:
            if self.phi_trochoide < 0 and self.rack.a > 0:
                phi0_t = a_t*(self.z/abs(self.z))/((self.dff)/2)
                phi0_r = a_r*(self.z/abs(self.z))/((self.dff)/2)

                theta_t = npy.linspace(1*self.phi_trochoide, phi0_r, discret)
                for t in theta_t:
                    point = self._trochoide(t, 'T')
                    list_2D.append(vm.Point2D(point[0], point[1]))
                p1 = list_2D[-1]
                theta_r = npy.linspace(phi0_t, -1*self.phi_trochoide, discret)
                for t in theta_r:
                    point = self._trochoide(t, 'R')
                    list_2D.append(vm.Point2D(point[0], point[1]))

                p2 = list_2D[-len(theta_r)]
            else:
                phi0_t = a_t*2*(self.z/abs(self.z))/((self.dff)/2)
                phi0_r = a_r*2*(self.z/abs(self.z))/((self.dff)/2)

                theta_t = npy.linspace(1*self.phi_trochoide, phi0_t, discret)
                for t in theta_t:
                    point = self._trochoide(t, 'T')
                    list_2D.append(vm.Point2D(point[0], point[1]))
                p1 = list_2D[-1]
                theta_r = npy.linspace(phi0_r, -1*self.phi_trochoide, discret)
                for t in theta_r:
                    point = self._trochoide(t, 'R')
                    list_2D.append(vm.Point2D(point[0], point[1]))

                p2 = list_2D[-len(theta_r)]


            space_2 = p2[0]-p1[0]
            for i, (point_1, point_2) in enumerate(zip(list_2D[:-1], list_2D[1:])):
                space_2 = point_2[0]-point_1[0]
                # if space_2<0:
                #     print(i)
                #     print(len(list_2D))
                #     print(space_2)
                #     print(space)



        list_2D = primitives2D.OpenedRoundedLineSegments2D(list_2D, {}, False)

        list_2D = list_2D.rotation(vm.Point2D(0, 0), -self.root_angle/2)


        export_2D = list_2D.rotation(vm.Point2D(0, 0), -number*2*math.pi/self.z)
        export_2D.points[0] = self.rac
        self.rac = export_2D.points[-1]

        return export_2D
    def _root_circle_trace(self, number):
        # 2D trace of the connection between the two trochoide

        # on the drive flank
        indice_flank = 1
        a = indice_flank*self.rack.a
        phi0 = a*(self.z/abs(self.z))/(self.dff/2)
        trochoide = self._trochoide(phi0, 'T')
        point = vm.Point2D(trochoide[0], trochoide[1])
        p1 = vm.Point2D(point[0], point[1])
        p1 = p1.rotation(vm.Point2D(0, 0), -self.root_angle/2)

        # on the coast flank

        a = indice_flank*self.rack.a
        phi0 = a*(self.z/abs(self.z))/(self.dff/2)
        trochoide = (self._trochoide(phi0, 'R'))
        p2 = vm.Point2D(trochoide[0], trochoide[1])
        p2 = p2.rotation(vm.Point2D(0, 0), -self.root_angle/2)


        list_2D = primitives2D.OpenedRoundedLineSegments2D([p1, p2], {}, False)

        export_2D = list_2D.rotation(vm.Point2D(0, 0), -number*2*math.pi/self.z)
        # if self.z<0:
        #     x=[]
        #     y=[]
        #     for point in export_2D.points:
        #         x.append(point.vector[0])
        #         y.append(point.vector[1])
        #     plt.plot(x,y)


        return export_2D

    def _outside_trace(self, number):
        # Trace of the top of the gear mesh
        theta4 = math.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter
        p1 = vm.Point2D(self.outside_diameter/2*math.cos(theta4), self.outside_diameter/2*math.sin(theta4))
        p2 = p1.rotation(vm.Point2D(0, 0), self.outside_active_angle/2)
        p3 = p2.rotation(vm.Point2D(0, 0), self.outside_active_angle/2)
        list_2D = primitives2D.OpenedRoundedLineSegments2D([p3, p2, p1], {}, False)

        export_2D = list_2D.rotation(vm.Point2D(0, 0), -number*2*math.pi/self.z)
        x = []
        y = []
        # for point in export_2D.points:
        #     x.append(point.vector[0])
        #     y.append(point.vector[1])
        # plt.plot(x,y)
        return export_2D



    def _involute(self, tan_alpha):
        """ Involute function estimation

        :param tan_alpha: tan of the pressure angle
        :results: tuple of the involute point (x,y) with the origin of the involute position at the point (x=base_diameter/2, y=0) and pressure angle is positive in the direction counter clockwise
        """
        x, y = [], []
        for ta in tan_alpha:
            x.append(self.db/2*math.cos(ta)+self.db/2*ta*math.sin(ta))
            y.append(self.db/2*math.sin(ta)-self.db/2*ta*math.cos(ta))



        return (x, y)

    def _trochoide(self, phi, type_flank='T'):
        # function generation of trochoide point

        if type_flank == 'T':
            indice_flank = 1
        else:
            indice_flank = -1
        a = indice_flank*self.rack.a
        b = self.rack.b-self.rack.module*self.coefficient_profile_shift
        r = self.dff/2

        rho = self.rack.root_radius
        x2 = rho*math.sin(math.atan((a-r*phi)/b)-phi)+a*math.cos(phi)-b*math.sin(phi)+r*(math.sin(phi)-phi*math.cos(phi))

        y2 = -rho*math.cos(math.atan((a-r*phi)/b)-phi)-a*math.sin(phi)-b*math.cos(phi)+r*(math.cos(phi)+phi*math.sin(phi))
        export_point = (y2, x2)
        return export_point

    ### Method for ISO stress

    def _iso_YS(self, s_thickness_iso):
        # stress concentration factor for ISO approach
        rho_f = self.rack.root_radius+self.rack.b**2/(self.dff/2+self.rack.b)
        coeff_ys_iso = 1+0.15*s_thickness_iso/rho_f
        return coeff_ys_iso

    def gear_iso_section(self, angle):
        """ Calculation of the ISO section

        :param angle: pressure angle of the ISO section calculation
        :type angle: radian

        :results: ISO section and ISO height
        """
        a = self.rack.a
        b = self.rack.b-self.rack.module*self.coefficient_profile_shift
        r = self.dff/2
        theta0 = fsolve((lambda theta: a + b*math.tan(theta) + r*(-angle - self.root_angle/2 - theta + math.pi/2)), 0)[0]
        phi0 = (a-b*math.tan(theta0))/r
        pt_iso = self._trochoide(phi0)
        angle0 = math.atan(pt_iso[1]/pt_iso[0])-self.root_angle/2
        angle_iso = self.root_gear_angle-2*angle0
        diameter_iso = 2*norm(pt_iso)
        s_thickness_iso = diameter_iso*math.sin(angle_iso/2)
        h_height_iso = (s_thickness_iso/2)/math.tan(angle)
        return s_thickness_iso, h_height_iso

    ### Export and graph method

    def _plot_data(self, x, heights, y, z, labels=True):
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



    def plot_data(self, centers={}, axis=(1, 0, 0), name=''):

        plot_datas = []



        Gears3D = self.contour(3)



        # L = []
        # L_vector = []
        # i = 0
        # for element in Gears3D:
        #         for point in element.points:
        #             if not point in L_vector:
        #                 # if i==100:
        #                     L_vector.append(point)
        #                     L.append(point)
                        #     i=0
                        # else:
                        #     i+=1

                    # else:
                    #     # print(point.vector)
                    # print(point)
        # L.append(L[0])

        C1 = vm.wires.Contour2D(Gears3D, {})
        # C1 = vm.wires.Contour2D([bezier_curve])
        surface_style = vmp.SurfaceStyle(color_fill=None, opacity=0)
        vect_position_1 = vm.Vector3D(0, 0, 0)
        circle_db = vm.wires.Circle2D(center=vect_position_1, radius=self.db/2)
        circle_dff = vm.wires.Circle2D(center=vect_position_1, radius=self.dff/2)
        circle_root_diameter = vm.wires.Circle2D(center=vect_position_1, radius=self.root_diameter/2)
        circle_outside_diameter = vm.wires.Circle2D(center=vect_position_1, radius=self.outside_diameter/2)
        circle_root_diameter_active = vm.wires.Circle2D(center=vect_position_1, radius=self.root_diameter_active/2)





        edge_style = vmp.EdgeStyle(line_width=2, color_stroke=vmp.colors.BLUE)
        circle_db_plot_data = circle_db.plot_data(edge_style=edge_style)

        text_style = vmp.TextStyle(text_color=vmp.colors.BLUE, text_align_x='center', font_size=1)
        text_db = vmp.Text(comment='db', position_x=0, position_y=self.db/2, text_style=text_style, text_scaling=True)


        edge_style = vmp.EdgeStyle(line_width=2, color_stroke=vmp.colors.RED)
        circle_dff_plot_data = circle_dff.plot_data(edge_style=edge_style)

        text_style = vmp.TextStyle(text_color=vmp.colors.RED, text_align_x='center', font_size=1)
        text_dff = vmp.Text(comment='dff', position_x=0, position_y=self.dff/2, text_style=text_style, text_scaling=True)




        edge_style = vmp.EdgeStyle(line_width=2, color_stroke=vmp.colors.CYAN)
        circle_root_diameter_plot_data = circle_root_diameter.plot_data(edge_style=edge_style)

        text_style = vmp.TextStyle(text_color=vmp.colors.CYAN, text_align_x='center', font_size=1)
        text_root_diameter = vmp.Text(comment='root_diameter', position_x=0,
                                      position_y=self.root_diameter/2, text_style=text_style, text_scaling=True)


        edge_style = vmp.EdgeStyle(line_width=2, color_stroke=vmp.colors.VIOLET)
        circle_outside_diameter_plot_data = circle_outside_diameter.plot_data(edge_style=edge_style)

        text_style = vmp.TextStyle(text_color=vmp.colors.VIOLET, text_align_x='center', font_size=1)
        text_outside_diameter = vmp.Text(comment='outside_diameter', position_x=0,
                                         position_y=self.outside_diameter/2, text_style=text_style, text_scaling=True)


        edge_style = vmp.EdgeStyle(line_width=2, color_stroke=vmp.colors.ORANGE)
        circle_root_diameter_active_plot_data = circle_root_diameter_active.plot_data(edge_style=edge_style)

        text_style = vmp.TextStyle(text_color=vmp.colors.ORANGE, text_align_x='center', font_size=1)
        text_root_diameter_active = vmp.Text(comment='root_diameter_active', position_x=0,
                                             position_y=self.root_diameter_active/2, text_style=text_style, text_scaling=True)





        edge_style = vmp.EdgeStyle(line_width=2, color_stroke=vmp.colors.BLACK)
        C1_plot_data = C1.plot_data(surface_style=surface_style, edge_style=edge_style)

        plot_datas.extend([circle_outside_diameter_plot_data, circle_dff_plot_data,
                           circle_root_diameter_active_plot_data, circle_db_plot_data,
                           C1_plot_data, circle_root_diameter_plot_data, text_db, text_dff,
                           text_root_diameter, text_root_diameter_active, text_outside_diameter])
        return [vmp.PrimitiveGroup(primitives=plot_datas)]

    def z_number_position_gears(self, vector, position, first_gear=True):

        vector_trochoide_gear_1 = vm.Vector2D(self.reference_point_trochoide[0]-position[1],
                                              self.reference_point_trochoide[1]-position[2])


        angle_1 = math.acos(vector.dot(vector_trochoide_gear_1)/(vector.norm()*vector_trochoide_gear_1.norm()))

        sign_angle_1 = npy.sign(vector_trochoide_gear_1.x*vector.y-vector.x*vector_trochoide_gear_1.y)

        if sign_angle_1 < 0:
            if first_gear:
                estimate_z = round(angle_1*self.z/(2*math.pi))
            else:
                estimate_z = round(angle_1*self.z/(2*math.pi))+int(self.z/2)
        else:
            if first_gear:
                estimate_z = self.z-round(angle_1*self.z/(2*math.pi))
            else:
                estimate_z = self.z-round(angle_1*self.z/(2*math.pi))-int(self.z/2)


        return estimate_z





class MeshCombination(DessiaObject):
    """
    Gear Mesh Combination definition
    :param center_distance: List of float defining the center_distance for each connection
    :type center_distance: List[float]
    :param connections: List of tuples of int  defining gear mesh connections [(1,2),(2,3)...] (the int corresponding to the index of the mesh in the List meshes)
    :type connections: List[Tuple[int, int]]
    :param meshes: List of class Mesh objects define each mesh
    :type meshes: List[meshe.Mesh]
    :param safety_factor: Safety factor used for the ISO design
    :type safety_factor: float
    :param transverse_pressure_angle_ini: float defining the first transverse angle of the connections
    :type transverse_pressure_angle_ini: float
    """

    _standalone_in_db = True
    _eq_is_data_eq = True

    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _non_serializable_attributes = ['internal_torque', 'normal_load', 'tangential_load', 'radial_load',
                                    'linear_backlash', 'total_contact_ratio', 'sigma_iso', 'sigma_lim',
                                    'cycle', 'external_torque', 'gear_graph', 'axial_load']

    def __init__(self, center_distance: List[float], connections: List[Tuple[int, int]],
                 meshes: List[Mesh],
                 safety_factor: float = 1.2, transverse_pressure_angle_ini: float = None,
                 name: str = '', infos: str = ''):

        self.center_distance = center_distance
        self.transverse_pressure_angle_ini = transverse_pressure_angle_ini
        self.connections = connections
        for i, connection in enumerate(self.connections):

            if connection.__class__.__name__ == 'list':
               self.connections[i] = (connection[0], connection[1])
        self.meshes = meshes
        self.meshes_dico = []
        self.infos = infos
        for i, meshe in enumerate(meshes):
            self.meshes_dico.append(meshe)



        self.safety_factor = safety_factor

        self.minimum_gear_width = 10e-3
        self.helix_angle = []
        self.external_torque = {}
        self.cycle = {}

        for i, meshe in enumerate(meshes):
            if meshe:
                if meshe.external_torque != None:
                    self.external_torque[i] = meshe.external_torque
                # if meshe.cycle != None:
                self.cycle[i] = meshe.cycle
                self.helix_angle.append(meshe.rack.helix_angle)



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

        transverse_pressure_angle = []

        for num_gear, (num1, num2) in enumerate(self.connections):
            mesh_first = self.meshes_dico[num1]
            mesh_second = self.meshes_dico[num2]
            df_first = 2*self.center_distance[num_gear]*mesh_first.z/mesh_second.z/(1+mesh_first.z/mesh_second.z)
            transverse_pressure_angle.append(math.acos(mesh_first.db/df_first))
        if not transverse_pressure_angle_ini:
            transverse_pressure_angle_ini = transverse_pressure_angle[0]
        self.Z = []
        self.material = []
        for i, mesh in enumerate(meshes):
            self.Z.append(mesh.z)
            self.material.append(mesh.material)

        self.gear_width = []
        self.DB = []
        for  i, mesh in enumerate(meshes):
            self.gear_width.append(mesh.gear_width)
            self.DB.append(mesh.db)

        self.DF, DB_new, self.connections_dfs, self.transverse_pressure_angle\
            = MeshCombination.gear_geometry_parameter(self.Z, transverse_pressure_angle_ini, center_distance,
                                                      connections, gear_graph)
        if len(self.cycle.keys()) < len(list_gear): # the gear mesh optimizer calculate this dictionary
            self.cycle = MeshCombination.cycle_parameter(self.cycle, self.Z, list_gear)

        self.internal_torque, self.normal_load, self.tangential_load, self.radial_load, self.axial_load = MeshCombination.gear_torque(self.Z, self.external_torque, self.DB,
                                                                                                                                      gear_graph, list_gear, connections, self.DF,
                                                                                                                                      self.transverse_pressure_angle, self.helix_angle)

        self.linear_backlash, self.total_contact_ratio, self.transverse_contact_ratio, self.axial_contact_ratio = \
            MeshCombination.gear_contact_ratio_parameter(self.Z, self.DF, self.transverse_pressure_angle,
                                                         center_distance,
                                                         self.meshes_dico, self.connections_dfs,
                                                         connections, self.helix_angle,
                                                         self.gear_width)

        gear_width_new, self.sigma_iso, self.sigma_lim = MeshCombination.gear_width_definition(self.safety_factor,
                                                                                               self.minimum_gear_width,
                                                                                               list_gear, self.tangential_load, self.meshes_dico,
                                                                                               connections,
                                                                                               self.material, self.cycle, self.total_contact_ratio,
                                                                                               self.helix_angle, self.transverse_pressure_angle)
        self.check()
        self._SAP_diameter()
        DessiaObject.__init__(self, name=name)

    def check(self):
        valid = True
        gear_width, _, _ = MeshCombination.gear_width_definition(self.safety_factor,
                                                                 self.minimum_gear_width,
                                                                 self.list_gear, self.tangential_load, self.meshes_dico,
                                                                 self.connections,
                                                                 self.material, self.cycle, self.total_contact_ratio,
                                                                 self.helix_angle,
                                                                 self.transverse_pressure_angle)
        for i, mesh in enumerate(self.meshes):
            if abs(gear_width[i] - mesh.gear_width) > 1e-6:
                valid = False
        self.DF, DB_new, self.connections_dfs, transverse_pressure_angle_new\
            = MeshCombination.gear_geometry_parameter(self.Z, self.transverse_pressure_angle[0], self.center_distance,
                                                      self.connections, self.gear_graph)
        for i, mesh in enumerate(self.meshes):
            if abs(DB_new[i] - mesh.db) > 1e-6:
                valid = False
        for num_gear, connection in enumerate(self.connections):
            if abs(transverse_pressure_angle_new[num_gear] - self.transverse_pressure_angle[num_gear]) > 1e-6:
                valid = False
        return valid

    @classmethod
    def create(cls, Z, center_distance, connections, transverse_pressure_angle_ini,
               coefficient_profile_shift, transverse_pressure_angle_rack,
               coeff_gear_addendum, coeff_gear_dedendum, coeff_root_radius,
               coeff_circular_tooth_thickness, helix_angle, total_contact_ratio_min, transverse_contact_ratio_min,
               percentage_width_difference_pinion_gear, max_width_difference_pinion_gear,
               material=None, external_torque=None, cycle=None,
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


        if material == None:
            material = {list_gear[0]:hardened_alloy_steel}
        for ne in list_gear:
            if ne not in material.keys():
                material[ne] = hardened_alloy_steel

        if external_torque == None:
            external_torque = [{list_gear[0]:100, list_gear[1]:'output'}]

        if cycle == None:
            cycle = {list_gear[0]:1e6}

        DF, DB, connections_dfs, transverse_pressure_angle\
            = cls.gear_geometry_parameter(Z, transverse_pressure_angle_ini, center_distance,
                                          connections, gear_graph)

        if len(cycle.keys()) < len(list_gear): # the gear mesh optimizer calculate this dictionary
            cycle = cls.cycle_parameter(cycle, Z, list_gear)

        internal_torque, normal_load, tangential_load, radial_load, axial_load = cls.gear_torque(Z, external_torque, DB,
                                                                                                 gear_graph, list_gear,
                                                                                                 connections, DF, transverse_pressure_angle, helix_angle)

        meshes = [0]*(len(list_gear))#TODO
        meshes_dico = {}
        for i, num_engr in enumerate(list_gear):
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
            helix_ang = helix_angle[num_engr]
            rack = Rack(transverse_pressure_angle_0=tpa,
                        coeff_gear_addendum=cga, coeff_gear_dedendum=cgd,
                        coeff_root_radius=crr, coeff_circular_tooth_thickness=cct,
                        helix_angle=helix_ang)
            meshes[i] = Mesh(z, db, cp, rack, mat)
            meshes_dico[num_engr] = meshes[i]

        gear_width, sigma_iso, sigma_lim, infos = cls.function_solve_width_definition(safety_factor,
                                                                                      minimum_gear_width,
                                                                                      list_gear, tangential_load, meshes_dico,
                                                                                      connections,
                                                                                      material, cycle, helix_angle,
                                                                                      transverse_pressure_angle, DF,
                                                                                      center_distance, connections_dfs, total_contact_ratio_min,
                                                                                      transverse_contact_ratio_min, percentage_width_difference_pinion_gear,
                                                                                      max_width_difference_pinion_gear)



        linear_backlash, total_contact_ratio, transverse_contact_ratio, axial_contact_ratio = \
            cls.gear_contact_ratio_parameter(Z, DF, transverse_pressure_angle,
                                             center_distance,
                                             meshes_dico, connections_dfs, connections,
                                             helix_angle, gear_width)


        for num_gear in list_gear:
            meshes_dico[num_gear].gear_width = gear_width[num_gear]
            if num_gear in external_torque.keys():
                meshes_dico[num_gear].external_torque = external_torque[num_gear]
            if num_gear in cycle.keys():
                meshes_dico[num_gear].cycle = cycle[num_gear]

        mesh_combination = cls(center_distance, connections, meshes, safety_factor, transverse_pressure_angle_ini, infos=infos)
        return mesh_combination

    def update(self, Z, center_distance, connections, transverse_pressure_angle_ini,
               coefficient_profile_shift,
               transverse_pressure_angle_rack, coeff_gear_addendum,
               coeff_gear_dedendum, coeff_root_radius, coeff_circular_tooth_thickness,
               material, internal_torque, cycle, safety_factor, total_contact_ratio_min, transverse_contact_ratio_min,
               percentage_width_difference_pinion_gear,
               max_width_difference_pinion_gear):
        """ update of the gear mesh assembly

        :param all: same parameters of this class initialisation

        >>> Z={1:13,2:46,4:38}
        >>> center_distance=[0.118,0.125]
        >>> mesh_assembly1.update(Z=Z,center_distance=center_distance)
        """
        self.center_distance = center_distance
        self.transverse_pressure_angle_ini = transverse_pressure_angle_ini
        self.DF, self.DB, self.connections_dfs, self.transverse_pressure_angle\
            = MeshCombination.gear_geometry_parameter(Z, transverse_pressure_angle_ini, center_distance,
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
            self.meshes_dico[num_engr].update(z, db, cp, tpa, cga, cgd,
                                              crr, cct, mat)

        self.gear_width, self.sigma_iso, self.sigma_lim, self.infos = self.function_solve_width_definition(self.safety_factor,
                                                                                                           self.minimum_gear_width,
                                                                                                           self.list_gear, self.tangential_load, self.meshes_dico,
                                                                                                           self.connections,
                                                                                                           self.material, self.cycle, self.helix_angle,
                                                                                                           self.transverse_pressure_angle, self.DF, self.center_distance,
                                                                                                           self.connections_dfs, total_contact_ratio_min, transverse_contact_ratio_min,
                                                                                                           percentage_width_difference_pinion_gear, max_width_difference_pinion_gear)


        self.linear_backlash, self.total_contact_ratio, self.transverse_contact_ratio, self.axial_contact_ratio = \
            MeshCombination.gear_contact_ratio_parameter(Z, self.DF, self.transverse_pressure_angle,
                                                         center_distance,
                                                         self.meshes_dico, self.connections_dfs,
                                                         connections, self.helix_angle,
                                                         self.gear_width)

        self._SAP_diameter()

    def update_helix_angle(self, helix_angle, total_contact_ratio_min, transverse_contact_ratio_min, percentage_width_difference_pinion_gear,
                           max_width_difference_pinion_gear):
        """ update of the gear mesh assembly

        :param all: same parameters of this class initialisation

        >>> Z={1:13,2:46,4:38}
        >>> center_distance=[0.118,0.125]
        >>> mesh_assembly1.update(Z=Z,center_distance=center_distance)
        """
        if helix_angle:
            self.helix_angle = helix_angle

        self.internal_torque, self.normal_load, self.tangential_load,\
        self.radial_load, self.axial_load = MeshCombination.gear_torque(self.Z, self.external_torque, self.DB,
                                                                        self.gear_graph, self.list_gear,
                                                                        self.connections, self.DF,
                                                                        self.transverse_pressure_angle, self.helix_angle)

        self.gear_width, sigma_iso, sigma_lim, self.infos = self.function_solve_width_definition(self.safety_factor,
                                                                                                 self.minimum_gear_width,
                                                                                                 self.list_gear, self.tangential_load, self.meshes_dico,
                                                                                                 self.connections,
                                                                                                 self.material, self.cycle, self.helix_angle,
                                                                                                 self.transverse_pressure_angle, self.DF,
                                                                                                 self.center_distance, self.connections_dfs,
                                                                                                 total_contact_ratio_min, transverse_contact_ratio_min,
                                                                                                 percentage_width_difference_pinion_gear, max_width_difference_pinion_gear)

        for num_engr in self.list_gear:
            helix_angle = self.helix_angle[num_engr]
            gear_width = self.gear_width[num_engr]
            self.meshes_dico[num_engr].update_helix_angle(helix_angle, gear_width)

        self.linear_backlash, self.total_contact_ratio, self.transverse_contact_ratio, self.axial_contact_ratio = \
            MeshCombination.gear_contact_ratio_parameter(self.Z, self.DF, self.transverse_pressure_angle,
                                                         self.center_distance,
                                                         self.meshes_dico, self.connections_dfs, self.connections,
                                                         self.helix_angle, self.gear_width)

    ### Optimization Method
    def check_minimum_backlash(self, backlash_min=2*1e-4):
        """ Define constraint and functional for the optimizer on backlash

        :param backlash_min: maximum backlash available
        :results:
            * check is a boolean (True if 0<backlash<backlash_min)
            * list_ineq a list of element that should be positive for the optimizer
            * obj is a functional on the backlash used for the optimizer
        """
        list_ineq = [] # liste of value to evaluate backlash
        obj = 0

        for lb in self.linear_backlash:
            list_ineq.append(lb) # backlash > 0
            list_ineq.append(backlash_min-lb) # backlash < backlash_min so (backlash_min-backlash)>0
            obj += 10*(lb-backlash_min)**2
        check = False
        if min(list_ineq) > 0:
            check = True
        return check, list_ineq, obj

    def check_total_contact_ratio(self, total_contact_ratio_min, transverse_contact_ratio_min):
        """ Define constraint and functional for the optimizer on radial contact ratio

        :param transverse_contact_ratio_min: minimum radial contact ratio available
        :results:
            * check is a boolean (True if transverse_contact_ratio_min<total_contact_ratio)
            * list_ineq a list of element that should be positive for the optimizer
            * obj is a functional on the backlash used for the optimizer
        """
        list_ineq = []
        obj = 0
        for num_mesh, (eng1, eng2) in enumerate(self.connections):
            rca = self.total_contact_ratio[num_mesh]
            list_ineq.append(rca-total_contact_ratio_min[(eng1, eng2)])
            transverse_contact_ratio = self.transverse_contact_ratio[num_mesh]
            list_ineq.append(transverse_contact_ratio-transverse_contact_ratio_min[(eng1, eng2)])
            if rca > total_contact_ratio_min[(eng1, eng2)]:
                obj += 0.001*(rca-total_contact_ratio_min[(eng1, eng2)])
            else:
                obj += 1000*(total_contact_ratio_min[(eng1, eng2)]-rca)
            if transverse_contact_ratio > transverse_contact_ratio_min[(eng1, eng2)]:
                obj += 0.001*(transverse_contact_ratio-transverse_contact_ratio_min[(eng1, eng2)])
            else:
                obj += 1000*(transverse_contact_ratio_min[(eng1, eng2)]-transverse_contact_ratio)
        check = False
        if min(list_ineq) > 0:
            check = True
        return check, list_ineq, obj

    def liste_ineq(self, total_contact_ratio_min, transverse_contact_ratio_min):
        """ Compilation method for inequality list used by the optimizer

        :results: vector of data that should be positive
        """
        _, ineq, _ = self.check_minimum_backlash(4*1e-4)

        _, list_ineq, _ = self.check_total_contact_ratio(total_contact_ratio_min, transverse_contact_ratio_min)

        ineq.extend(list_ineq)

        for mesh in self.meshes:
            list_ineq = mesh.liste_ineq()
            ineq.extend(list_ineq)


        return ineq

    def functional(self, total_contact_ratio_min=1, transverse_contact_ratio_min=1):
        """ Compilation method for a part of the functional used by the optimizer

        :results: scalar add to the global functional of the optimizer
        """
        check1, ineq1, obj1 = self.check_minimum_backlash(4*1e-4)
        check2, ineq2, obj2 = self.check_total_contact_ratio(total_contact_ratio_min, transverse_contact_ratio_min)
        obj = obj1 + obj2
        return obj

    ### Method gear mesh calculation
    @classmethod
    def gear_contact_ratio_parameter(cls, Z, DF, transverse_pressure_angle,
                                     center_distance,
                                     meshes, connections_dfs, connections, helix_angle, gear_width):


        linear_backlash = []
        total_contact_ratio = []
        axial_contact_ratio = []
        transverse_contact_ratio = []
        total_contact_ratio = []
        for engr1, engr2 in connections_dfs:
            if (engr1, engr2) in connections:
                num_mesh = connections.index((engr1, engr2))
                engr1_position = 0
                engr2_position = 1

            elif (engr2, engr1) in connections:
                num_mesh = connections.index((engr2, engr1))
                engr1_position = 1
                engr2_position = 0


            else:
                raise RuntimeError
            circular_tooth_thickness1 = meshes[engr1].gear_section(DF[num_mesh][engr1_position])
            circular_tooth_thickness2 = meshes[engr2].gear_section(DF[num_mesh][engr2_position])

            transverse_radial_pitch1 = math.pi*DF[num_mesh][engr1_position]/abs(meshes[engr1].z)
            space_width1 = transverse_radial_pitch1-circular_tooth_thickness1
            space_width2 = transverse_radial_pitch1-circular_tooth_thickness2

            linear_backlash.append(min(space_width1-circular_tooth_thickness2, space_width2-circular_tooth_thickness1))
            transverse_pressure_angle1 = transverse_pressure_angle[num_mesh]
            center_distance1 = abs(center_distance[num_mesh])
            axial_contact_ratio.append(abs(math.sin(helix_angle[engr1])*gear_width[engr1]/(math.pi*meshes[engr1].rack.module)))



            transverse_contact_ratio.append(1/2.*(math.sqrt(meshes[engr1].outside_diameter**2
                                                            - meshes[engr1].db**2)
                                                  + math.sqrt(meshes[engr2].outside_diameter**2
                                                              - meshes[engr2].db**2)
                                                  - 2*center_distance1*math.sin(transverse_pressure_angle1))
                                            /(transverse_radial_pitch1*math.cos(transverse_pressure_angle1)))#TODO


            total_contact_ratio.append(axial_contact_ratio[-1]+transverse_contact_ratio[-1])



        return linear_backlash, total_contact_ratio, transverse_contact_ratio, axial_contact_ratio

    @classmethod
    def gear_geometry_parameter(cls, Z, transverse_pressure_angle_ini, center_distance, connections, gear_graph):
        # Construction of pitch and base diameter
        DF = [0]*len(connections)
        db = {}


        dict_transverse_pressure_angle = {0: transverse_pressure_angle_ini}
        connections_dfs = list(nx.edge_dfs(gear_graph,
                                           [connections[0][0], connections[0][1]]))

        for num_dfs, ((engr1, engr2), cd) in enumerate(zip(connections_dfs, center_distance)):

            if (engr1, engr2) in connections:
                num_mesh = connections.index((engr1, engr2))
                engr1_position = 0
                engr2_position = 1
            else:

                num_mesh = connections.index((engr2, engr1))
                engr1_position = 1
                engr2_position = 0
            Z1 = Z[engr1]
            Z2 = Z[engr2]
            if Z1 < 0:
                DF1 = abs(2*cd*abs(Z1/Z2)/(-1+abs(Z1/Z2)))
                DF2 = abs(-2*cd+DF1)

            elif Z2 < 0:

                DF1 = abs(2*cd*abs(Z1/Z2)/(1-abs(Z1/Z2)))
                DF2 = abs(2*cd+DF1)
            else:
                DF1 = abs(2*cd*Z1/Z2/(1+Z1/Z2))
                DF2 = abs(2*cd-DF1)

            DF[num_mesh] = [0]*2
            DF[num_mesh][engr1_position] = DF1
            DF[num_mesh][engr2_position] = DF2
            if num_mesh == 0:
                db1 = float(DF1*math.cos(transverse_pressure_angle_ini))
                db2 = float(DF2*math.cos(transverse_pressure_angle_ini))
            else:
                db1 = db[engr1]
                try:

                    dict_transverse_pressure_angle[num_mesh] = math.acos(db1/DF1)
                except:
                    print('Error Diameter DB {}, DF {}, Z1 {}, Z2 {}, pa {}'.format(db1, DF1, Z1, Z2, transverse_pressure_angle_ini))
                    raise ValidGearDiameterError()
                db2 = DF2*math.cos(dict_transverse_pressure_angle[num_mesh])
            db[engr1] = db1
            db[engr2] = db2


        transverse_pressure_angle = []
        for num_mesh in sorted(dict_transverse_pressure_angle.keys()):
            tpa = dict_transverse_pressure_angle[num_mesh]
            transverse_pressure_angle.append(tpa)


        return DF, db, connections_dfs, transverse_pressure_angle


    def _SAP_diameter(self):

        self.SAP = [0]*len(self.connections)
        self.SAP_diameter = [0]*len(self.connections)
        for num_mesh, (engr1, engr2) in enumerate(self.connections):


            self.SAP[num_mesh] = [0]*2

            self.SAP_diameter[num_mesh] = [0]*2

            teta_om = math.acos(self.DB[engr2]/self.meshes[engr2].outside_diameter)
            self.SAP[num_mesh][0] = 180*((self.Z[engr1]+self.Z[engr2])*math.tan(self.transverse_pressure_angle[num_mesh])-self.Z[engr2]*math.tan(teta_om))/(math.pi*self.Z[engr1])

            self.SAP_diameter[num_mesh][0] = self.DB[engr1]*math.sqrt((math.pi*self.SAP[num_mesh][0]/180)**2+1)


            teta_om = math.acos(self.DB[engr1]/self.meshes[engr1].outside_diameter)
            self.SAP[num_mesh][1] = 180*((self.Z[engr1]+self.Z[engr2])*math.tan(self.transverse_pressure_angle[num_mesh])-self.Z[engr1]*math.tan(teta_om))/(math.pi*self.Z[engr2])

            self.SAP_diameter[num_mesh][1] = self.DB[engr2]*math.sqrt((math.pi*self.SAP[num_mesh][1]/180)**2+1)

    @classmethod
    def gear_torque(cls, Z, external_torque, db, gear_graph, list_gear,
                    connections, DF, transverse_pressure_angle, helix_angle):
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
        if 'output' in external_torque.values():
            for num_gear, tq in external_torque.items():
                if tq == 'output':
                    node_output = num_gear
            torque_graph_dfs = list(nx.dfs_edges(gear_graph, node_output))


            order_torque_calculation = [(eng2, eng1) for (eng1, eng2) in torque_graph_dfs[::-1]]

            # calculation torque distribution
            temp_torque = {}
            for eng1 in list_gear:
                temp_torque[eng1] = 0
            for num_mesh_tq, (eng1, eng2) in enumerate(order_torque_calculation):
                if eng1 in external_torque.keys():
                    temp_torque[eng1] += external_torque[eng1]
                temp_torque[eng2] += -temp_torque[eng1]*Z[eng2]/float(Z[eng1])
            dic_torque = {}
            for num_mesh_tq, (eng1, eng2) in enumerate(order_torque_calculation):
                dic_torque[(eng1, eng2)] = temp_torque[eng1]

        else:#TODO

            external_torque[list(external_torque.keys())[0]] = 'output'
            for num_gear, tq in external_torque.items():
                if tq == 'output':
                    node_output = num_gear
            torque_graph_dfs = list(nx.dfs_edges(gear_graph, node_output))


            order_torque_calculation = [(eng2, eng1) for (eng1, eng2) in torque_graph_dfs[::-1]]

            # calculation torque distribution
            temp_torque = {}
            for eng1 in list_gear:
                temp_torque[eng1] = 0
            for num_mesh_tq, (eng1, eng2) in enumerate(order_torque_calculation):
                if eng1 in external_torque.keys():
                    temp_torque[eng1] += external_torque[eng1]
                temp_torque[eng2] += -temp_torque[eng1]*Z[eng2]/float(Z[eng1])
            dic_torque = {}
            for num_mesh_tq, (eng1, eng2) in enumerate(order_torque_calculation):
                dic_torque[(eng1, eng2)] = temp_torque[eng1]

        normal_load = {}
        tangential_load = {}
        radial_load = {}
        axial_load = {}
        for num_mesh, (eng1, eng2) in enumerate(connections):
            # if 'output' not in external_torque.values():
            #     dic_torque=external_torque
            try:

                tq = dic_torque[(eng1, eng2)]
                eng1_position = 0
                eng2_position = 1
                normal_load[num_mesh] = abs(tq)*2/(db[eng1])
                tangential_load[num_mesh] = abs(tq)*2/(DF[num_mesh][eng1_position])

                axial_load[num_mesh] = tangential_load[num_mesh]*math.tan(helix_angle[eng1])
                radial_load[num_mesh] = math.tan(transverse_pressure_angle[num_mesh])*tangential_load[num_mesh]/math.cos(helix_angle[eng1])
            except:
                tq = dic_torque[(eng2, eng1)]
                eng1_position = 1
                eng2_position = 0

                normal_load[num_mesh] = abs(tq)*2/(db[eng2])
                tangential_load[num_mesh] = abs(tq)*2/(DF[num_mesh][eng2_position])

                axial_load[num_mesh] = tangential_load[num_mesh]*math.tan(helix_angle[eng2])

                radial_load[num_mesh] = math.tan(transverse_pressure_angle[num_mesh])*tangential_load[num_mesh]/math.cos(helix_angle[eng2])
        return dic_torque, normal_load, tangential_load, radial_load, axial_load

    @classmethod
    def cycle_parameter(cls, cycle, Z, list_gear):
        """ Calculation of the gear mesh cycle

        :param Z: dictionary define the number of teeth {node1:Z1, node2:Z2, node3:Z3 ...}
        :param cycle: Dictionary defining the number of cycle for one node {node3: number_cycle3}

        :results: dictionary define the number of cycle for each gear mesh {node1:cycle1, node2:cycle2, node3:cycle3 ...}
        """

        eng_init = list(cycle.keys())[0]
        for eng in list_gear:
            if eng not in cycle.keys():
                cycle[eng] = cycle[eng_init]*Z[eng_init]/float(Z[eng])
        return cycle

    @classmethod
    def gear_width_definition(cls, safety_factor, minimum_gear_width,
                              list_gear, tangential_load, meshes, connections,
                              material, cycle, total_contact_ratio, helix_angle,
                              transverse_pressure_angle):
        """ Calculation of the gear width

        :param safety_factor: Safety factor used for the ISO design

        :results:
            * **gear_width** - dictionary define the gear mesh width {node1 : gw1, node2 : gw2, node3 : gw3 ...}
            * **sigma_iso** - dictionary define the ISO stress {mesh1 : {node1 sig_iso1: , node2 : sig_iso2_1}, mesh2 : {node2 : sig_iso2_2, node3 : sig_iso3} ...}
            * **sigma_lim** - dictionary define the limit material stress {mesh1 : {node1 sig_lim1: , node2 : sig_lim2}, mesh2 : {node2 : sig_lim2, node3 : sig_lim3} ...}

        in this function, we define the gear width for each gear mesh to respect sig_lim = sig_iso for each gear mesh
        """
        coeff_yf_iso = cls._coeff_YF_iso(connections, meshes, transverse_pressure_angle)
        coeff_ye_iso = cls._coeff_YE_iso(connections, total_contact_ratio)
        coeff_yb_iso = cls._coeff_YB_iso(connections, material, helix_angle)

        sigma_lim = cls.sigma_material_iso(safety_factor, connections,
                                           material, cycle, meshes)

        gear_width = {}
        for eng in list_gear:
            gear_width[eng] = minimum_gear_width

        for num_mesh, (eng1, eng2) in enumerate(connections):

            gear_width1 = abs((tangential_load[num_mesh]
                               / (sigma_lim[num_mesh][eng1]
                                  * meshes[eng1].rack.module))
                              *coeff_yf_iso[num_mesh][eng1]
                              *coeff_ye_iso[num_mesh]
                              *coeff_yb_iso[num_mesh][eng1])

            gear_width2 = abs((tangential_load[num_mesh]
                               /(sigma_lim[num_mesh][eng2]
                                 *meshes[eng2].rack.module))
                              *coeff_yf_iso[num_mesh][eng2]
                              *coeff_ye_iso[num_mesh]
                              *coeff_yb_iso[num_mesh][eng2])

            gear_width_set = max(gear_width1, gear_width2)
            gear_width[eng1] = max(gear_width[eng1], gear_width_set)
            gear_width[eng2] = max(gear_width[eng2], gear_width_set)

        sigma_iso = sigma_lim
        return gear_width, sigma_iso, sigma_lim


    @classmethod
    def function_solve_width_definition(cls, safety_factor, minimum_gear_width,
                                        list_gear, tangential_load, meshes, connections,
                                        material, cycle, helix_angle,
                                        transverse_pressure_angle, DF, center_distance, connections_dfs,
                                        total_contact_ratio_min, transverse_contact_ratio_min, percentage_width_difference_pinion_gear,
                                        max_width_difference_pinion_gear):
        """ Calculation of the gear width

        :param safety_factor: Safety factor used for the ISO design

        :results:
            * **gear_width** - dictionary define the gear mesh width {node1 : gw1, node2 : gw2, node3 : gw3 ...}
            * **sigma_iso** - dictionary define the ISO stress {mesh1 : {node1 sig_iso1: , node2 : sig_iso2_1}, mesh2 : {node2 : sig_iso2_2, node3 : sig_iso3} ...}
            * **sigma_lim** - dictionary define the limit material stress {mesh1 : {node1 sig_lim1: , node2 : sig_lim2}, mesh2 : {node2 : sig_lim2, node3 : sig_lim3} ...}

        in this function, we define the gear width for each gear mesh to respect sig_lim = sig_iso for each gear mesh
        """


        infos = ''
        coeff_yf_iso = cls._coeff_YF_iso(connections, meshes, transverse_pressure_angle)

        coeff_yb_iso = cls._coeff_YB_iso(connections, material, helix_angle)

        sigma_lim = cls.sigma_material_iso(safety_factor, connections,
                                           material, cycle, meshes)
        max_gear_width = minimum_gear_width
        gear_width = {}
        for eng in list_gear:
            gear_width[eng] = minimum_gear_width

        def f(x, tangential_load, sigma_lim, meshes, coeff_yf_iso,
              coeff_yb_iso, DF, connection, transverse_pressure_angle,
              center_distance, connections_dfs, helix_angle, total_contact_ratio_min, transverse_contact_ratio_min):

            contact_ratio = cls.gear_contact_ratio_parameter(Z=[meshes[0].z, meshes[1].z], DF=DF, transverse_pressure_angle=transverse_pressure_angle,
                                                             center_distance=center_distance, meshes=meshes, connections_dfs=connections_dfs,
                                                             connections=connection, helix_angle=helix_angle, gear_width=[abs(x[0]), abs(x[1])])


            f_eng1 = abs(abs(x[0])-abs((tangential_load[0]
                                        / (sigma_lim[0][0]
                                           * meshes[0].rack.module))
                                       *coeff_yf_iso[0][0]
                                       *1/contact_ratio[1][0]
                                       *coeff_yb_iso[0][0]))

            f_eng2 = abs(abs(x[1])-abs((tangential_load[0]
                                        / (sigma_lim[0][1]
                                           * meshes[1].rack.module))
                                       *coeff_yf_iso[0][1]
                                       *1/contact_ratio[1][0]
                                       *coeff_yb_iso[0][1]))


            # print(x)
            # print(contact_ratio[1][0])
            # print(abs((tangential_load[0]
            #                         / (sigma_lim[0][1]
            #                            * meshes[1].rack.module))
            #                        *coeff_yf_iso[0][1]
            #                        *1/contact_ratio[1][0]
            #                        *coeff_yb_iso[0][1]))


            if contact_ratio[1][0] < total_contact_ratio_min:
                f_contact_ratio_min = abs(total_contact_ratio_min-transverse_contact_ratio_min-contact_ratio[3][0])
            else:
                f_contact_ratio_min = 0



            return [f_eng1+f_contact_ratio_min, f_eng2+f_contact_ratio_min]

        for num_mesh, (eng1, eng2) in enumerate(connections):


            xs = fsolve(f, npy.zeros(2), args=([tangential_load[num_mesh]],
                                               [[sigma_lim[num_mesh][eng1], sigma_lim[num_mesh][eng2]]],
                                               [meshes[eng1], meshes[eng2]],
                                               [[coeff_yf_iso[num_mesh][eng1], coeff_yf_iso[num_mesh][eng2]]],
                                               [[coeff_yb_iso[num_mesh][eng1], coeff_yb_iso[num_mesh][eng2]]],
                                               [[DF[num_mesh][0], DF[num_mesh][1]]],
                                               [(0, 1)],
                                               [transverse_pressure_angle[num_mesh]],
                                               [center_distance[num_mesh]],
                                               [(0, 1)],
                                               [helix_angle[eng1], helix_angle[eng2]],
                                               total_contact_ratio_min[(eng1, eng2)],
                                               transverse_contact_ratio_min[(eng1, eng2)]), full_output=0)



            gear_width1 = abs(xs[0])
            gear_width2 = abs(xs[1])

            contact_ratio = cls.gear_contact_ratio_parameter(Z=[meshes[eng1].z, meshes[eng2].z], DF=[[DF[num_mesh][0], DF[num_mesh][1]]], transverse_pressure_angle=[transverse_pressure_angle[num_mesh]],
                                                             center_distance=[center_distance[num_mesh]], meshes=[meshes[eng1], meshes[eng2]], connections_dfs=[(0, 1)],
                                                             connections=[(0, 1)], helix_angle=[helix_angle[eng1], helix_angle[eng2]], gear_width=[gear_width1, gear_width2])

            width_torque_gear_1 = abs((tangential_load[num_mesh]
                                       / (sigma_lim[num_mesh][eng1]
                                          * meshes[eng1].rack.module))
                                      *coeff_yf_iso[num_mesh][eng1]
                                      *1/contact_ratio[1][0]
                                      *coeff_yb_iso[num_mesh][eng1])


            width_torque_gear_2 = abs((tangential_load[num_mesh]
                                       / (sigma_lim[num_mesh][eng2]
                                          * meshes[eng2].rack.module))
                                      *coeff_yf_iso[num_mesh][eng2]
                                      *1/contact_ratio[1][0]
                                      *coeff_yb_iso[num_mesh][eng2])

            axial_contact_ratio_min = total_contact_ratio_min[(eng1, eng2)]-contact_ratio[2][0]
            if helix_angle[eng1]:
                width_contact_ratio = abs(axial_contact_ratio_min*math.pi*meshes[eng1].rack.module/(math.sin(helix_angle[eng1])))
            else:
                width_contact_ratio = width_torque_gear_1
            gear_width_set = max(gear_width1, gear_width2, minimum_gear_width)
            
            if len(connections) < 2:
                if meshes[eng1].z < meshes[eng2].z:
                    if gear_width_set*percentage_width_difference_pinion_gear[(eng1, eng2)] > max_width_difference_pinion_gear[(eng1, eng2)]:
                        gear_width[eng1] = gear_width_set + max_width_difference_pinion_gear[(eng1, eng2)]
                        gear_width[eng2] = gear_width_set
                    else:
                        gear_width[eng1] = gear_width_set*(1+percentage_width_difference_pinion_gear[(eng1, eng2)])
                        gear_width[eng2] = gear_width_set
                else:
                    if gear_width_set*percentage_width_difference_pinion_gear[(eng1, eng2)] > max_width_difference_pinion_gear[(eng1, eng2)]:
                        gear_width[eng1] = gear_width_set
                        gear_width[eng2] = gear_width_set + max_width_difference_pinion_gear[(eng1, eng2)]
                    else:
                        gear_width[eng1] = gear_width_set
                        gear_width[eng2] = gear_width_set*(1+percentage_width_difference_pinion_gear[(eng1, eng2)])
            else:
                if gear_width_set > max_gear_width:
                    max_gear_width = gear_width_set
                    for eng in gear_width.keys():
                        gear_width[eng] = gear_width_set
                
                    
                    

            angle = 30./180.*math.pi
            s_thickness_iso_1, h_height_iso_1 = meshes[eng1].gear_iso_section(angle)
            coeff_ys_iso_gear_1 = meshes[eng1]._iso_YS(s_thickness_iso_1)
            s_thickness_iso_2, h_height_iso_2 = meshes[eng2].gear_iso_section(angle)
            coeff_ys_iso_gear_2 = meshes[eng2]._iso_YS(s_thickness_iso_2)

            infos += 'For sizing a mesh, we need to determinate the minimum width to support the tangential load for the 2 gears and the minimum width of having the good contact ratio. \n\n'
            infos += 'To Begin, some general infos on the mesh that will serve us in all the different calculs: \n\n'
            infos += '|Module|Safety Factor|Tangential Load|Helix Angle|' + '\n'
            infos += '|:--------:|:-------------:|:---------------:|:---------------:|' + '\n'
            infos += '|'+str(round(meshes[eng1].rack.module*10**3, 3))+' mm'+'|'+str(round(safety_factor, 3))+'|'+str(round(tangential_load[num_mesh], 3))+' N'+'|'+str(round(helix_angle[eng1]*180/math.pi, 3))+' ° |' + '\n\n'
            infos += 'To calculate the minimum width for the tangential load, we need for that to use some factors and parameters which depend on the material and the tooth forms of the gear:' +'\n\n \n\n'
            infos += '|Gear|Coeff YS (stress concentration factor)|Coeff YB( helix angle factor )|Coeff YF (form factor)|Sigma Lim|' + '\n'
            infos += '|:--------:|:-------------:|:---------------:|:---------------:|:---------------:|' + '\n'
            infos += '| 1 |'+str(round(coeff_ys_iso_gear_1, 3))+'|'+str(round(coeff_yb_iso[num_mesh][eng1], 3))+'|'+str(round(coeff_yf_iso[num_mesh][eng1], 3))+'|' + str(round(coeff_ys_iso_gear_1*safety_factor*sigma_lim[num_mesh][eng1]*10**-6, 3))+' MPa'+'|\n'
            infos += '| 2 |'+str(round(coeff_ys_iso_gear_2, 3))+'|'+str(round(coeff_yb_iso[num_mesh][eng2], 3))+'|'+str(round(coeff_yf_iso[num_mesh][eng2], 3))+'|' + str(round(coeff_ys_iso_gear_1*safety_factor*sigma_lim[num_mesh][eng2]*10**-6, 3))+' MPa'+'|\n\n'

            infos += 'The minimum width to support the tangential load for the first gear  is **'+str(round(width_torque_gear_1*10**3, 3))+'** mm '
            infos += 'and for the second gear is **'+str(round(width_torque_gear_2*10**3, 3))+'** mm'+ '\n\n'
            infos += 'To calculate the minimum width for the contact ratio, we need to have the  axial contact ratio minimum require.\n\n'
            infos += 'The total contact ratio minimum is **' +str(round(total_contact_ratio_min[(eng1, eng2)], 3))+ '** and the actual transverse contact ratio is '+str(round(contact_ratio[2][0], 3))+ '\n\n'
            infos += 'So the axial contact ratio minimum is **' +str(round(axial_contact_ratio_min, 3))+ '**\n\n'
            infos += 'To have this requirement we need to have a width equal at **'+str(round(width_contact_ratio*10**3, 3))+'** mm'+'\n\n'
            infos += '|Width Tangential Load minimum|Width Contact Ratio minimum|' + '\n'
            infos += '|:--------:|:-------------:|' + '\n'
            infos += '|'+str(round(max([width_torque_gear_1*10**3, width_torque_gear_2*10**3]), 3))+' mm'+'|'+str(round(width_contact_ratio*10**3, 3))+'|' + '\n\n'
            if contact_ratio[3][0] > axial_contact_ratio_min*1.05:
                infos += 'The tangential load condition is most restritive \n\n'
            else:
                infos += 'The contact ratio condition is most restritive \n\n'
            infos += 'Percentage Difference Between Pinion And Gear: ' +str(round(percentage_width_difference_pinion_gear[(eng1, eng2)]*100, 3))+' %'+'\n\n'
            infos += 'Max Difference Between Pinion And Gear: ' +str(round(max_width_difference_pinion_gear[(eng1, eng2)]*1000, 3))+' mm'+'\n\n'
            infos += '####Solution: ' +'\n\n'

            infos += 'Width Gear 1: '+str(round(gear_width[eng1]*10**3, 2))+' mm'+'\n\n'+\
                    'Width Gear 2: '+str(round(gear_width[eng2]*10**3, 2))+' mm'+'\n\n'+\
                    'Total Contact Ratio: '+str(round(contact_ratio[1][0], 3))+'\n\n'+\
                    'Transverse Contact Ratio: '+str(round(contact_ratio[2][0], 3))+'\n\n'+\
                    'Axial Contact Ratio: '+str(round(contact_ratio[3][0], 3))+ '\n\n \n\n'







        sigma_iso = sigma_lim

        return gear_width, sigma_iso, sigma_lim, infos

    @classmethod
    def sigma_material_iso(cls, safety_factor, connections, material, cycle,
                           meshes):
        """ Calculation of the material limit stress

        :param safety_factor: Safety factor used for the ISO design

        :results:
            * **sigma_lim** - dictionary define the limit material stress {mesh1 : {node1 sig_lim1: , node2 : sig_lim2}, mesh2 : {node2 : sig_lim2, node3 : sig_lim3} ...}

        in this function, we use the FunCoeff function of the Material class to interpolate the material parameters
        """
        angle = 30./180.*math.pi
        sigma_lim = {}
        for num_mesh, (eng1, eng2) in enumerate(connections):
            sigma_lim[num_mesh] = {}

            matrice_wohler = material[eng1].data_wohler_curve
            matrice_material = material[eng1].data_gear_material
            sgla = material[eng1].FunCoeff(cycle[eng1], npy.array(matrice_wohler.data), matrice_wohler.x, matrice_wohler.y)
            # sgl1 = material[eng1].FunCoeff(sgla,npy.array(matrice_material.data), matrice_material.x, matrice_material.y)
            s_thickness_iso_1, h_height_iso_1 = meshes[eng1].gear_iso_section(angle)
            coeff_ys_iso = meshes[eng1]._iso_YS(s_thickness_iso_1)

            sigma_lim[num_mesh][eng1] = float((sgla/(safety_factor*coeff_ys_iso))*10**6)

            matrice_wohler = material[eng2].data_wohler_curve
            matrice_material = material[eng2].data_gear_material
            sglb = material[eng2].FunCoeff(cycle[eng2], npy.array(matrice_wohler.data), matrice_wohler.x, matrice_wohler.y)
            # sgl2 = material[eng2].FunCoeff(sglb, npy.array(matrice_material.data), matrice_material.x, matrice_material.y)
            s_thickness_iso_2, h_height_iso_2 = meshes[eng2].gear_iso_section(angle)
            coeff_ys_iso = meshes[eng2]._iso_YS(s_thickness_iso_2)

            sigma_lim[num_mesh][eng2] = float((sglb/(safety_factor*coeff_ys_iso))*10**6)

        return sigma_lim

    @classmethod
    def _coeff_YF_iso(cls, connections, meshes, transverse_pressure_angle):
        # shape factor for ISO stress calculation
        angle = 30./180.*math.pi
        coeff_yf_iso = {}
        for num_mesh, (eng1, eng2) in enumerate(connections):
            coeff_yf_iso[num_mesh] = {}
            s_thickness_iso_1, h_height_iso_1 = meshes[eng1].gear_iso_section(angle)
            s_thickness_iso_2, h_height_iso_2 = meshes[eng2].gear_iso_section(angle)
            coeff_yf_iso[num_mesh][eng1] = ((6*(h_height_iso_1/meshes[eng1].rack.module)*math.cos(transverse_pressure_angle[num_mesh]))
                                            /((s_thickness_iso_1/meshes[eng1].rack.module)**2
                                              *math.cos(meshes[eng1].rack.transverse_pressure_angle_0)))
            coeff_yf_iso[num_mesh][eng2] = ((6*(h_height_iso_2/meshes[eng2].rack.module)*math.cos(transverse_pressure_angle[num_mesh]))
                                            /((s_thickness_iso_2/meshes[eng2].rack.module)**2
                                              *math.cos(meshes[eng2].rack.transverse_pressure_angle_0)))


        return coeff_yf_iso

    @classmethod
    def _coeff_YE_iso(cls, connections, total_contact_ratio):
        #  radial contact ratio factor for ISO stress calculation
        coeff_ye_iso = []
        for ne, eng in enumerate(connections):
            coeff_ye_iso.append(1/total_contact_ratio[ne])

        return coeff_ye_iso

    @classmethod
    def _coeff_YB_iso(cls, connections, material, helix_angle):
        # gear widht factor impact for ISO stress calculation
        coeff_yb_iso = {}
        for num_mesh, (eng1, eng2) in enumerate(connections):
            coeff_yb_iso[num_mesh] = {}

            matrice_YB = material[eng1].data_coeff_YB_Iso
            coeff_yb_iso[num_mesh][eng1] = material[eng1].FunCoeff(helix_angle[eng1]*180/3.14, npy.array(matrice_YB.data), matrice_YB.x, matrice_YB.y)
            matrice_YB = material[eng2].data_coeff_YB_Iso
            coeff_yb_iso[num_mesh][eng2] = material[eng2].FunCoeff(helix_angle[eng2]*180/3.14, npy.array(matrice_YB.data), matrice_YB.x, matrice_YB.y)

        return coeff_yb_iso

    ### Function graph and export

    def gear_rotate(self, gear_index, list_gear, list_center, list_rot):
        """ Displacement of the volmdlr gear profile (rotation and translation)

        :param list_gear: list of volmdlr contour [meshes.Contour, meshes.Contour ...], each contour is centered on the origin
        :param list_center: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
        :param list_rot: list of rotation for each gear mesh [node1 : rot1, node2 : rot2 ...]

        :results: list of volmdlr component
        """
        export = []

        for (index, i, center, k) in zip(gear_index, list_gear, list_center, list_rot):
            model_export = []



            for m in i:
                center = vm.Point2D(center[0], center[1])
                model_trans = m.translation(center)
                model_trans_rot = model_trans.rotation(center, k)

                model_export.append(model_trans_rot)
            export.append(model_export)
        return export

    def gear_rotate_reference_point(self, gear_index, list_center, list_trans, list_rot):
        """ Displacement of the volmdlr gear profile (rotation and translation)

        :param list_gear: list of volmdlr contour [meshes.Contour, meshes.Contour ...], each contour is centered on the origin
        :param list_center: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
        :param list_rot: list of rotation for each gear mesh [node1 : rot1, node2 : rot2 ...]

        :results: list of volmdlr component
        """
        export = []

        for (index, center, trans, k) in zip(gear_index, list_center, list_trans, list_rot):

            position = [self.meshes_dico[index].reference_point_trochoide[0]+trans[0],
                        self.meshes_dico[index].reference_point_trochoide[1]+trans[1]]


            u = [position[0] - center[0], position[1] - center[1]]
            v2x = math.cos(k)*u[0] - math.sin(k)*u[1] + center[0]
            v2y = math.sin(k)*u[0] + math.cos(k)*u[1] + center[1]

            self.meshes_dico[index].reference_point_trochoide.x = v2x
            self.meshes_dico[index].reference_point_trochoide.y = v2y


            position = [self.meshes_dico[index].reference_point_outside[0]+trans[0],
                        self.meshes_dico[index].reference_point_outside[1]+trans[1]]

            u = [position[0] - center[0], position[1] - center[1]]
            v2x = math.cos(k)*u[0] - math.sin(k)*u[1] + center[0]
            v2y = math.sin(k)*u[0] + math.cos(k)*u[1] + center[1]
            self.meshes_dico[index].reference_point_outside.x = v2x
            self.meshes_dico[index].reference_point_outside.y = v2y



        return export

    def initial_position(self, positions, liste_eng=()):
        """ Calculation of the rotation for two gear mesh to initiate the contact

        :param list_gear: list of volmdlr contour [meshes.Contour, meshes.Contour ...], each contour is centered on the origin
        :param list_center: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
        :param list_rot: list of rotation for each gear mesh [node1 : rot1, node2 : rot2 ...]

        :results: list of volmdlr component
        """


        reference_point_trochoide_gear_1 = self.meshes_dico[liste_eng[0]].reference_point_trochoide



        reference_point_outside_gear_2 = self.meshes_dico[liste_eng[1]].reference_point_outside


        center_distance = vm.Vector2D((positions[1][1]-positions[0][1]), (positions[1][2]-positions[0][2]))

        vector_trochoide_gear_1 = vm.Vector2D(reference_point_trochoide_gear_1[0]-positions[0][1],
                                              reference_point_trochoide_gear_1[1]-positions[0][2])

        vector_outside_gear_2 = vm.Vector2D(reference_point_outside_gear_2[0]-positions[1][1],
                                            reference_point_outside_gear_2[1]-positions[1][2])

        angle_1 = math.acos(center_distance.dot(vector_trochoide_gear_1)/(center_distance.norm()*vector_trochoide_gear_1.norm()))
        angle_2 = math.acos(center_distance.dot(vector_outside_gear_2)/(center_distance.norm()*vector_outside_gear_2.norm()))

        sign_angle_1 = npy.sign(vector_trochoide_gear_1.x*center_distance.y-center_distance.x*vector_trochoide_gear_1.y)
        sign_angle_2 = npy.sign(vector_outside_gear_2.x*center_distance.y-center_distance.x*vector_outside_gear_2.y)



        rotation_gear_2 = sign_angle_2*(angle_2-math.pi)+sign_angle_1*angle_1*self.meshes_dico[liste_eng[0]].z/self.meshes_dico[liste_eng[1]].z
        return 0, rotation_gear_2

    # TODO: use volmdlr Vector and points

    def plot_data(self, centers={}, axis=(1, 0, 0), name=''):
        """
        2D mesh combination visualization

        :param centers: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
        :param axis: direction of gear mesh rotation
        :returns: List of Primitives groups for the data_plot

        """
        x = vm.Vector3D(axis[0], axis[1], axis[2])

        # y = x.RandomUnitNormalVector()
        # y= vm.Vector3D((0,1,0))
        y = x.deterministic_unit_normal_vector()

        z = x.cross(y)
        if len(centers) == 0:
            centers = {}
            center_var = self.pos_axis({self.list_gear[0]:[0, 0]})

            for engr_num in center_var.keys():
                centers[engr_num] = [0, center_var[engr_num][0], center_var[engr_num][1]]
        else:
            center_var = {}
            for engr_num in centers.keys():

                center_var[engr_num] = npy.dot(centers[engr_num], (x[0], x[1], x[2]))*x+npy.dot(centers[engr_num], (y[0], y[1], y[2]))*y+npy.dot(centers[engr_num], (z[0], z[1], z[2]))*z
                center_var[engr_num] = (center_var[engr_num][0], center_var[engr_num][1], center_var[engr_num][2])
            centers = center_var



        Gears3D = {}
        Struct = []
        Rotation = {}
        plot_datas = []
        # plt.figure()
        # plt.axis('equal')
        list_rot = [0]*len(self.meshes)
        for set_pos_dfs, (eng1, eng2) in enumerate(self.connections_dfs):

            position1 = centers[eng1]
            position2 = centers[eng2]
            if set_pos_dfs == 0:
                Gears3D[eng1] = self.meshes_dico[eng1].contour(3)

            Gears3D[eng2] = self.meshes_dico[eng2].contour(3)

            vect_position_1 = vm.Vector3D(position1[0], position1[1], position1[2])
            vect_position_2 = vm.Vector3D(position2[0], position2[1], position2[2])
            self.meshes_dico[eng1].update_reference_point()
            self.meshes_dico[eng2].update_reference_point()
            self.gear_rotate_reference_point([eng1, eng2],
                                             [([vect_position_1.dot(y), vect_position_1.dot(z)]),
                                              ([vect_position_2.dot(y), vect_position_2.dot(z)])],
                                             [([vect_position_1.dot(y), vect_position_1.dot(z)]),
                                              ([vect_position_2.dot(y), vect_position_2.dot(z)])],
                                             list_rot=[list_rot[eng1], list_rot[eng2]])

            if (eng1, eng2) in self.connections:
                set_pos = self.connections.index((eng1, eng2))
                rot_gear_2 = self.initial_position([position1, position2], (eng1, eng2))
                eng1_position = 0
                eng2_position = 1

            elif (eng2, eng1) in self.connections:
                set_pos = self.connections.index((eng2, eng1))
                rot_gear_2 = self.initial_position([position2, position1], (eng2, eng1))
                eng2_position = 0
                eng1_position = 1

            Rotation[set_pos] = {}

            Struct.append(vm.wires.Circle2D(vm.Point2D(position1[0], position1[1]), self.DF[set_pos][eng1_position]/2.))
            Struct.append(vm.wires.Circle2D(vm.Point2D(position2[0], position2[1]), self.DF[set_pos][eng2_position]/2.))

            list_rot[eng1] += rot_gear_2[0]

            list_rot[eng2] += rot_gear_2[1]


            if set_pos_dfs == 0:
                Gears3D_Rotate = self.gear_rotate([eng1, eng2],
                                                  [Gears3D[eng1], Gears3D[eng2]],
                                                  [([vect_position_1.dot(y), vect_position_1.dot(z)]),
                                                   ([vect_position_2.dot(y), vect_position_2.dot(z)])],
                                                  list_rot=[list_rot[eng1], list_rot[eng2]])
            else:
                Gears3D_Rotate = self.gear_rotate([eng2],
                                                  [Gears3D[eng2]],
                                                  [([vect_position_2.dot(y), vect_position_2.dot(z)])],
                                                  list_rot=[list_rot[eng2]])








            # L = []
            # L_vector = []
            # i = 0
            # for element in Gears3D_Rotate[0]:
            #         for point in element.points:
            #            if not point in L_vector:
            #                # if i==100:
            #                    L_vector.append(point)
            #                    L.append(point)
                           #     i=0
                           # else:
                           #     i+=1

                       # else:
                       #     # print(point.vector)
                       # print(point)
            # L.append(L[0])


            # vmp.plot([C1.plot_data('contour')])
            # L2 = []
            # L2_vector = []
            # i = 0
            # for element in Gears3D_Rotate[1]:
            #         for point in element.points:
            #            if not point in L2_vector:
            #                # if i==100:
            #                    L2_vector.append(point)
            #                    L2.append(point)
            #                    i = 0
                           # else:
                           #     i+=1

            # L2.append(L2[0])

            # L2=set(L2)





            if set_pos_dfs == 0:
                C1 = vm.wires.Contour2D(Gears3D_Rotate[0], {})
                C2 = vm.wires.Contour2D(Gears3D_Rotate[1], {})

                circle_DF = vm.wires.Circle2D(center=vect_position_1, radius=self.DF[0][eng1_position]/2)

                circle_SAP_diameter = vm.wires.Circle2D(center=vect_position_1, radius=self.SAP_diameter[0][eng1_position]/2)


                edge_style = vmp.EdgeStyle(line_width=2, color_stroke=vmp.colors.GREEN)
                circle_DF_plot_data = circle_DF.plot_data(edge_style=edge_style)

                text_style = vmp.TextStyle(text_color=vmp.colors.GREEN, text_align_x='center', font_size=0.7)
                text_DF = vmp.Text(comment='DF', position_x=0, position_y=self.DF[0][eng1_position]/2, text_style=text_style)




                edge_style = vmp.EdgeStyle(line_width=2, color_stroke=vmp.colors.ROSE)
                circle_SAP_diameter_plot_data = circle_SAP_diameter.plot_data(edge_style=edge_style)



                text_style = vmp.TextStyle(text_color=vmp.colors.ROSE, text_align_x='center', font_size=0.7)
                text_SAP_diameter = vmp.Text(comment='SAP_diameter', position_x=0,
                                             position_y=self.SAP_diameter[0][eng1_position]/2, text_style=text_style, text_scaling=True)

                surface_style = vmp.SurfaceStyle(color_fill=vmp.colors.WHITE)
                edge_style = vmp.EdgeStyle(line_width=2, color_stroke=vmp.colors.BLACK)
                C1_plot_data = C1.plot_data(surface_style=surface_style, edge_style=edge_style)
                C2_plot_data = C2.plot_data(surface_style=surface_style, edge_style=edge_style)
                plot_datas.extend([circle_DF_plot_data, circle_SAP_diameter_plot_data, C1_plot_data, C2_plot_data, text_SAP_diameter, text_DF])


            else:
                C2 = vm.wires.Contour2D(Gears3D_Rotate[0], {})

                C2_plot_data = C2.plot_data(surface_style=surface_style, edge_style=edge_style)
                plot_datas.extend([C2_plot_data])

        return [vmp.PrimitiveGroup(primitives=plot_datas)]






    # def volmdlr_primitives_2(self, centers={}, axis=(1, 0, 0), name='', z_number=10):
    #     """ Generation of the 3D volume for all the gear mesh


    #     :param center: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
    #     :param axis: direction of gear mesh rotation


    #     :results: list of 3D volmdlr component
    #     """
    #     primitives=[]

    #     x = vm.Vector3D(axis[0],axis[1],axis[2])
    #     # y = x.RandomUnitNormalVector()
    #     # y= vm.Vector3D((0,1,0))
    #     y = x.deterministic_unit_normal_vector()

    #     z = x.cross(y)
    #     if len(centers) == 0:
    #         centers = {}
    #         center_var = self.pos_axis({self.list_gear[0]:[0, 0]})

    #         for engr_num in center_var.keys():
    #             centers[engr_num]=[0, center_var[engr_num][0], center_var[engr_num][1]]
    #     else:
    #         center_var = {}
    #         for engr_num in centers.keys():

    #             center_var[engr_num] = npy.dot(centers[engr_num],(x[0],x[1],x[2]))*x+npy.dot(centers[engr_num],(y[0],y[1],y[2]))*y+npy.dot(centers[engr_num],(z[0],z[1],z[2]))*z
    #             center_var[engr_num] = (center_var[engr_num][0],center_var[engr_num][1],center_var[engr_num][2])
    #         centers = center_var



    #     Gears3D = {}
    #     Struct = []
    #     Rotation = {}
    #     plot_datas = []
    #     # plt.figure()
    #     # plt.axis('equal')

    #     for set_pos_dfs, (eng1, eng2) in enumerate(self.connections_dfs):

    #         position1 = centers[eng1]
    #         position2 = centers[eng2]
    #         if set_pos_dfs == 0:
    #             Gears3D[eng1] = self.meshes_dico[eng1].contour(3)

    #         Gears3D[eng2] = self.meshes_dico[eng2].contour(3)


    #         if (eng1, eng2) in self.connections:
    #             set_pos = self.connections.index((eng1, eng2))
    #             rot_gear_2= self.initial_position([position1,position2], (eng1, eng2))
    #             eng1_position=0
    #             eng2_position=1


    #         elif (eng2, eng1) in self.connections:
    #             set_pos = self.connections.index((eng2, eng1))
    #             rot_gear_2= self.initial_position([position2,position1], (eng2, eng1))
    #             eng1_position=1
    #             eng2_position=0


    #         Rotation[set_pos] = {}
    #         Struct.append(vm.wires.Circle2D(vm.Point2D(position1[0],position1[1]),self.DF[set_pos][eng1_position]/2.))
    #         Struct.append(vm.wires.Circle2D(vm.Point2D(position2[0],position2[1]),self.DF[set_pos][eng2_position]/2.))






    #         vect_position_1 = vm.Vector3D(position1[0],position1[1],position1[2])
    #         vect_position_2 = vm.Vector3D(position2[0],position2[1],position2[2])
    #         Gears3D_Rotate = self.gear_rotate([eng1, eng2],
    #                                             [Gears3D[eng1],Gears3D[eng2]],
    #                                             [([vect_position_1.dot(y),vect_position_1.dot(z)]),([vect_position_2.dot(y),vect_position_2.dot(z)])],
    #                                             list_rot=[rot_gear_2[0],rot_gear_2[1]])



    #         # for Gears in Gears3D_Rotate:
    #         #     for element in Gears:
    #         #         for point in element.points:
    #         #             x2.append(point.vector[0])
    #         #             y2.append(point.vector[1])
    #         # plt.plot(x2,y2)


    #         L = []
    #         L_vector = []
    #         i=0
    #         for element in Gears3D_Rotate[0]:
    #                 for point in element.points:
    #                     if not point in L_vector:
    #                         # if i==100:
    #                             L_vector.append(point)
    #                             L.append(point)
    #                         #     i=0
    #                         # else:
    #                         #     i+=1

    #                     # else:
    #                     #     # print(point.vector)
    #                     # print(point)
    #         # L.append(L[0])
    #         bezier_curve=vm.edges.BezierCurve2D(3, L)
    #         C1 = vm.wires.ClosedPolygon2D(L,{})
    #         # C1 = vm.wires.Contour2D([bezier_curve])
    #         # vmp.plot([C1.plot_data('contour')])
    #         L2 = []
    #         L2_vector = []
    #         i=0
    #         for element in Gears3D_Rotate[1]:
    #                 for point in element.points:
    #                     if not point in L2_vector:
    #                         # if i==100:
    #                             L2_vector.append(point)
    #                             L2.append(point)
    #                             i=0
    #                         # else:
    #                         #     i+=1

    #         # L2.append(L2[0])

    #         # L2=set(L2)
    #         C2 = vm.wires.ClosedPolygon2D(L2, {})

    #     #     C1=vm.Contour2D(Gears3D_Rotate[0])
    #     #     # print(Gears3D_Rotate[0])
    #     #     C2=vm.Contour2D(Gears3D_Rotate[1])

    #         extrusion_vector1 = (self.gear_width[eng1]*x)
    #         extrusion_vector2 = (self.gear_width[eng2]*x)


    #         if set_pos_dfs == 0:
    #             vect_x = -0.5*self.gear_width[eng1]*x + x.dot(vm.Vector3D(centers[eng1][0],centers[eng1][1],centers[eng1][2]))*x

    #             if self.Z[eng1] < 0:
    #                 vect_center = vm.Vector3D(centers[eng1][0],centers[eng1][1],centers[eng1][2])
    #                 vector=vm.Vector2D(vect_center.dot(y),vect_center.dot(z))
    #                 circle = vm.wires.Circle2D(vm.Point2D(vector[0],vector[1]),(self.DB[eng1]*1.3)/2)
    #                 t1 = primitives3D.ExtrudedProfile(vm.Vector3D(vect_x[0],vect_x[1],vect_x[2]), y, z,circle , [C1], vm.Vector3D(extrusion_vector1[0],extrusion_vector1[1],extrusion_vector1[2]))
    #             else:
    #                 try:
    #                     t1 = primitives3D.ExtrudedProfile(vm.Vector3D(vect_x[0],vect_x[1],vect_x[2]), y, z, C1, [], vm.Vector3D(extrusion_vector1[0],extrusion_vector1[1],extrusion_vector1[2]))
    #                 except ZeroDivisionError or ValueError:
    #                     vector=vm.Vector2D(vect_center.dot(y),vect_center.dot(z))
    #                     circle = vm.wires.Circle2D(vm.Point2D(vector[0],vector[1]),(self.DB[eng1])/2)
    #                     t1 = primitives3D.ExtrudedProfile(vm.Vector3D(vect_x[0],vect_x[1],vect_x[2]), y, z, circle, [], vm.Vector3D(extrusion_vector1[0],extrusion_vector1[1],extrusion_vector1[2]))

    #             primitives.append(t1)
    #         vect_x = -0.5*self.gear_width[eng2]*x + x.dot(vm.Vector3D(centers[eng2][0],centers[eng2][1],centers[eng2][2]))*x

    #         if self.Z[eng2] < 0:
    #                 vect_center = vm.Vector3D(centers[eng2][0],centers[eng2][1],centers[eng2][2])

    #                 circle = vm.wires.Circle2D(vm.Point2D(vect_center.dot(y),vect_center.dot(z)),(self.DB[eng2]*1.3)/2)

    #                 t2 = primitives3D.ExtrudedProfile(vm.Vector3D(vect_x[0],vect_x[1],vect_x[2]), y, z,circle , [C2], vm.Vector3D(extrusion_vector2[0],extrusion_vector2[1],extrusion_vector2[2]))

    #         else:
    #             t2 = primitives3D.ExtrudedProfile(vm.Vector3D(vect_x[0],vect_x[1],vect_x[2]), y, z, C2, [], vm.Vector3D(extrusion_vector2[0],extrusion_vector2[1],extrusion_vector2[2]))


    #         primitives.append(t2)





    #     return primitives



    def volmdlr_primitives(self, centers={}, axis=(1, 0, 0), name='', z_number=10):
        """ Generation of the 3D volume for all the gear mesh

        :param center: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
        :param axis: direction of gear mesh rotation

        :results: list of 3D volmdlr component
        """
        primitives = []

        x = vm.Vector3D(axis[0], axis[1], axis[2])
        y = x.deterministic_unit_normal_vector()

        z = x.cross(y)
        if len(centers) == 0:
            centers = {}
            center_var = self.pos_axis({self.list_gear[0]:[0, 0]})

            for engr_num in center_var.keys():
                centers[engr_num] = [0, center_var[engr_num][0], center_var[engr_num][1]]
        else:
            center_var = {}
            for engr_num in centers.keys():

                center_var[engr_num] = npy.dot(centers[engr_num], (x[0], x[1], x[2]))*x+\
                npy.dot(centers[engr_num], (y[0], y[1], y[2]))*y+\
                npy.dot(centers[engr_num], (z[0], z[1], z[2]))*z
                center_var[engr_num] = (center_var[engr_num][0], center_var[engr_num][1], center_var[engr_num][2])
            centers = center_var



        Gears3D = {}
        Struct = []
        Rotation = {}
        plot_datas = []
        # plt.figure()
        # plt.axis('equal')
        list_z_gear = [0]*len(self.meshes)
        primitive_plot_data = []
        list_rot = [0]*len(self.meshes)

        for set_pos_dfs, (eng1, eng2) in enumerate(self.connections_dfs):
            if not list_z_gear[eng1]:
                list_z_gear[eng1] = []
            if not list_z_gear[eng2]:
                list_z_gear[eng2] = []


            position1 = centers[eng1]
            position2 = centers[eng2]

            self.meshes_dico[eng1].update_reference_point()
            self.meshes_dico[eng2].update_reference_point()
            vect_position_1 = vm.Vector3D(position1[0], position1[1], position1[2])
            vect_position_2 = vm.Vector3D(position2[0], position2[1], position2[2])
            x = vm.Vector3D(axis[0], axis[1], axis[2])
            y = x.deterministic_unit_normal_vector()

            z = x.cross(y)

            self.gear_rotate_reference_point([eng1, eng2],
                                             [([vect_position_1.dot(y), vect_position_1.dot(z)]), ([vect_position_2.dot(y), vect_position_2.dot(z)])],
                                             [([vect_position_1.dot(y), vect_position_1.dot(z)]), ([vect_position_2.dot(y), vect_position_2.dot(z)])],
                                             list_rot=[list_rot[eng1], list_rot[eng2]])



            if (eng1, eng2) in self.connections:
                set_pos = self.connections.index((eng1, eng2))
                rot_gear_2 = list(self.initial_position([position1, position2], (eng1, eng2)))
                eng1_position = 0
                eng2_position = 1

            elif (eng2, eng1) in self.connections:
                set_pos = self.connections.index((eng2, eng1))
                rot_gear_2 = list(self.initial_position([position2, position1], (eng2, eng1)))
                eng1_position = 1
                eng2_position = 0


            list_rot[eng1] += rot_gear_2[0]


            list_rot[eng2] += rot_gear_2[1]








            vect_position_1 = vm.Vector3D(position1[0], position1[1], position1[2])
            vect_position_2 = vm.Vector3D(position2[0], position2[1], position2[2])

            vect_center_1 = vm.Point2D(position1[1], position1[2])
            vect_center_2 = vm.Point2D(position2[1], position2[2])

            x = vm.Vector3D(axis[0], axis[1], axis[2])
            y = x.deterministic_unit_normal_vector()

            z = x.cross(y)
            Gears3D_Rotate = self.gear_rotate_reference_point([eng1, eng2],
                                                              [([vect_position_1.dot(y), vect_position_1.dot(z)]), ([vect_position_2.dot(y), vect_position_2.dot(z)])],
                                                              [([0, 0]), ([0, 0])],
                                                              list_rot=[rot_gear_2[0], rot_gear_2[1]])

            center_distance = vm.Vector2D((position2[1]-position1[1]), (position2[2]-position1[2]))

            estimate_z_1 = self.meshes_dico[eng1].z_number_position_gears(vector=center_distance, position=position1)
            estimate_z_2 = self.meshes_dico[eng2].z_number_position_gears(vector=center_distance, position=position2, first_gear=False)




            z_num = z_number/2-1

            if z_number > self.meshes_dico[eng1].z:
                z_num = self.meshes_dico[eng1].z/2


            list_number_1 = list(npy.arange(estimate_z_1-z_num, estimate_z_1+z_num+1))

            for z in list_number_1:
                 if z < 0:
                     z += self.meshes_dico[eng1].z
                 if z >= self.meshes_dico[eng1].z:
                     z -= self.meshes_dico[eng1].z
                 if z not in list_z_gear[eng1]:
                     list_z_gear[eng1].append(z)
            # model_trans_rot_1=[]
            # Gears3D[eng1] = self.meshes_dico[eng1].contour(3,list_number=list_number_1)

            # for element in Gears3D[eng1]:
            #     model_trans_1 = element.translation(vect_center_1)
            #     model_trans_rot_1.append(model_trans_1.rotation(vect_center_1, rot_gear_2[0]))


            z_num = z_number/2-1
            if z_number > self.meshes_dico[eng2].z:
                z_num = self.meshes_dico[eng2].z/2


            list_number_2 = list(npy.arange(estimate_z_2-z_num, estimate_z_2+z_num+1))

            for z in list_number_2:
                 if z < 0:
                     z += self.meshes_dico[eng2].z

                 if z >= self.meshes_dico[eng2].z:
                     z -= self.meshes_dico[eng2].z
                 if z not in list_z_gear[eng2]:
                     list_z_gear[eng2].append(z)


        for num_gear, l in enumerate(list_z_gear):

            L_total = []
            list_z = sorted(l)
            l_contour_z = []

            if list_z[0] > 0:
                list_z_circle = list(npy.arange(0, list_z[0]+1))

                L_total.extend(self.meshes_dico[num_gear].contour_circle(list_number=list_z_circle))
            for i, z in enumerate(list_z):
                if z-list_z[i-1] > 1:
                    L_total.extend(self.meshes_dico[num_gear].contour(3, list_number=l_contour_z))
                    l_contour_z = []
                    list_z_circle = list(npy.arange(list_z[i-1]+1, z+2))

                    L_total.extend(self.meshes_dico[num_gear].contour_circle(list_number=list_z_circle))
                else:
                    l_contour_z.append(z)
            if l_contour_z:
                L_total.extend(self.meshes_dico[num_gear].contour(3, list_number=l_contour_z))
            if list_z[-1] < self.meshes_dico[num_gear].z-1:
                list_z_circle = list(npy.arange(list_z[-1]+1, self.meshes_dico[num_gear].z+1))
                L_total.extend(self.meshes_dico[num_gear].contour_circle(list_number=list_z_circle))

            model_trans_rot_1 = []
            position1 = centers[num_gear]
            vect_center_1 = vm.Point2D(position1[1], position1[2])

            for element in L_total:
                model_trans_1 = element.translation(vect_center_1)
                model_trans_rot_1.append(model_trans_1.rotation(vect_center_1, list_rot[num_gear]))
            # L = []

            C1 = vm.wires.Contour2D(primitives=model_trans_rot_1)


            # print(C1.primitives)
            # for element in model_trans_rot_1:
            #         for point in element.points:
            #            if not point in L:
            #                    L.append(point)

            # C1 = vm.wires.ClosedPolygon2D(L, {})

            x = vm.Vector3D(axis[0], axis[1], axis[2])
            y = x.deterministic_unit_normal_vector()

            z = x.cross(y)




            extrusion_vector1 = (self.gear_width[num_gear]*x)



            vect_x = -0.5*self.gear_width[num_gear]*x + x.dot(vm.Vector3D(centers[num_gear][0], centers[num_gear][1], centers[num_gear][2]))*x

            if self.Z[num_gear] < 0:
                vect_center = vm.Vector3D(centers[num_gear][0], centers[num_gear][1], centers[num_gear][2])
                vector = vm.Vector2D(vect_center.dot(y), vect_center.dot(z))
                circle = vm.wires.Circle2D(vm.Point2D(vector[0], vector[1]), (self.DB[num_gear]*1.3)/2)
                t1 = primitives3D.ExtrudedProfile(vm.Vector3D(vect_x[0], vect_x[1], vect_x[2]), y, z, circle, [C1], vm.Vector3D(extrusion_vector1[0], extrusion_vector1[1], extrusion_vector1[2]))


            else:
                try:

                    t1 = primitives3D.ExtrudedProfile(vm.Vector3D(vect_x[0], vect_x[1], vect_x[2]), y, z, C1, [], vm.Vector3D(extrusion_vector1[0], extrusion_vector1[1], extrusion_vector1[2]))


                except ZeroDivisionError or ValueError:
                    vector = vm.Vector2D(vect_center.dot(y), vect_center.dot(z))
                    circle = vm.wires.Circle2D(vm.Point2D(vector[0], vector[1]), (self.DB[num_gear])/2)
                    t1 = primitives3D.ExtrudedProfile(vm.Vector3D(vect_x[0], vect_x[1], vect_x[2]), y, z, circle, [], vm.Vector3D(extrusion_vector1[0], extrusion_vector1[1], extrusion_vector1[2]))

            primitives.append(t1)








        return primitives

    def mass(self):
        """
        Estimation of gear mesh mass

        :results: mass of all gear mesh
        """
        DF = {}
        for i, (ic1, ic2) in enumerate(self.connections):
            DF[ic1] = self.DF[i][0]
            DF[ic2] = self.DF[i][1]

        mass = 0.
        for i, df in DF.items():
            mass += self.gear_width[i] * self.material[i].volumic_mass* math.pi * (0.5*DF[i])**2

        return mass

    # Waiting for meshes to know how to plot themselves
#    def plot_data(self, x, heights, ys, zs, labels = True):
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
#            t, a = mesh.plot_data(x, heights, ys[imesh], zs[imesh], labels)
#            transversal_plot_data.extend(t)
#            axial_plot_data.extend(a)
#
#        # Ploting axial because mesh doesn't know its width
#
#        return axial_plot_data, transversal_plot_data

    def FreeCADExport(self, fcstd_filepath, centers={}, axis=(1, 0, 0),
                      python_path='python', path_lib_freecad='/usr/lib/freecad/lib',
                      export_types=['fcstd']):
        """ Export 3D volume to FreeCAD

        :param file_path: file path for the freecad file
        :param center: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
        :param axis: direction of gear mesh rotation

        :results: export of a FreeCAD file
        """
        model = self.volume_model(centers, axis)
        model.FreeCADExport(fcstd_filepath, python_path, path_lib_freecad, export_types)

    # TODO change this function to make it like plot_data in PWT: output is dict of geometrical shapes
    def SVGExport(self, name, position):
        """ Export SVG graph of all gear mesh

        :param name: name of the svg file
        :param position: dictionary define some center position {node2 : [0,0], node4 : [0.12,0] ..}

        :results: SVG graph

        in the position dictionary, you have to be coherent with the center position

            * for exemple, if the center-distance of the mesh1 (node1, node2) is 0.117 m you can define position such as:

                * {node1 : [0,0], node2 : [0.117,0]}
                * {node1 : [0,0]}
        """
        x_opt = position
        TG = {}
        L1 = []
        Struct = []
        Rot = {}
        for num, en in enumerate(self.connections_dfs):
            position1 = (x_opt[en[0]][0], x_opt[en[0]][1])
            position2 = (x_opt[en[1]][0], x_opt[en[1]][1])
            #tuple1 et 2 correspondent a la position des centres
            ne = self.connections.index(en)
            Rot[ne] = {}
            if num == 0:
                TG[en[0]] = self.meshes[en[0]].Contour(5)
            if (en[0], en[1]) in self.connections:
                eng1_position = 0
                eng2_position = 1

            elif (en[1], en[0]) in self.connections:
                eng2_position = 0
                eng1_position = 1


            Struct.append(vm.Circle2D(vm.Point2D(position1), self.DF[ne][eng1_position]/2.))
            TG[en[1]] = self.meshes[en[1]].Contour(5)
            Struct.append(vm.Circle2D(vm.Point2D(position2), self.DF[ne][eng2_position]/2.))
            #Definition de la position angulaire initiale
            list_rot = self.initial_position(ne, en)
            if position2[0] == position1[0]:
                if position2[1]-position1[1] > 0:
                    angle = math.pi/2.
                else:
                    angle = -math.pi/2.
            else:
                angle = -math.atan((position2[1]-position1[1])/(position2[0]-position1[0]))
            if num == 0:
                Rot[ne][en[0]] = list_rot[0]-angle
                Rot[ne][en[1]] = list_rot[1]-angle
            else:
                for k1, v1 in Rot.items():
                    if en[0] in v1.keys():
                        Rot[ne][en[0]] = v1[en[0]]
                        delta_rot = Rot[ne][en[0]]-(list_rot[0]-angle)
                Rot[ne][en[1]] = list_rot[1]-angle-delta_rot*((self.meshes[en[0]].z)/(self.meshes[en[1]].z))
            sol = self.gear_rotate([TG[en[0]], TG[en[1]]], [position1, position2], list_rot=[Rot[ne][en[0]], Rot[ne][en[1]]])
            if num == 0:
                L1.extend(sol[0])
            L1.extend(sol[1])
        L1.extend(Struct)
#        G1=vm.Contour2D(L1)
#        G1.MPLPlot()
        return L1

    def pos_axis(self, position):
        """
        Definition of the initial center for all gear (when not given by the user)

        :param position: dictionary define some center position {node2 : [0,0], node4 : [0.12,0] ..}
        :returns: A dictionary where each value coresponds to the center of one gear
        """
        connections = []
        for connection in self.connections:
            connections.append([connection])
        gear_graph = nx.Graph()
        gear_graph.add_nodes_from(self.list_gear)
        # for num_cd, list_connections in enumerate(connections):

        #     (eng1_m, eng2_m) = list_connections[0]
        #     if len(list_connections) > 1:
        #         for (eng1, eng2) in list_connections[1:]:
        #             gear_graph.add_edges_from([(eng1_m, eng1),(eng2_m, eng2)])
        #             eng1_m = eng1
        #             eng2_m = eng2
#        list_line=list(nx.connected_component_subgraphs(gear_graph))
        list_line = [gear_graph.subgraph(c).copy() for c in nx.connected_components(gear_graph)]
        dict_line = {}
        for num_line, list_num_eng in enumerate(list_line):

            for num_eng in list_num_eng:
                dict_line[num_eng] = num_line

        def fun(x):
            obj = 0
            for num_cd, list_connections in enumerate(connections):
                eng1 = dict_line[list_connections[0][0]]
                eng2 = dict_line[list_connections[0][1]]
                obj += (((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5-abs(self.center_distance[num_cd]))**2
                # print((((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5-abs(self.center_distance[num_cd]))**2)
            return obj
        def eg(x):
            ine = []
            for k, val in position.items():
                key = dict_line[k]
                ine.append(x[2*int(key)]-val[0])
                ine.append(x[2*int(key)+1]-val[1])
            return ine
        def ineg(x):
            ine = []
            for num_cd, list_connections in enumerate(connections):
                eng1 = dict_line[list_connections[0][0]]
                eng2 = dict_line[list_connections[0][1]]
                ine.append(((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5-0.999*abs(self.center_distance[num_cd]))
                ine.append(1.001*abs(self.center_distance[num_cd])-((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5)
                # print(ine)
            return ine
        cons = ({'type': 'eq', 'fun' : eg}, {'type': 'ineq', 'fun' : ineg})
        drap = 1
        while drap == 1:
            x0 = tuple(npy.random.random(2*len(list_line))*1)
            Bound = [[0, 10]]*(len(list_line)*2)
            res = minimize(fun, x0, method='SLSQP', bounds=Bound, constraints=cons)


            if (min(ineg(res.x)) > 0) and (max(eg(res.x)) < 1e-7):
                drap = 0
        x_opt = res.x
        centers = {}
        for num_pos, num_eng in enumerate(self.list_gear):
            opt_pos = dict_line[num_eng]
            centers[num_eng] = [x_opt[2*opt_pos], x_opt[2*opt_pos+1]]

        return centers

    def to_markdown(self):

        if hasattr(self, 'infos'):
            return self.infos


class MeshAssembly(DessiaObject):

    """
    Gear Mesh Assembly definition

    :param connections: List of list of tuples defining gear mesh connections.
     Each connection which are in the same list have the same center distance.
     For example [[(1,2),(5,6)], [(2,3)],[(3,4)],[(6,7)],[(7,8)]...]: the connection (1,2) and (5,6) have the same center distance.
    :type connections: List[List[Tuple[int, int]]]
    :param mesh_combinaitons: List of class MechCombination objetcs defining each mesh combination
    :type mesh_combinaitons: List[meshes.MeshCombination]
    :param num_gear_match: List of tuple containing three integer values, each corresponding to one index value. The first corresponds to the gear index in the mesh assembly,
        the second to the gear index in the mesh combination and finally the third to the index of mesh combination.
    :type num_gear_match: List[Tuple[int, int, int]]
    :param safety_factor: Safety factor used for the ISO design
    :type safety_factor: float

    }

    """

    _standalone_in_db = True
    _eq_is_data_eq = True
    _non_serializable_attributes = ['cycle', 'internal_torque',
                                    'general_data', 'dico_gear_match',
                                    'dico_gear_match_inverse']
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']

    def __init__(self, connections: List[List[Tuple[int, int]]],
                 mesh_combinations: List[MeshCombination], num_gear_match: List[Tuple[int, int, int]],
                 strong_links=None, safety_factor: float = 1, name: str = ''):

        self.connections = connections
        self.mesh_combinations = mesh_combinations
        self.num_gear_match = num_gear_match
        self.internal_torque = {}
        self.cycle = {}

        self.dico_gear_match = {}
        self.dico_gear_match_inverse = {}
        for gear in num_gear_match:
            self.dico_gear_match[gear[0]] = (gear[1], gear[2])
            self.dico_gear_match_inverse[(gear[1], gear[2])] = gear[0]
        for num_mesh_combination, mesh_combination in enumerate(mesh_combinations):
            for element in mesh_combination.internal_torque.keys():
                element_mesh_assembly = (self.dico_gear_match_inverse[(element[0], num_mesh_combination)],
                                         self.dico_gear_match_inverse[(element[1], num_mesh_combination)])
                self.internal_torque[element_mesh_assembly] = mesh_combination.internal_torque[element]
            for element in mesh_combination.cycle.keys():
                num_gear_mesh_assembly = self.dico_gear_match_inverse[(element, num_mesh_combination)]
                self.cycle[num_gear_mesh_assembly] = mesh_combination.cycle[element]

        self.strong_links = strong_links
        self.safety_factor = safety_factor


        dict_num_gear_match = {}
        for (n_g_m_a, n_g_m_c, n_m_c) in self.num_gear_match:
            dict_num_gear_match[n_g_m_a] = (n_g_m_c, n_m_c)
        self.dict_num_gear_match = dict_num_gear_match



        self.center_distance = []



        for num_cd, list_connection in enumerate(self.connections):
            for num_mesh_iter, gs in enumerate(list_connection):
                valid = False
                gs_assignate_gear = (self.dico_gear_match[gs[0]][0], self.dico_gear_match[gs[1]][0])
                for mesh_combination in mesh_combinations:
                    for num_mesh_local, gs_local in enumerate(mesh_combination.connections):


                        if set(gs_assignate_gear) == set(gs_local):
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

        dict_gear = {}
        num_gear = 0
        for k, mesh_combination in enumerate(self.mesh_combinations):
            for i, mesh in enumerate(mesh_combination.meshes):
                for match in num_gear_match:
                    if match[1] == i and match[2] == k:
                        num_gear = match[0]
                        break

                dict_gear[num_gear] = mesh
                num_gear += 1
        self.dict_gear = dict_gear


        self.gear_list = [gear for gear in self.dict_gear.values()]
        self.list_gear_index = [gear for gear in self.dict_gear.keys()]



        # self.gear_list = [gear for gear in self.dict_gear.values()]



        coefficient_profile_shift = {}
        for num_mesh, mesh in dict_gear.items():
            coefficient_profile_shift[num_mesh] = mesh.coefficient_profile_shift
        Z = {}
        for num_mesh, mesh in dict_gear.items():
            Z[num_mesh] = mesh.z
        material = {}
        for num_mesh, mesh in dict_gear.items():
            material[num_mesh] = mesh.material
        transverse_pressure_angle_rack = {}
        for num_mesh, mesh in dict_gear.items():
            transverse_pressure_angle_rack[num_mesh] = mesh.rack.transverse_pressure_angle_0
        coeff_gear_addendum = {}
        for num_mesh, mesh in dict_gear.items():
            coeff_gear_addendum[num_mesh] = mesh.rack.coeff_gear_addendum
        coeff_gear_dedendum = {}
        for num_mesh, mesh in dict_gear.items():
            coeff_gear_dedendum[num_mesh] = mesh.rack.coeff_gear_dedendum
        coeff_root_radius = {}
        for num_mesh, mesh in dict_gear.items():
            coeff_root_radius[num_mesh] = mesh.rack.coeff_root_radius
        coeff_circular_tooth_thickness = {}
        for num_mesh, mesh in dict_gear.items():
            coeff_circular_tooth_thickness[num_mesh] = mesh.rack.coeff_circular_tooth_thickness

        self.general_data = []
        for num_graph, list_sub_graph in enumerate(self.sub_graph_dfs):
            num_mesh = 0
            general_data = {'Z': {}, 'connections': [],
                            'material':{}, 'internal_torque':{}, 'cycle':{},
                            'safety_factor':safety_factor}
            input_data = {'center_distance':[], 'transverse_pressure_angle_ini':0,
                          'coefficient_profile_shift':{}, 'transverse_pressure_angle_rack':{},
                          'coeff_gear_addendum':{}, 'coeff_gear_dedendum':{},
                          'coeff_root_radius':{}, 'coeff_circular_tooth_thickness':{}}
            li_connection = []
            for num_cd, list_connection in enumerate(connections):
                for num_mesh_iter, gs in enumerate(list_connection):
                    if (gs in list_sub_graph) or (gs[::-1] in list_sub_graph):
                        li_connection.append((self.dico_gear_match[gs[0]][0], self.dico_gear_match[gs[1]][0]))
                        for num_gear in gs:
                            num_gear_assignation = self.dico_gear_match[num_gear][0]
                            if num_gear in coefficient_profile_shift.keys():
                                input_data['coefficient_profile_shift'][num_gear_assignation] = coefficient_profile_shift[num_gear]
                            if num_gear in transverse_pressure_angle_rack.keys():
                                input_data['transverse_pressure_angle_rack'][num_gear_assignation] = transverse_pressure_angle_rack[num_gear]
                            if num_gear in coeff_gear_addendum.keys():
                                input_data['coeff_gear_addendum'][num_gear_assignation] = coeff_gear_addendum[num_gear]
                            if num_gear in coeff_gear_dedendum.keys():
                                input_data['coeff_gear_dedendum'][num_gear_assignation] = coeff_gear_dedendum[num_gear]
                            if num_gear in coeff_root_radius.keys():
                                input_data['coeff_root_radius'][num_gear_assignation] = coeff_root_radius[num_gear]
                            if num_gear in coeff_circular_tooth_thickness.keys():
                                input_data['coeff_circular_tooth_thickness'][num_gear_assignation] = coeff_circular_tooth_thickness[num_gear]
                            if num_gear in Z.keys():
                                general_data['Z'][num_gear_assignation] = Z[num_gear]
                            if num_gear in material.keys():
                                general_data['material'][num_gear_assignation] = material[num_gear]
                        if num_mesh == 0:
                            input_data['transverse_pressure_angle_ini'] = transverse_pressure_angle[num_mesh]
                    num_mesh += 1

                input_data['center_distance'].append(self.center_distance[num_cd])

            general_data['connections'] = li_connection

            for (eng1, eng2) in list_sub_graph:
                if (eng1, eng2) in self.internal_torque.keys():
                    general_data['internal_torque'][(eng1, eng2)] = self.internal_torque[(eng1, eng2)]
                if (eng2, eng1) in self.internal_torque.keys():
                    general_data['internal_torque'][(eng2, eng1)] = self.internal_torque[(eng2, eng1)]
                if eng1 not in general_data['cycle'].keys():
                    general_data['cycle'][eng1] = self.cycle[eng1]
                if eng2 not in general_data['cycle'].keys():
                    general_data['cycle'][eng2] = self.cycle[eng2]
            self.general_data.append(general_data)

        DessiaObject.__init__(self, name=name)

    def check(self):
        valid = True
        for mesh_combination in self.mesh_combinations:
            valid = (valid and mesh_combination.check())
        return valid

    @classmethod
    def create(cls, center_distance, connections, transverse_pressure_angle,
               coefficient_profile_shift, transverse_pressure_angle_rack,
               coeff_gear_addendum, coeff_gear_dedendum, coeff_root_radius,
               coeff_circular_tooth_thickness, Z, helix_angle, total_contact_ratio_min,
               transverse_contact_ratio_min, percentage_width_difference_pinion_gear,
               max_width_difference_pinion_gear,
               strong_links=None, material=None,
               internal_torque=None, external_torque=None, cycle=None,
               safety_factor=1):

        mesh_combinations = []
        output_data = []

        graph_dfs, _ = gear_graph_simple(connections)
        num_mesh = 0

        num_gear_match = []
        for num_graph, list_sub_graph in enumerate(graph_dfs):


            general_data = {'Z': {}, 'connections': [],
                            'material': {}, 'external_torque': {}, 'cycle': {},
                            'safety_factor': safety_factor}
            input_data = {'center_distance': [], 'transverse_pressure_angle_ini': 0,
                          'coefficient_profile_shift': {}, 'transverse_pressure_angle_rack': {},
                          'coeff_gear_addendum': {}, 'coeff_gear_dedendum': {},
                          'coeff_root_radius': {}, 'coeff_circular_tooth_thickness': {},
                          'helix_angle':{}, 'total_contact_ratio_min':{}, 'transverse_contact_ratio_min':{},
                          'percentage_width_difference_pinion_gear':{}, 'max_width_difference_pinion_gear':{}}


            li_connection = []
            num_mesh_assignation = 0
            num_gear_mesh = 0
            num_gear_assignation = {}
            num_gear_assignation_inverse = {}
            for num_cd, list_connection in enumerate(connections):
                for num_mesh_iter, gs in enumerate(list_connection):
                    if (gs in list_sub_graph) or (gs[::-1] in list_sub_graph):


                        gs_mesh = [0]*len(gs)
                        for i, num_gear in enumerate(gs):

                            if not num_gear in num_gear_assignation.keys():

                                num_gear_assignation[num_gear] = num_gear_mesh

                                num_gear_match.append((num_gear, num_gear_mesh, num_graph))
                                num_gear_mesh += 1
                            gs_mesh[i] = num_gear_assignation[num_gear]

                            if num_gear in coefficient_profile_shift.keys():
                                input_data['coefficient_profile_shift'][num_gear_assignation[num_gear]] = coefficient_profile_shift[num_gear]
                            if num_gear in transverse_pressure_angle_rack.keys():
                                input_data['transverse_pressure_angle_rack'][num_gear_assignation[num_gear]] = transverse_pressure_angle_rack[num_gear]
                            if num_gear in coeff_gear_addendum.keys():
                                input_data['coeff_gear_addendum'][num_gear_assignation[num_gear]] = coeff_gear_addendum[num_gear]
                            if num_gear in coeff_gear_dedendum.keys():
                                input_data['coeff_gear_dedendum'][num_gear_assignation[num_gear]] = coeff_gear_dedendum[num_gear]
                            if num_gear in coeff_root_radius.keys():
                                input_data['coeff_root_radius'][num_gear_assignation[num_gear]] = coeff_root_radius[num_gear]
                            if num_gear in coeff_circular_tooth_thickness.keys():
                                input_data['coeff_circular_tooth_thickness'][num_gear_assignation[num_gear]] = coeff_circular_tooth_thickness[num_gear]
                            if num_gear in helix_angle.keys():
                                input_data['helix_angle'][num_gear_assignation[num_gear]] = helix_angle[num_gear]
                            if num_gear in Z.keys():
                                general_data['Z'][num_gear_assignation[num_gear]] = Z[num_gear]
                            if num_gear in material.keys():
                                general_data['material'][num_gear_assignation[num_gear]] = material[num_gear]
                        li_connection.append((gs_mesh[0], gs_mesh[1]))
                        if num_mesh_assignation == 0:

                            input_data['transverse_pressure_angle_ini'] = transverse_pressure_angle[num_mesh]
                        input_data['total_contact_ratio_min'][li_connection[-1]] = total_contact_ratio_min[gs]
                        input_data['transverse_contact_ratio_min'][li_connection[-1]] = transverse_contact_ratio_min[gs]
                        input_data['percentage_width_difference_pinion_gear'][li_connection[-1]] = percentage_width_difference_pinion_gear[gs]
                        input_data['max_width_difference_pinion_gear'][li_connection[-1]] = max_width_difference_pinion_gear[gs]
                        num_mesh += 1
                        num_mesh_assignation += 1
                input_data['center_distance'].append(center_distance[num_cd])
            general_data['connections'] = li_connection
            for (eng1, eng2) in list_sub_graph:
                # if (eng1,eng2) in internal_torque.keys():
                #     general_data['external_torque'][(eng1,eng2)]=internal_torque[(eng1,eng2)]
                # if (eng2,eng1) in internal_torque.keys():
                #     general_data['external_torque'][(eng2,eng1)]=internal_torque[(eng2,eng1)]

                if not eng1 in general_data['external_torque'].keys():
                    if eng1 in external_torque.keys():
                        general_data['external_torque'][num_gear_assignation[eng1]] = external_torque[eng1]
                if not eng2 in general_data['external_torque'].keys():
                    if eng2 in external_torque.keys():
                        general_data['external_torque'][num_gear_assignation[eng2]] = external_torque[eng2]
                if eng1 not in general_data['cycle'].keys():

                    general_data['cycle'][num_gear_assignation[eng1]] = cycle[eng1]
                if eng2 not in general_data['cycle'].keys():
                    general_data['cycle'][num_gear_assignation[eng2]] = cycle[eng2]


            output_data.append(general_data)
            xt = dict(list(input_data.items()) + list(general_data.items()))

            mesh_combinations.append(MeshCombination.create(**xt))

        mesh_assembly = cls(connections, mesh_combinations, num_gear_match,
                            strong_links, safety_factor)
        return mesh_assembly

    def _get_graph_dfs(self):

        """
        :returns:
        """
        _graph_dfs, _ = gear_graph_simple(self.connections)

        return _graph_dfs
    sub_graph_dfs = property(_get_graph_dfs)

    def _get_list_gear(self):

        """
        :returns: List with the gear indexes of the mesh assembly
        """
        _, _list_gear = gear_graph_simple(self.connections)

        return _list_gear
    list_gear = property(_get_list_gear)

    def SVGExport(self, name, position):
        """ Export SVG graph of all gear mesh combinations

        :param name: name of the svg file
        :param position: dictionary define some center position {node2 : [0,0], node4 : [0.12,0] ..}

        :results: SVG graph

        in the position dictionary, you have to be coherent with the center position

            * for exemple, if the center-distance of the mesh1 (node1, node2) is 0.117 m you can define position such as:

                * {node1 : [0,0], node2 : [0.117,0]}
                * {node1 : [0,0]}
        """

        centers = self.pos_axis(position)
        L = []
        for mesh_assembly_iter in self.mesh_combinations:
            position_svg = {}
            for num_gear, pos in centers.items():
                if num_gear in mesh_assembly_iter.Z.keys():
                    position_svg[num_gear] = pos
            L.extend(mesh_assembly_iter.SVGExport('gear', position_svg))
        G1 = vm.Contour2D(L)
        G1.MPLPlot()

    def FreeCADExport(self, fcstd_filepath, centers={}, axis=(1, 0, 0), export_types=['fcstd'], python_path='python',
                      path_lib_freecad='/usr/lib/freecad/lib'):
        """ Export 3D volume to FreeCAD

        :param file_path: file path for the freecad file
        :param center: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
        :param axis: direction of gear mesh rotation

        :results: export of a FreeCAD file
        """
        for ma in self.mesh_combinations:
            ma.FreeCADExport(fcstd_filepath, centers, axis, python_path, path_lib_freecad, export_types)

    def update(self, optimizer_data):

        output_x = []
        for num_graph, list_sub_graph in enumerate(self.sub_graph_dfs):
            num_mesh = 0
            input_data = {'center_distance':[], 'transverse_pressure_angle_ini':[],
                          'coefficient_profile_shift':{}, 'transverse_pressure_angle_rack':{},
                          'coeff_gear_addendum':{}, 'coeff_gear_dedendum':{},
                          'coeff_root_radius':{}, 'coeff_circular_tooth_thickness':{}, 'total_contact_ratio_min':{},
                          'transverse_contact_ratio_min':{}, 'percentage_width_difference_pinion_gear':{},
                          'max_width_difference_pinion_gear':{}}
            li_connection = []
            for num_cd, list_connection in enumerate(self.connections):
                for num_mesh_iter, (eng1, eng2) in enumerate(list_connection):
                    if ((eng1, eng2) in list_sub_graph) or ((eng2, eng1) in list_sub_graph):
                        li_connection.append((eng1, eng2))
                        num_gear_assignation_1 = self.dico_gear_match[eng1][0]
                        num_gear_assignation_2 = self.dico_gear_match[eng2][0]
                        for key, list_value in optimizer_data.items():

                            if key in ['coefficient_profile_shift',
                                       'transverse_pressure_angle_rack',
                                       'coeff_gear_addendum', 'coeff_gear_dedendum',
                                       'coeff_root_radius', 'coeff_circular_tooth_thickness']:

                                input_data[key][num_gear_assignation_1] = optimizer_data[key][eng1]
                                input_data[key][num_gear_assignation_2] = optimizer_data[key][eng2]
                            elif key in ['center_distance']:
                                input_data[key].append(optimizer_data[key][num_cd])
                            elif key in ['transverse_pressure_angle']:
                                input_data['transverse_pressure_angle_ini'].append(optimizer_data[key][num_mesh])
                            input_data['total_contact_ratio_min'][(num_gear_assignation_1, num_gear_assignation_2)] \
                            = optimizer_data['total_contact_ratio_min'][(eng1, eng2)]
                            input_data['transverse_contact_ratio_min'][(num_gear_assignation_1, num_gear_assignation_2)] \
                            = optimizer_data['transverse_contact_ratio_min'][(eng1, eng2)]
                            input_data['percentage_width_difference_pinion_gear'][(num_gear_assignation_1, num_gear_assignation_2)] \
                            = optimizer_data['percentage_width_difference_pinion_gear'][(eng1, eng2)]
                            input_data['max_width_difference_pinion_gear'][(num_gear_assignation_1, num_gear_assignation_2)] \
                            = optimizer_data['max_width_difference_pinion_gear'][(eng1, eng2)]
                    num_mesh += 1
            input_data['transverse_pressure_angle_ini'] = input_data['transverse_pressure_angle_ini'][0]
            xt = dict(list(input_data.items())+list(self.general_data[num_graph].items()))
            output_x.append(xt)

#            if self.save!=optimizer_data:
            self.mesh_combinations[num_graph].update(**xt)
        return output_x

    def update_helix_angle(self, optimizer_data):
        output_x = []
        for num_graph, list_sub_graph in enumerate(self.sub_graph_dfs):
            num_mesh = 0
            input_data = {'helix_angle':{}, 'total_contact_ratio_min':{}, 'transverse_contact_ratio_min':{},
                          'percentage_width_difference_pinion_gear':{},
                          'max_width_difference_pinion_gear':{}}
            li_connection = []
            for num_cd, list_connection in enumerate(self.connections):
                for num_mesh_iter, (eng1, eng2) in enumerate(list_connection):
                    if ((eng1, eng2) in list_sub_graph) or ((eng2, eng1) in list_sub_graph):
                        li_connection.append((eng1, eng2))
                        num_gear_assignation_1 = self.dico_gear_match[eng1][0]
                        num_gear_assignation_2 = self.dico_gear_match[eng2][0]
                        for key, list_value in optimizer_data.items():
                            if key in ['helix_angle']:
                                if optimizer_data[key]:
                                    input_data[key][num_gear_assignation_1] = optimizer_data[key][eng1]
                                    input_data[key][num_gear_assignation_2] = optimizer_data[key][eng2]
                            input_data['total_contact_ratio_min'][(num_gear_assignation_1, num_gear_assignation_2)] = \
                            optimizer_data['total_contact_ratio_min'][(eng1, eng2)]
                            input_data['transverse_contact_ratio_min'][(num_gear_assignation_1, num_gear_assignation_2)] = \
                            optimizer_data['transverse_contact_ratio_min'][(eng1, eng2)]
                            input_data['percentage_width_difference_pinion_gear'][(num_gear_assignation_1, num_gear_assignation_2)] = \
                            optimizer_data['percentage_width_difference_pinion_gear'][(eng1, eng2)]
                            input_data['max_width_difference_pinion_gear'][(num_gear_assignation_1, num_gear_assignation_2)] = \
                            optimizer_data['max_width_difference_pinion_gear'][(eng1, eng2)]
                    num_mesh += 1

            xt = dict(list(input_data.items()))
            output_x.append(xt)

#            if self.save!=optimizer_data:
            self.mesh_combinations[num_graph].update_helix_angle(**xt)
        return output_x
    
    def list_gear_per_mesh_combinations(self):
        list_gear_per_mesh_combinations = []
        for i,mesh_combination in enumerate(self.mesh_combinations):
            gear_per_mesh_combinations=[]
            for gear_match in self.num_gear_match:
                if gear_match[2] == i:
                    gear_per_mesh_combinations.append(gear_match[0])
            list_gear_per_mesh_combinations.append(gear_per_mesh_combinations)
        
        return list_gear_per_mesh_combinations
                    
    def plot_data(self):
        """
        2D mesh combination visualization

        :returns: List of Primitives groups for the data_plot
        """
        list_colors = [vmp.colors.BLACK, vmp.colors.RED, vmp.colors.BLUE]
        export_data = []
        x_position = 0
        count_previous_mesh = 0
        centers_yz = self.pos_axis({self.list_gear[0]:[0, 0]})
        for i, mesh_combination in enumerate(self.mesh_combinations):
            centers = {}
            for j, mesh in enumerate(mesh_combination.meshes):
                for gear_index_mesh_assembly, (gear_index_mesh_comb, index_mesh_comb) in zip(self.dico_gear_match.keys(), self.dico_gear_match.values()):
                    if i == index_mesh_comb and j == gear_index_mesh_comb:
                        gear_index = gear_index_mesh_assembly
                        center_yz = centers_yz[gear_index]

                        if j == 0:
                            if i == 0:
                                centers[j] = (x_position, center_yz[0], center_yz[1])
                            else:
                                centers[j] = (x_position, center_yz[0] + 3*max([mesh.outside_diameter for mesh in mesh_combination.meshes])/2, center_yz[1])
                                # centers[j] = (x_position, previous_center_yz[0], previous_center_yz[0])

                        else:
                            if i == count_previous_mesh:
                                centers[j] = (x_position, center_yz[0], center_yz[1])
                            else:
                                centers[j] = (x_position, center_yz[0] + 3*max([mesh.outside_diameter for mesh in mesh_combination.meshes])/2, center_yz[1])


            primitives = mesh_combination.plot_data(centers)[0].primitives
            for primitive in primitives:
                if type(primitive) is not vmp.Text and type(primitive) is not vmp.Circle2D:
                    primitive.edge_style.color_stroke = list_colors[i]

            export_data.extend(primitives)

            count_previous_mesh += i
        return [vmp.PrimitiveGroup(primitives=export_data)]

    def volmdlr_primitives(self):
        """
        Generation of the 3D volume for all the gear mesh

        :results: list of 3D volmdlr component
        """

        primitives = []


        offset = 0.01
        x_position = 0
        centers_yz = self.pos_axis({self.list_gear[0]:[0, 0]})
        # count = 0


        for i, mesh_combination in enumerate(self.mesh_combinations):
            centers = {}
            for j, mesh in enumerate(mesh_combination.meshes):
                for gear_index_mesh_assembly, (gear_index_mesh_comb, index_mesh_comb) in zip(self.dico_gear_match.keys(), self.dico_gear_match.values()):
                    if i == index_mesh_comb and j == gear_index_mesh_comb:
                        gear_index = gear_index_mesh_assembly
                        center_yz = centers_yz[gear_index]

                        if j == 0:
                            if i == 0:
                                centers[j] = (x_position, center_yz[0], center_yz[1])
                            else:
                                x_position += max(mesh_combination.gear_width)/2
                                centers[j] = (x_position, center_yz[0], center_yz[1])
                                # centers[j] = (x_position, previous_center_yz[0], previous_center_yz[0])
                        else:
                            centers[j] = (x_position, center_yz[0], center_yz[1])
            primitives.extend(mesh_combination.volmdlr_primitives(centers=centers))
            # previous_center_yz = center_yz
            x_position += max(mesh_combination.gear_width)/2 + offset
        return primitives



    def pos_axis(self, position):
        """
        Definition of the initial center for all gear (when not given by the user)

        :param position: dictionary define some center position {node2 : [0,0], node4 : [0.12,0] ..}
        :returns: A dictionary where each value coresponds to the center of one gear
        """
        gear_graph = nx.Graph()

        gear_graph.add_nodes_from(self.list_gear)

        for num_cd, list_connections in enumerate(self.connections):
            (eng1_m, eng2_m) = list_connections[0]
            if len(list_connections) > 1:
                for (eng1, eng2) in list_connections[1:]:
                    gear_graph.add_edges_from([(eng1_m, eng1), (eng2_m, eng2)])
                    eng1_m = eng1
                    eng2_m = eng2
#        list_line=list(nx.connected_component_subgraphs(gear_graph))
        list_line = [gear_graph.subgraph(c).copy() for c in nx.connected_components(gear_graph)]
        dict_line = {}

        for num_line, list_num_eng in enumerate(list_line):
            for num_eng in list_num_eng:
                dict_line[num_eng] = num_line

        def fun(x):
            obj = 0
            for num_cd, list_connections in enumerate(self.connections):
                eng1 = dict_line[list_connections[0][0]]
                eng2 = dict_line[list_connections[0][1]]
                obj += (((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5-self.center_distance[num_cd])**2
            return obj
        def eg(x):
            ine = []
            for k, val in position.items():
                key = dict_line[k]
                ine.append(x[2*int(key)]-val[0])
                ine.append(x[2*int(key)+1]-val[1])
            return ine
        def ineg(x):
            ine = []
            for num_cd, list_connections in enumerate(self.connections):
                eng1 = dict_line[list_connections[0][0]]
                eng2 = dict_line[list_connections[0][1]]

                ine.append(((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5-0.999*abs(self.center_distance[num_cd]))
                ine.append(1.001*abs(self.center_distance[num_cd])-((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5)
            return ine
        cons = ({'type': 'eq', 'fun' : eg}, {'type': 'ineq', 'fun' : ineg})
        drap = 1
        while drap == 1:
            x0 = tuple(npy.random.random(2*len(list_line))*1)
            Bound = [[0, 1]]*(len(list_line)*2)
            res = minimize(fun, x0, method='SLSQP', bounds=Bound, constraints=cons)

            if (min(ineg(res.x)) > 0) and (max(eg(res.x)) < 1e-7):
                drap = 0
        x_opt = res.x
        centers = {}
        for num_pos, num_eng in enumerate(self.list_gear):
            opt_pos = dict_line[num_eng]
            centers[num_eng] = [x_opt[2*opt_pos], x_opt[2*opt_pos+1]]
        return centers

def gear_graph_simple(connections):
    """
    NetworkX graph construction

    :param connections : List of tuples defining gear mesh connection [[(node1,node2)], [(node2,node3)]...]

    :returns:
        *
        *List with the gear indexes of the mesh assembly
    """
    # NetworkX graph construction
    list_gear = [] # list of all gears
    compt_mesh = 0 # number of gear mesh
    for gs in connections:
        for (eng1, eng2) in gs:
            compt_mesh += 1
            if eng1 not in list_gear:
                list_gear.append(eng1)
            if eng2 not in list_gear:
                list_gear.append(eng2)
    # Construction of one graph including all different connection types (gear_mesh, same_speed, same_shaft)
    gear_graph = nx.Graph()
    gear_graph.add_nodes_from(list_gear)
    for list_edge in connections:
        gear_graph.add_edges_from(list_edge)
#    sub_graph=list(nx.connected_component_subgraphs(gear_graph))
    sub_graph = [gear_graph.subgraph(c).copy() for c in nx.connected_components(gear_graph)]

    sub_graph_dfs = []
    for s_graph in sub_graph:
        node_init = list(s_graph.nodes())[0]
        sub_graph_dfs.append(list(nx.dfs_edges(s_graph, node_init)))
    return sub_graph_dfs, list_gear

def gear_graph_complex(connections, strong_link):
    # Construction of one graph include all different connection type (gear_mesh, same_speed, same_shaft)
    _, list_gear = gear_graph_simple(connections)
    gear_graph = nx.Graph()
    gear_graph.add_nodes_from(list_gear)
    for list_edge in connections:
        gear_graph.add_edges_from(list_edge, typ='gear_mesh')
        li_shaft1 = []
        li_shaft2 = []
        for eng1, eng2 in list_edge:
            li_shaft1.append(eng1)
            li_shaft2.append(eng2)
        if len(li_shaft1) > 1:
            for pos_gear, num_gear in enumerate(li_shaft1[1:]):
                valid_strong_ling = False
                for list_strong_link in strong_link:
                    if (num_gear in list_strong_link) and (li_shaft1[pos_gear] in list_strong_link):
                        valid_strong_ling = True
                if valid_strong_ling:
                    gear_graph.add_edges_from([(num_gear, li_shaft1[pos_gear])], typ='same_speed')
                else:
                    gear_graph.add_edges_from([(num_gear, li_shaft1[pos_gear])], typ='same_shaft')
        if len(li_shaft2) > 1:
            for pos_gear, num_gear in enumerate(li_shaft2[1:]):
                valid_strong_ling = False
                for list_strong_link in strong_link:
                    if (num_gear in list_strong_link) and (li_shaft2[pos_gear] in list_strong_link):
                        valid_strong_ling = True
                if valid_strong_ling:
                    gear_graph.add_edges_from([(num_gear, li_shaft2[pos_gear])], typ='same_speed')
                else:
                    gear_graph.add_edges_from([(num_gear, li_shaft2[pos_gear])], typ='same_shaft')
    connections_dfs = list(nx.dfs_edges(gear_graph, list_gear[0]))
    # construction of a graph without same_shaft attribute
    gear_graph_kinematic = copy.deepcopy(gear_graph)
    for edge, typ in nx.get_edge_attributes(gear_graph_kinematic, 'typ').items():
        if typ == 'same_shaft':
            gear_graph_kinematic.remove_edges_from([edge])
    connections_kinematic_dfs = list(nx.dfs_edges(gear_graph_kinematic, list_gear[0]))
    return connections_dfs, connections_kinematic_dfs, gear_graph

class ValidGearDiameterError(Exception):
    def __init__(self):
        super().__init__('Fail base diameter is greater than pitch diameter')