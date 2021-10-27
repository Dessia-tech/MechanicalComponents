#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Bearings module

"""

__title__ = 'Bearings'
__description__ = 'Ball and roller bearings'

import numpy as npy
npy.seterr(divide='raise', over='ignore', under='ignore')
#import math as mt
from scipy import interpolate
#import os
import volmdlr as vm
import plot_data
import volmdlr.primitives3d as primitives3d
import volmdlr.primitives2d as primitives2d
import math
from dessia_common import DessiaObject, dict_merge, Evolution
from mechanical_components import shafts_assembly

#import copy
import json

from scipy.optimize import fsolve
import networkx as nx
import matplotlib.pyplot as plt
from typing import TypeVar, List
from dataclasses import dataclass
import inspect

#import pandas
import pkg_resources

#from mechanical_components.bearings_snr import RadialRollerBearingSNR
#from mechanical_components.catalogs.ISO_bearings \
#    import iso_bearings, iso_rollers, iso_radial_clearances, bearing_rules

import genmechanics
import genmechanics.linkages as linkages
import genmechanics.loads as gm_loads
import genmechanics.unidimensional as unidimensional

from mechanical_components.tools import StringifyDictKeys

import matplotlib.colors


class Mounting(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, left: bool = False, right: bool = False,
                 name: str = ''):

        self.left = left
        self.right = right
        DessiaObject.__init__(self, name=name)

    @property
    def both(self):
        if self.left and self.right:
            return True
        else:
            return False

    @property
    def free(self):
        if not self.left and not self.right:
            return True
        else:
            return False

class CombinationMounting(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, mountings: List[Mounting],
                 name: str = ''):

        self.mountings = mountings
        DessiaObject.__init__(self, name=name)

class Linkage(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, ball_joint: bool = False, cylindric_joint: bool = False,
                 name: str = ''):

        self.ball_joint = ball_joint
        self.cylindric_joint = cylindric_joint
        DessiaObject.__init__(self, name=name)

class SelectionLinkage(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, linkages: List[Linkage],
                 name: str = ''):

        self.linkages = linkages
        DessiaObject.__init__(self, name=name)

#oil_kinematic_viscosity
#iso_vg_1500={'data':[[47.21238870380181,922.5481847223729],
#                     [76.41592953982855,191.5471560481642],
#                     [110.70796589064605,54.90426918109079]],'x':'Linear','y':'Log'}
#iso_vg_1000={'data':[[41.68141577143845,877.2173704075102],
#                     [62.477877333256444,261.78400435754804],
#                     [100.53097486106454,57.74149074755608]],'x':'Linear','y':'Log'}
#iso_vg_680={'data':[[38.80530942959304,777.3038394524846],
#                    [57.168142067138206,251.441902731903],
#                    [89.46902691125547,61.96159004047747]],'x':'Linear','y':'Log'}
#iso_vg_460={'data':[[36.15044283907511,580.3315115122488],
#                    [59.159291488756054,159.77215436382392],
#                    [85.48672598293739,53.80881018274548]],'x':'Linear','y':'Log'}
#iso_vg_320={'data':[[32.16814191075703,551.8160309554283],
#                    [57.38937973338331,131.93199565920736],
#                    [81.06194763704671,48.65076991453211]],'x':'Linear','y':'Log'}
#iso_vg_220={'data':[[29.95575273781169,407.8526478060212],
#                    [56.725664649565616,96.53454045940936],
#                    [83.05309705866458,35.24081769455843]],'x':'Linear','y':'Log'}
#iso_vg_150={'data':[[27.964601231111455,307.58489032135725],
#                    [50.97345196587479,89.9597273365993],
#                    [87.25663773831015,23.313898800373792]],'x':'Linear','y':'Log'}
#iso_vg_100={'data':[[33.05309674590221,148.89034266049572],
#                    [60.7079655778837,41.82579554586569],
#                    [91.23893866662823,15.115797660575524]],'x':'Linear','y':'Log'}
#iso_vg_68={'data':[[29.95575273781169,113.42390186278217],
#                   [56.94690231581072,34.19139868782362],
#                   [89.91150432882806,11.749567781125915]],'x':'Linear','y':'Log'}
#iso_vg_46={'data':[[29.070795817584123,76.56429373059628],
#                   [60.48672582655621,21.946066333596434],
#                   [99.4247781895095,7.316992362396092]],'x':'Linear','y':'Log'}
#iso_vg_32={'data':[[27.52212381353886,56.0220352459023],
#                   [58.27433665361086,17.058761946001017],
#                   [82.16814222351938,8.510952177980519]],'x':'Linear','y':'Log'}
#iso_vg_22={'data':[[30.619469906711767,32.840621976693456],
#                   [57.38937973338331,12.481883087082235],
#                   [90.79646124905564,4.939173694948054]],'x':'Linear','y':'Log'}
#iso_vg_15={'data':[[23.982300928318086,28.519522512213047],
#                   [44.115044695711276,13.126893028298408],
#                   [77.07964670872863,4.889651598606255]],'x':'Linear','y':'Log'}
#iso_vg_10={'data':[[25.088495514790754,17.231530421142708],
#                   [46.548673619984115,8.092752773048398],
#                   [66.01769875891955,4.649391098015179]],'x':'Linear','y':'Log'}


#Ordre de rangement du coefficient de contamination de l'huile: Dpw mini/ Dpw maxi/ grade/ coeff de contamination
dict_oil_contamination = {0:{0.1:{1:1, 2:0.7, 3:0.55, 4:0.4, 5:0.2, 6:0.05, 7:0}},
                          0.1:{math.inf:{1:1, 2:0.85, 3:0.7, 4:0.5, 5:0.3, 6:0.05, 7:0}}}

class TypeEvolution(DessiaObject):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, _linear: bool = True, _log: bool = False, name: str = ''):
        self._linear = _linear
        self._log = _log

        DessiaObject.__init__(self, name=name)

class Oil(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, oil_kinematic_viscosity: Evolution, oil_temperature: Evolution,
                 oil_kinematic_viscosity_type: TypeEvolution, oil_temperature_type: TypeEvolution,
                 name=''):
        self.oil_kinematic_viscosity = oil_kinematic_viscosity
        self.oil_temperature = oil_temperature
        self.oil_kinematic_viscosity_type = oil_kinematic_viscosity_type
        self.oil_temperature_type = oil_temperature_type

        self.A, self.B = self.kinematic_viscosity(oil_temperature, oil_kinematic_viscosity)

        DessiaObject.__init__(self, name=name)

    def kinematic_viscosity(self, oil_temperature, oil_kinematic_viscosity):
        if self.oil_kinematic_viscosity_type._log and self.oil_temperature_type._linear:
            evol_temp = oil_temperature.evolution
            evol_kc = oil_kinematic_viscosity.evolution
            A = (math.log10(math.log10(0.6+evol_kc[0]))-math.log10(math.log10(0.6+evol_kc[-1])))/(math.log10(evol_temp[0])-math.log10(evol_temp[-1]))
            B = math.log10(math.log10(0.6+evol_kc[0]))-A*math.log10(evol_temp[0])
            return A, B
#        oil_kinematic_viscosity_curve = {}
#        for (temp, kv) in zip(oil_temperature, oil_kinematic_viscosity)
#            val_np = npy.array(val)
#            if key not in ['x','y']:
#                oil_kinematic_viscosity_curve[key] = {}
#                A = (math.log10(math.log10(0.6+val_np[0,1]))-math.log10(math.log10(0.6+val_np[-1,1])))/(math.log10(val_np[0,0])-math.log10(val_np[-1,0]))
#                B = math.log10(math.log10(0.6+val_np[0,1]))-A*math.log10(val_np[0,0])
#                oil_kinematic_viscosity_curve[key]['A'] = A
#                oil_kinematic_viscosity_curve[key]['B'] = B
#        return oil_kinematic_viscosity_curve['data']

    def oil_parameter_contamination(self, Dpw, grade):
        for k, v in dict_oil_contamination.items():
            if (Dpw >= k) and (Dpw < list(v.keys())[0]):
                return list(v.values())[0][grade]


oil_iso_vg_1500 = Oil(oil_kinematic_viscosity=Evolution([47.21238870380181, 76.41592953982855, 110.70796589064605]),
                      oil_temperature=Evolution([922.5481847223729, 191.5471560481642, 54.90426918109079]),
                      oil_kinematic_viscosity_type=TypeEvolution(False, True),
                      oil_temperature_type=TypeEvolution(True, False))

#oil_iso_vg_1000=Oil(iso_vg_1000)
#oil_iso_vg_680=Oil(iso_vg_680)
#oil_iso_vg_460=Oil(iso_vg_460)
#oil_iso_vg_320=Oil(iso_vg_320)
#oil_iso_vg_220=Oil(iso_vg_220)
#oil_iso_vg_150=Oil(iso_vg_150)
#oil_iso_vg_100=Oil(iso_vg_100)
#oil_iso_vg_68=Oil(iso_vg_68)
#oil_iso_vg_46=Oil(iso_vg_46)
#oil_iso_vg_32=Oil(iso_vg_32)
#oil_iso_vg_22=Oil(iso_vg_22)
#oil_iso_vg_15=Oil(iso_vg_15)
#oil_iso_vg_10=Oil(iso_vg_10)

class Material(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, weibull_e: float, weibull_c: float, weibull_h: float,
                 B1: float, mu_delta: float, c_gamma: float, name: str = ''):

        """
        Definition of the object material for ring

        :param weibull_e: weibull parameter e, 10/9 for point contact and 9/8 for linear contact
        :param weibull_c: weibull parameter c
        :param weibull_h: weibull parameter h
        """
        self.weibull_e = weibull_e
        self.weibull_c = weibull_c
        self.weibull_h = weibull_h
        self.B1 = B1
        self.mu_delta = mu_delta
        self.c_gamma = c_gamma

        DessiaObject.__init__(self, name=name)

material_iso = Material(9/8., 31/3., 7/3., 551.13373/0.483, 0.83, 0.05)

class RadialBearing(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _non_hash_attributes = ['E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _generic_eq = True

    symmetric = None
    taking_loads = None
    generate_axial_load = None
    linkage = None

    def __init__(self, d: float, D: float, B: float, alpha: float, i: int, Z: int, Dw: float, Cr: float = None,
                 C0r: float = None, material: Material = material_iso,
                 contact_type_point: bool = True, contact_type_linear: bool = False, contact_type_mixed: bool = False,
                 width: float = 0.02, mass: float = None,
                 oil_speed_limit: float = None,
                 grease_speed_limit: float = None,
                 name: str = '', infos: List[str] = None, speed_infos: List[float] = None):


        self.d = d
        self.D = D
        self.B = B
        self.Dpw = (d + D)/2.
        self.i = i
        self.width = width
        if Dw is None:
            self.Dw = min((D - d)/4, 0.85*B)
        else:
            self.Dw = Dw
        if Z is None:
            self.Z = int(math.pi*self.Dpw/self.Dw)
        else:
            self.Z = Z

        self.alpha = alpha
        self.material = material
        self.contact_type_point = contact_type_point
        self.contact_type_linear = contact_type_linear
        self.contact_type_mixed = contact_type_mixed
        if not infos:
            self.infos = []
        else:
            self.infos = infos
        if not speed_infos:
            self.speed_infos = []
        else:
            self.speed_infos = speed_infos


        self.oil_speed_limit = oil_speed_limit
        self.grease_speed_limit = grease_speed_limit


        if Cr is not None:
            self.Cr = Cr
        if C0r is not None:
            self.C0r = C0r

        # estimation parameter for plot (define exactly on geometry object)
        self.E = self.Dpw + self.Dw
        self.F = self.Dpw - self.Dw
        self.d1 = self.F + 0.6*self.Dw
        self.D1 = self.E - 0.6*self.Dw
        self.radius = 5e-4
        self.slack = (self.E-self.F-2*self.Dw)/4.
        self.name = name
        if mass is None:
            self.mass = self.mass()
        else:
            self.mass = mass
        self.cost = self.mass*self.cost_coefficient + self.cost_constant

        DessiaObject.__init__(self, name=name)

    @property
    def speed_limit(self):
        if self.grease_speed_limit:
            speed_limit = self.grease_speed_limit
        else:
            speed_limit = None

        if self.oil_speed_limit:
            if speed_limit:
                return max(speed_limit, self.oil_speed_limit)
            else:
                return self.oil_speed_limit

        else:
            return speed_limit

    def check(self):
        if self.d <= 0.:
            return False
        if self.d >= self.D:
            return False
        return True

    @classmethod
    def estimate_base_life_time(cls, Fr, N, t, Cr):
        total_cycles = 0.
        Pr = 0.
        for fr, ni, ti in zip(Fr, N, t):
            C = fr**(cls.coeff_baselife)

            if C != 0.:
                cycles = ni * ti * 2 * math.pi
                Pr += cycles * C
                total_cycles += cycles

        if total_cycles == 0.:
            return math.inf #TODO
        else:
            Pr = (Pr / total_cycles) ** (1/cls.coeff_baselife)

            L10 = (Cr/Pr)**(cls.coeff_baselife)
            return L10

    def base_life_time(self, Fr, Fa, N, t, Cr, verbose=False):
        """
        Lifetime in millions of cycles for 90% fiability

        :param Fr: a list of radial forces for each usecase
        :param Fa: a list of axial forces for each usecase
        :param N: a list of rotating speeds in rad/s for each usecase
        :param t: a list of operating times in seconds for each usecase
        :param Cr: a float define the base dynamic load of the bearing
        """
        total_cycles = 0.
        Pr = 0.

      
        for fr, fa, ni, ti in zip(Fr, Fa, N, t):



            C = self.equivalent_dynamic_load(fr, fa)**(self.coeff_baselife)

           
            if ti == 0:
                ti = 1
            if C != 0.:
                cycles = ni * ti * 2 * math.pi

                Pr += cycles * C
                total_cycles += cycles


       
        if total_cycles != 0:

            Pr = (Pr / total_cycles) ** (1/self.coeff_baselife)

            L10 = (Cr/Pr)**(self.coeff_baselife)
            if verbose:
                print(Cr)
                print(Pr)


            return L10
        else:
            raise BearingL10Error()

    def adjusted_life_time(self, Fr, Fa, N, t, T=None, Cr=None, C0r=None, S=0.9):
        """
        Adjusted Lifetime in millions of cycles for a 100* S % fiability

        :param Fr: a list of radial forces for each usecase
        :param Fa: a list of axial forces for each usecase
        :param N: a list of rotating speeds in rad/s for each usecase
        :param t: a list of operating times in seconds for each usecase
        :param T: a list of operating temperature in celcius degree for each usecase
        :param Cr: a float define the base dynamic load of the bearing
        :param C0r: a float define the base static load of the bearing
        :param S: fiability between 0 and 1

        """
        if Cr is not None:
            self.Cr = Cr
        if C0r is not None:
            self.C0r = C0r
        total_cycles = 0.
        nci_Lpi = 0.

        for fr, fa, n, ti   in zip(Fr, Fa, N, t, ):
            if (((fr != 0.) or (fa != 0.)) and (n > 0.)):
                cycles = n * ti * 2 * math.pi
                total_cycles += cycles

                a1 = ((1-self.material.c_gamma)
                      * (math.log(1/S)/math.log(100/90.))**(1/self.material.weibull_e)
                      + self.material.c_gamma)
                L10 = self.base_life_time([fr], [fa], [n], [ti], self.Cr)
                # Pr = self.equivalent_dynamic_load(fr, fa)
                # viscosité cinématique de référence
                # if n < (1000*2*math.pi/60.):
                #     nu1 = 45000*(n*60/(2*math.pi))**(-0.83)*(self.Dpw*1e3)**(-0.5)
                # else:
                #     nu1 = 4500*(n*60/(2*math.pi))**(-0.5)*(self.Dpw*1e3)**(-0.5)

#                coeff_oil = self.oil.oil_kinematic_viscosity_curve
                # nu = 10**(10**(self.oil.A*math.log10(Ti)+self.oil.B))-0.6
                # kappa = nu/nu1
                # # Oil Contamination
                # ec = self.oil.OilParameterContamination(self.Dpw,3)
                # # Wear limit load
                # if self.Dpw<0.1:
                #     Cu = self.C0r/8.2
                # else:
                #     Cu = self.C0r/8.2*(100/(self.Dpw*1e3))**0.3

                # kappa = min(kappa, 4)
                # a_iso = self.a_iso(kappa, ec, Cu, Pr)
                # a_iso = min(50., a_iso)

                Lpi = a1*L10 # Corrected lifetime
                nci_Lpi += cycles/Lpi
        if nci_Lpi > 0:
            return total_cycles / nci_Lpi
        else:
            return math.inf

#    def CheckFNRRules(self, Fr, Fa, N):
#        check_rules = True
#        val_rules = math.inf
#        for fr, fa, n  in zip(Fr, Fa, N):
#            if self.typ_bearing == 'radial_roller_bearing':
#                rules_snr = RadialRollerBearingSNR(self.d, self.D, self.B, self.Z, self.alpha, self.Dpw)
#                check_rules_iter, val_rules_iter = rules_snr.Ruleaxial_load(fr, fa, n, level_axial_load='constant_load')
#                val_rules = min(val_rules_iter, val_rules)
#                if check_rules_iter == False:
#                    check_rules = False
#        return check_rules, val_rules

    def mass(self):
        # TODO: enhance this but without querying CAD volumes!
        return 7800 * math.pi*self.B*(self.D-self.d) * (self.d+self.D)

    def volmdlr_primitives(self, center=vm.O3D, axis=vm.X3D):
        # TODO: mutualization of this in parent class?
        axis.normalize()

        y = axis.random_unit_normal_vector()
        z = axis.cross(y)

        #Internal Ring
        cbi1, bi1 = self.internal_ring_contour()
        if cbi1:
            IRC = vm.wires.Contour2D([cbi1]+bi1.primitives)
        else:
            IRC = vm.wires.Contour2D(bi1.primitives)
        irc = primitives3d.RevolvedProfile(center, axis, z, IRC, center,
                                           axis, angle=2*math.pi, name='Internal Ring')
        #External Ring

        cbe1, be1 = self.external_ring_contour()
        if cbe1:
            ERC = vm.wires.Contour2D([cbe1]+be1.primitives)
        else:
            ERC = vm.wires.Contour2D(be1.primitives)
        erc = primitives3d.RevolvedProfile(center, axis, z, ERC, center,
                                           axis, angle=2*math.pi, name='External Ring')
        #roller
        ROL = self.rolling_contour_cad()

        # radius=self.F/2.+self.slack+self.Dw/2.
        # rollers=[]
        # theta=2*math.pi/self.Z

        # for zi in range(int(self.Z)):
        #     center_roller = center + radius*math.cos(zi*theta) * y + radius*math.sin(zi*theta) * z

        #     axis_2=vm.Vector3D(1,0,0)
        #     rollers.append(primitives3d.RevolvedProfile(center_roller, axis_2, z, ROL,
        #                                             center_roller, axis_2,
        #                                             angle=2*math.pi,name='Roller {}'.format(zi+1)))

        volumes = [irc, erc]
        return volumes

    def FreeCADExport(self, fcstd_filepath, python_path='python',
                      freecad_lib_path='/usr/lib/freecad/lib', export_types=['fcstd']):
        model = self.VolumeModel()
#        tolerance = self.D/130.
        model.FreeCADExport(fcstd_filepath, python_path=python_path,
                            freecad_lib_path=freecad_lib_path,
                            export_types=export_types)

    def plot_data_quote(self, pos=0):
        delta_quote = 0.05*self.B
        plot_data = []
        #internal diameter
        quote_x = 1.1*self.B/2.
        line1 = vm.edges.LineSegment2D(vm.Point2D(0, self.d/2.), vm.Point2D(quote_x + delta_quote, self.d/2.))
        line1.translation(vm.Vector2D(pos, 0))
        li_data = [line1.plot_data(color=(0, 0, 0), stroke_width=0.1, dash=True)]
        line2 = vm.edges.LineSegment2D(vm.Point2D(0, -self.d/2.), vm.Point2D(quote_x + delta_quote, -self.d/2.))
        line2.translation(vm.Vector2D(pos, 0))
        li_data.append(line2.plot_data(color=(0, 0, 0), stroke_width=0.1, dash=True))
        line3 = vm.edges.LineSegment2D(vm.Point2D(quote_x, self.d/2.), vm.Point2D(quote_x, -self.d/2.))
        line3.translation(vm.Vector2D(pos, 0))
        li_data.append(line3.plot_data(color=(0, 0, 0), stroke_width=0.1, dash=False, marker='triangle_quote'))

        pt_data = {}
        pt_data['fill'] = None
        pt_data['name'] = 'internal diameter'
        pt_data['type'] = 'quote'
        pt_data['label'] = str(round(self.d * 1000, 2)) + ' mm'
        pt_data['x_label'] = quote_x + pos
        pt_data['y_label'] = 0.
        pt_data['rot_label'] = 90
        pt_data['orient_label'] = 'v'
        pt_data['plot_data'] = li_data
        plot_data.append(pt_data)
        #external diameter
        quote_x = -1.3*self.B/2.
        line1 = vm.edgesLineSegment2D(vm.Point2D(0, self.D/2.), vm.Point2D(quote_x - delta_quote, self.D/2.))
        line1.translation(vm.Vector2D(pos, 0))
        li_data = [line1.plot_data(color=(0, 0, 0), stroke_width=0.1, dash=True)]
        line2 = vm.edges.LineSegment2D(vm.Point2D(0, -self.D/2.), vm.Point2D(quote_x - delta_quote, -self.D/2.))
        line2.translation(vm.Vector2D(pos, 0))
        li_data.append(line2.plot_data(color=(0, 0, 0), stroke_width=0.1, dash=True))
        line3 = vm.edges.LineSegment2D(vm.Point2D(quote_x, self.D/2.), vm.Point2D(quote_x, -self.D/2.))
        line3.translation(vm.Vector2D(pos, 0))
        li_data.append(line3.plot_data(color=(0, 0, 0), stroke_width=0.1, dash=False, marker='triangle_quote'))

        pt_data = {}
        pt_data['fill'] = None
        pt_data['name'] = 'external diameter'
        pt_data['type'] = 'quote'
        pt_data['label'] = str(round(self.D * 1000, 2)) + ' mm'
        pt_data['x_label'] = quote_x + pos
        pt_data['y_label'] = 0.
        pt_data['rot_label'] = -90
        pt_data['orient_label'] = 'v'
        pt_data['plot_data'] = li_data
        plot_data.append(pt_data)
        #width
        quote_x = 1.1*self.B/2.
        line1 = vm.edges.LineSegment2D(vm.Point2D(-self.B/2., -self.D/2.), vm.Point2D(-self.B/2., -self.D/2. - quote_x - delta_quote))
        line1.translation(vm.Vector2D(pos, 0))
        li_data = [line1.plot_data(color=(0, 0, 0), stroke_width=0.1, dash=True)]
        line2 = vm.edges.LineSegment2D(vm.Point2D(self.B/2., -self.D/2.), vm.Point2D(self.B/2., -self.D/2. - quote_x - delta_quote))
        line2.translation(vm.Vector2D(pos, 0))
        li_data.append(line2.plot_data(color=(0, 0, 0), stroke_width=0.1, dash=True))
        line3 = vm.edges.LineSegment2D(vm.Point2D(self.B/2., -self.D/2. - quote_x), vm.Point2D(-self.B/2., -self.D/2. - quote_x))
        line3.translation(vm.Vector2D(pos, 0))
        li_data.append(line3.plot_data(color=(0, 0, 0), stroke_width=0.1, dash=False, marker='triangle_quote'))

        pt_data = {}
        pt_data['fill'] = None
        pt_data['name'] = 'width'
        pt_data['type'] = 'quote'
        pt_data['label'] = str(round(self.B * 1000, 2)) + ' mm'
        pt_data['x_label'] = 0 + pos
        pt_data['y_label'] = -self.D/2. - quote_x
        pt_data['rot_label'] = 0
        pt_data['orient_label'] = 'h'
        pt_data['plot_data'] = li_data
        plot_data.append(pt_data)
        return plot_data

    # def plot(self, direction=1, a=None, typ=None):
    #     bg = self.plot_contour(direction)
    #     if a is None:
    #         f, a = bg.MPLPlot(color = 'k')
    #     else:
    #         bg.MPLPlot(a,color = '-k')

    #     if typ == 'Graph':
    #         graph = self.PlotGraph()
    #         graph.MPLPlot(a, 'b', True)

    #     elif typ == 'Load':
    #         self.PlotLoad(a)

    def to_shaft(self):
        return shafts_assembly.Shaft(self.plot_contour(), name=self.name)

    def to_markdown(self):
        infos = ''
        for info in self.infos:
            if info:
                infos += info
        return info



class RadialBallBearing(RadialBearing):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['h', 'E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _non_hash_attributes = ['h', 'E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _generic_eq = True

    symmetric = True
    taking_loads = 'both'
    generate_axial_load = False
    linkage = 'ball'
    coeff_baselife = 3.
    class_name = 'RadialBallBearing'
    cost_coefficient = 0.2
    cost_constant = 1

    def __init__(self, d: float, D: float, B: float, i: int = 1, Z: int = None, Dw: float = None,
                 Cr: float = None, C0r: float = None,
                 material: Material = material_iso,
                 contact_type_point: bool = True, contact_type_linear: bool = False, contact_type_mixed: bool = False,
                 mass: float = None, width: float = 0.02,
                 oil_speed_limit: float = None,
                 grease_speed_limit: float = None,
                 name: str = '', infos: List[str] = None, speed_infos: List[float] = None):

        RadialBearing.__init__(self, d, D, B, alpha=0, i=i, Z=Z, Dw=Dw, Cr=Cr,
                               C0r=C0r, material=material,
                               contact_type_point=contact_type_point, contact_type_linear=contact_type_linear, contact_type_mixed=contact_type_mixed,
                               mass=mass, width=width,
                               oil_speed_limit=oil_speed_limit,
                               grease_speed_limit=grease_speed_limit,
                               name=name, infos=infos, speed_infos=speed_infos)


        # estimation for the graph 2D description
        h1 = self.Dw/2. - (self.E - self.D1)/2.
        self.h = self.B/2. - self.Dw/2.*math.sin(math.acos(h1/(self.Dw/2.))) - 1e-4


    def equivalent_static_load(self, fr, fa=None):
        #Charge radiale statique équivalente
        X0 = 0.6
        Y0 = 0.5
        P0r = max(fr, X0*fr+Y0*fa)
        return P0r

    @classmethod
    def estimate_equivalent_dynamic_load(cls, fr):
        Pr = fr
        return Pr

    def equivalent_dynamic_load(self, fr, fa=0):

        alphap = fsolve((lambda alphap: math.cos(5/180.*math.pi)/math.cos(alphap) \
                        -(1.+0.012534*(abs(fa)/(self.i*self.Z*((self.Dw*1e3)**2) \
                        *math.sin(alphap)))**(2/3.))), self.alpha + 1.)[0]


#        alphap = 0.001
        ksi = 1.05
        nu = 1-math.sin(5/180.*math.pi)/2.5
        e = ksi*math.tan(alphap)
        X1 = 1-0.4*ksi/nu
        X3 = 1
        Y1 = 0.4/nu*1/math.tan(alphap)
        Y3 = 0
        X2 = 1-0.4*ksi/nu
        Y2 = Y1
        if self.i == 1:
            if fa <= e*fr:
                Pr = fr
            else:
                Pr = X1*fr+Y1*fa
        elif self.i == 2:
            if fa <= e*fr:
                Pr = X3*fr+Y3*fa
            else:
                Pr = X2*fr+Y2*fa

        return Pr

    def a_iso(self, kappa, ec, Cu, Pr):
        if kappa < 0.4:
            f = lambda coeff: (1-(2.5671-2.2649/(kappa**0.054381))**(0.83)*((coeff)**(1/3.)))
        elif kappa < 1:
            f = lambda coeff: (1-(2.5671-1.9987/(kappa**0.019087))**(0.83)*((coeff)**(1/3.)))
        else:
            f = lambda coeff: (1-(2.5671-1.9987/(kappa**0.071739))**(0.83)*((coeff)**(1/3.)))
        coeff0 = fsolve(f, ec*Cu/Pr)[0]
        coeff = min(coeff0, ec*Cu/Pr)
        a_iso = 0.1*(f(coeff)**(-9.3))
        return a_iso

    def internal_ring_contour(self):

        pbi2 = vm.Point2D(-self.B/2., self.d1/2.)
        pbi2v = pbi2.to_vector()
        pbi1v = pbi2v.translation(vm.Vector2D(self.h, 0))
        pbi1 = pbi1v.to_point()
        pbi3 = vm.Point2D(-self.B/2., self.d/2.)
        pbi4 = vm.Point2D(self.B/2., self.d/2.)
        pbi5 = vm.Point2D(self.B/2., self.d1/2.)
        pbi6 = pbi5.to_vector().translation(vm.Vector2D(-self.h, 0)).to_point()
        bi1 = primitives2d.OpenedRoundedLineSegments2D([pbi6, pbi5, pbi4, pbi3, pbi2, pbi1],
                                                       {1: self.radius,
                                                        2: self.radius,
                                                        3: self.radius,
                                                        4: self.radius},
                                                       adapt_radius=True)
        cbi1 = vm.edges.Arc2D(pbi1, vm.Point2D(0, self.F/2), pbi6)

        return cbi1, bi1

    def external_ring_contour(self):

        pbe2 = vm.Point2D(-self.B/2., self.D1/2.)
        pbe1 = pbe2.to_vector().translation(vm.Vector2D(self.h, 0)).to_point()
        pbe3 = vm.Point2D(-self.B/2., self.D/2.)
        pbe4 = vm.Point2D(self.B/2., self.D/2.)
        pbe5 = vm.Point2D(self.B/2., self.D1/2.)
        pbe6 = pbe5.to_vector().translation(vm.Vector2D(-self.h, 0)).to_point()


        be1 = primitives2d.OpenedRoundedLineSegments2D([pbe6, pbe5, pbe4, pbe3, pbe2, pbe1],
                                                       {1: self.radius,
                                                        2: self.radius,
                                                        3: self.radius,
                                                        4: self.radius},
                                                       adapt_radius=True)
        cbe1 = vm.edges.Arc2D(pbe1, vm.Point2D(0, self.E/2), pbe6)


        return cbe1, be1

    def rolling_contour(self):

        p0 = vm.Point2D(0, 0)
        c1 = vm.wires.Circle2D(p0, self.Dw/2.)
        return c1

    def rolling_contour_cad(self):
        p0 = vm.Point2D(-self.Dw/2., 0)
        p1 = vm.Point2D(0, self.Dw/2.)
        p2 = vm.Point2D(self.Dw/2., 0)
        a1 = vm.edges.Arc2D(p0, p1, p2)
        l1 = vm.edges.LineSegment2D(p2, p0)
#        c1 = vm.Circle2D(p0, self.Dw/2.)
        return a1

    def plot_contour(self, direction=1):

        be_sup = self.external_ring_contour()
        bi_sup = self.internal_ring_contour()
        ball_sup = self.rolling_contour()
        ball_sup.translation(vm.Vector2D(0, self.Dpw/2.))

        bearing_sup = vm.wires.Contour2D([be_sup, bi_sup, ball_sup])
        bearing_inf = bearing_sup.rotation(vm.Point2D(0, 0), math.pi, True)

        bg = vm.wires.Contour2D([bearing_sup, bearing_inf])
        return bg

    def plot_data(self, pos=0, pos_y=0, stroke_width=1):
        plot_datas = []


        be_arc2d, be_arc2d_2 = self.external_ring_contour()
        be_arc2d1 = be_arc2d.translation((pos, pos_y), True)
        be_arc2d1_2 = be_arc2d_2.translation((pos, pos_y), True)
        plot_datas.append(be_arc2d1.plot_data())
        plot_datas.extend(be_arc2d1_2.plot_data())
        bi_arc2d, bi_arc2d_2 = self.internal_ring_contour()

        bi_arc2d1 = bi_arc2d.translation((pos, pos_y), True)
        bi_arc2d1_2 = bi_arc2d_2.translation((pos, pos_y), True)
        plot_datas.append(bi_arc2d1.plot_data())
        plot_datas.extend(bi_arc2d1_2.plot_data())
        ball_sup = self.rolling_contour()
        ball_sup1 = ball_sup.translation((pos, self.Dpw/2.+pos_y), True)
        plot_datas.append(ball_sup1.plot_data())

        be_arc2d_inf = be_arc2d.rotation(vm.Point2D(0, 0), math.pi, True)
        be_arc2d_inf1 = be_arc2d_inf.translation(vm.Vector2D(pos, pos_y), True)
        plot_datas.append(be_arc2d_inf1.plot_data())
        be_arc2d_inf_2 = be_arc2d_2.rotation(vm.Point2D(0, 0), math.pi, True)
        be_arc2d_inf1_2 = be_arc2d_inf_2.translation(vm.Vector2D(pos, pos_y), True)
        plot_datas.extend(be_arc2d_inf1_2.plot_data())

        bi_arc2d_inf = bi_arc2d.rotation(vm.Point2D(0, 0), math.pi, True)
        bi_arc2d_inf1 = bi_arc2d_inf.translation(vm.Vector2D(pos, pos_y), True)
        plot_datas.append(bi_arc2d_inf1.plot_data())

        bi_arc2d_inf_2 = bi_arc2d_2.rotation(vm.Point2D(0, 0), math.pi, True)
        bi_arc2d_inf1_2 = bi_arc2d_inf_2.translation(vm.Vector2D(pos, pos_y), True)
        plot_datas.extend(bi_arc2d_inf1_2.plot_data())


        ball_inf1 = ball_sup1.rotation(vm.Point2D(pos, pos_y), math.pi, True)
        plot_datas.append(ball_inf1.plot_data())

        return [plot_data.PrimitiveGroup(plot_datas)]

    @classmethod
    def graph(cls, list_node, direction=1):

        graph = nx.DiGraph()
        graph.add_edges_from([(list_node[4], list_node[0])])
        graph.add_edges_from([(list_node[1], list_node[5])])
        graph.add_edges_from([(list_node[6], list_node[2])])
        graph.add_edges_from([(list_node[3], list_node[7])])

        graph.add_edges_from([(list_node[4], list_node[2])])
        graph.add_edges_from([(list_node[3], list_node[5])])
        graph.add_edges_from([(list_node[1], list_node[7])])
        graph.add_edges_from([(list_node[6], list_node[0])])

        return graph

#    @classmethod
#    def dict_to_object(cls, d):
#        if 'Cr' not in d.keys():
#            d['Cr'] = None
#        if 'C0r' not in d.keys():
#            d['C0r'] = None
#        obj = cls(d = d['d'], D = d['D'], B = d['B'], i = d['i'], Z = d['Z'],
#                  Dw = d['Dw'], Cr = d['Cr'], C0r = d['C0r'],
#                  material = Material.dict_to_object(d['material']),
#                  contact_type = d['contact_type'],
#                  name=d['name'], mass=d['mass'])
#        return obj

#    def Copy(self):
#        if not hasattr(self, 'Cr'):
#            Cr = None
#        else:
#            Cr = self.Cr
#        if not hasattr(self, 'C0r'):
#            C0r = None
#        else:
#            C0r = self.C0r
#        obj = RadialBallBearing(d = self.d, D = self.D, B = self.B, i = self.i, Z = self.Z,
#                  Dw = self.Dw, Cr = Cr, C0r = C0r,
#                  oil = self.oil, material = self.material,
#                  contact_type = self.contact_type,
#                  name = self.name)
#        return obj

class AngularBallBearing(RadialBearing):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['h1', 'h2', 'D2', 'd2', 'E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _non_hash_attributes = ['h1', 'h2', 'D2', 'd2', 'E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _generic_eq = True

    symmetric = False
    taking_loads = 'right'
    generate_axial_load = True
    linkage = 'ball_joint'
    coeff_baselife = 3.
    class_name = 'AngularBallBearing'
    cost_coefficient = 0.4
    cost_constant = 1.5

    def __init__(self, d: float, D: float, B: float, alpha: float, i: int = 1, Z: int = None,
                 Dw: float = None, Cr: float = None, C0r: float = None,
                 material: Material = material_iso,
                 contact_type_point: bool = True, contact_type_linear: bool = False,
                 contact_type_mixed: bool = False,
                 mass: float = None,
                 width: float = None,
                 oil_speed_limit: float = None,
                 grease_speed_limit: float = None,
                 name: str = '', infos: List[str] = None, speed_infos: List[float] = None):
        RadialBearing.__init__(self, d, D, B, alpha=alpha, i=1, Z=Z, Dw=Dw, Cr=Cr,
                               C0r=C0r, material=material,
                               contact_type_point=contact_type_point, contact_type_linear=contact_type_linear, contact_type_mixed=contact_type_mixed,
                               mass=mass, width=width,
                               oil_speed_limit=oil_speed_limit,
                               grease_speed_limit=grease_speed_limit,
                               name=name, infos=infos, speed_infos=speed_infos)



        # estimation for the graph 2D description
        h1 = self.Dw/2. - (self.E - self.D1)/2.
        self.h1 = self.B/2. - self.Dw/2.*math.sin(math.acos(h1/(self.Dw/2.))) - 1e-4
        h2 = 0.95*self.Dw/2.
        self.h2 = self.B/2. - self.Dw/2.*math.sin(math.acos(h2/(self.Dw/2.))) - 1e-4
        self.D2 = 0.6*(self.D - self.E) + self.E
        self.d2 = 0.6*(self.F - self.d) + self.d

    def equivalent_static_load(self, fr, fa=None):
        #Charge radiale statique équivalente
        if self.i == 1:
            X0 = 0.5
            evol_Y0 = [0.52, 0.5, 0.46, 0.42, 0.38, 0.33, 0.29, 0.26, 0.22]
            evol_alpha = npy.array([5, 10, 15, 20, 25, 30, 35, 40, 45])/180*math.pi
            f = interpolate.interp1d(list(evol_alpha), evol_Y0, fill_value='extrapolate')
            Y0 = float(f(self.alpha))
        elif self.i == 2:
            X0 = 1
            evol_Y0 = [1.04, 1, 0.92, 0.84, 0.76, 0.66, 0.58, 0.52, 0.44]
            evol_alpha = npy.array([5, 10, 15, 20, 25, 30, 35, 40, 45])/180*math.pi
            f = interpolate.interp1d(list(evol_alpha), evol_Y0, fill_value='extrapolate')
            Y0 = float(f(self.alpha))
        P0r = max(fr, X0*fr+Y0*fa)
        return P0r

    @classmethod
    def estimate_equivalent_dynamic_load(cls, fr):
        Pr = fr
        return Pr


    def equivalent_dynamic_load(self, fr, fa=0):


        alphap = fsolve((lambda alphap: math.cos(self.alpha)/math.cos(alphap) \
                    -(1+0.012534*(abs(fa)/(self.i*self.Z*((self.Dw*1e3)**2)*math.sin(alphap)))**(2/3.))), self.alpha + 1)[0]


        if self.alpha <= 5/180.*math.pi:
            if self.i == 1:
                ksi = 1.05
            elif self.i == 2:
                ksi = 1.25
        else:
            ksi = 1.25
        if self.alpha <= 5/180.*math.pi:
            nu = 1-math.sin(5/180.*math.pi)/2.5
        elif self.alpha <= 15/180.*math.pi:
            nu = 1-math.sin(self.alpha)/2.5
        else:
            nu = 1-math.sin(self.alpha)/2.75
        if self.alpha <= 15/180.*math.pi:
            e = ksi*math.tan(alphap)
        X1 = 1-0.4*ksi/nu
        X3 = 1
        if self.alpha <= 15/180.*math.pi:
            Y1 = 0.4/nu*1/math.tan(alphap)
        else:
            alpha1 = math.acos(math.cos(self.alpha)*0.9724)
            Y1 = fsolve((lambda Y1: Y1-(0.4/math.tan(alpha1))/(1-(1/3.)*math.sin(alpha1))), 1)[0]
        if self.alpha > 15/180.*math.pi:
            e = (1-X1)/Y1
            Y3 = 0.625/e
        elif self.alpha <= 15/180.*math.pi:
            Y3 = 0.625/ksi*(1/math.tan(alphap))
        X2 = 1.625*X1
        Y2 = 1.625*Y1


        if self.i == 1:
            if abs(fa) <= e*fr:

                Pr = fr
            else:


                Pr = X1*fr+Y1*abs(fa)
        elif self.i == 2:
            if abs(fa) <= e*fr:

                Pr = X3*fr+Y3*abs(fa)
            else:

                Pr = X2*fr+Y2*abs(fa)

        return Pr

    def a_iso(self, kappa, ec, Cu, Pr):
        if kappa < 0.4:
            f = lambda coeff: (1-(2.5671-2.2649/(kappa**0.054381))**(0.83)*((coeff)**(1/3.)))
        elif kappa < 1:
            f = lambda coeff: (1-(2.5671-1.9987/(kappa**0.019087))**(0.83)*((coeff)**(1/3.)))
        else:
            f = lambda coeff: (1-(2.5671-1.9987/(kappa**0.071739))**(0.83)*((coeff)**(1/3.)))
        coeff0 = fsolve(f, ec*Cu/Pr)[0]
        coeff = min(coeff0, ec*Cu/Pr)
        a_iso = 0.1*(f(coeff)**(-9.3))
        return a_iso

    def internal_ring_contour(self, direction=1, sign_V=1):

        pbi2 = vm.Point2D(direction*self.B/2., sign_V*self.d2/2.)
        pbi1 = vm.Point2D(direction*(self.B/2. - self.h2), sign_V*(self.Dpw/2. - self.Dw/2.*0.95))
        pbi3 = vm.Point2D(direction*self.B/2., sign_V*self.d/2.)
        pbi4 = vm.Point2D(-direction*self.B/2., sign_V*self.d/2.)
        pbi5 = vm.Point2D(-direction*self.B/2., sign_V*self.d1/2.)
        pbi6 = pbi5.to_vector().translation(vm.Vector2D(direction*self.h1, 0)).to_point()
        bi1 = primitives2d.OpenedRoundedLineSegments2D([pbi6, pbi5, pbi4, pbi3, pbi2, pbi1],
                                                       {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius},
                                                       adapt_radius=True)

        cbi1 = vm.edges.Arc2D(pbi1, vm.Point2D(0, sign_V*self.F/2), pbi6)

        return cbi1, bi1

    def external_ring_contour(self, direction=1, sign_V=1):

        pbe2 = vm.Point2D(direction*self.B/2., sign_V*self.D1/2.)
        pbe1 = pbe2.to_vector().translation(vm.Vector2D(-direction*self.h1, 0)).to_point()
        pbe3 = vm.Point2D(direction*self.B/2., sign_V*self.D/2.)
        pbe4 = vm.Point2D(-direction*self.B/2., sign_V*self.D/2.)
        pbe5 = vm.Point2D(-direction*self.B/2., sign_V*self.D2/2.)
        pbe6 = vm.Point2D(-direction*(self.B/2. - self.h2), sign_V*(self.Dpw/2. + self.Dw/2.*0.95))
        be1 = primitives2d.OpenedRoundedLineSegments2D([pbe1, pbe2, pbe3, pbe4, pbe5, pbe6],
                                                       {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius},
                                                       adapt_radius=True)

        cbe1 = vm.edges.Arc2D(pbe6, vm.Point2D(0, sign_V*self.E/2), pbe1)
        return cbe1, be1


    def rolling_contour(self):

        p0 = vm.Point2D(0, 0)
        c1 = vm.wires.Circle2D(p0, self.Dw/2.)
        return c1

    def rolling_contour_cad(self):

        p0 = vm.Point2D(-self.Dw/2., 0)
        p1 = vm.Point2D(0, self.Dw/2.)
        p2 = vm.Point2D(self.Dw/2., 0)
        a1 = vm.edges.Arc2D(p0, p1, p2)
        l1 = vm.edges.LineSegment2D(p2, p0)
#        c1 = vm.Circle2D(p0, self.Dw/2.)
        return vm.wires.Contour2D([a1, l1])

    def plot_contour(self, direction=1):

        be_sup = self.external_ring_contour(direction=direction, sign_V=1)
        be_inf = self.external_ring_contour(direction=direction, sign_V=-1)
        bi_sup = self.internal_ring_contour(direction=direction, sign_V=1)
        bi_inf = self.internal_ring_contour(direction=direction, sign_V=-1)
        ball = self.rolling_contour()
        ball_sup = ball.translation(vm.Vector2D(0, self.Dpw/2.), True)
        ball_inf = ball.translation(vm.Vector2D(0, -self.Dpw/2.), True)
        bg = vm.wires.Contour2D([be_sup, bi_sup, ball_sup, be_inf, bi_inf, ball_inf])
        return bg

    def plot_data(self, pos=0, pos_y=0, quote=True, constructor=True, direction=1):

        plot_datas = []
        be_sup, be_sup_2 = self.external_ring_contour(direction=direction, sign_V=1)
        be_sup1 = be_sup.translation((pos, pos_y), True)
        plot_datas.append(be_sup1.plot_data())
        be_sup1_2 = be_sup_2.translation((pos, pos_y), True)
        plot_datas.extend(be_sup1_2.plot_data())


        bi_sup, bi_sup_2 = self.internal_ring_contour(direction=direction, sign_V=1)
        bi_sup1 = bi_sup.translation((pos, pos_y), True)
        plot_datas.append(bi_sup1.plot_data())
        bi_sup1_2 = bi_sup_2.translation((pos, pos_y), True)
        plot_datas.extend(bi_sup1_2.plot_data())

        ball = self.rolling_contour()
        ball_sup = ball.translation((0, self.Dpw/2.), True)
        ball_sup1 = ball_sup.translation((pos, pos_y), True)
        plot_datas.append(ball_sup1.plot_data())

        be_inf, be_inf_2 = self.external_ring_contour(direction=direction, sign_V=-1)
        be_inf1 = be_inf.translation((pos, pos_y), True)
        plot_datas.append(be_inf1.plot_data())

        be_inf1_2 = be_inf_2.translation((pos, pos_y), True)
        plot_datas.extend(be_inf1_2.plot_data())

        bi_inf, bi_inf_2 = self.internal_ring_contour(direction=direction, sign_V=-1)
        bi_inf1 = bi_inf.translation((pos, pos_y), True)
        plot_datas.append(bi_inf1.plot_data())

        bi_inf1_2 = bi_inf_2.translation((pos, pos_y), True)
        plot_datas.extend(bi_inf1_2.plot_data())

        ball_inf = ball.translation((0, -self.Dpw/2.), True)
        ball_inf1 = ball_inf.translation((pos, pos_y), True)
        plot_datas.append(ball_inf1.plot_data())

        return [plot_data.PrimitiveGroup(plot_datas)]

    def volmdlr_primitives(self, center=vm.O3D, axis=vm.X3D, direction=1):
        # TODO: mutualization of this in parent class?
        axis.normalize()

        y = axis.random_unit_normal_vector()
        z = axis.cross(y)

        #Internal Ring
        cbi1, bi1 = self.internal_ring_contour(direction=direction)
        if cbi1:
            IRC = vm.wires.Contour2D([cbi1]+bi1.primitives)
        else:
            IRC = vm.wires.Contour2D(bi1.primitives)
        irc = primitives3d.RevolvedProfile(center, axis, z, IRC, center,
                                           axis, angle=2*math.pi, name='Internal Ring')
        #External Ring

        cbe1, be1 = self.external_ring_contour(direction=direction)
        if cbe1:
            ERC = vm.wires.Contour2D([cbe1]+be1.primitives)
        else:
            ERC = vm.wires.Contour2D(be1.primitives)
        erc = primitives3d.RevolvedProfile(center, axis, z, ERC, center,
                                           axis, angle=2*math.pi, name='External Ring')
        #roller
        ROL = self.rolling_contour_cad()

        # radius=self.F/2.+self.slack+self.Dw/2.
        # rollers=[]
        # theta=2*math.pi/self.Z

        # for zi in range(int(self.Z)):
        #     center_roller = center + radius*math.cos(zi*theta) * y + radius*math.sin(zi*theta) * z

        #     axis_2=vm.Vector3D(1,0,0)
        #     rollers.append(primitives3d.RevolvedProfile(center_roller, axis_2, z, ROL,
        #                                             center_roller, axis_2,
        #                                             angle=2*math.pi,name='Roller {}'.format(zi+1)))

        volumes = [irc, erc]
        return volumes

    @classmethod
    def graph(cls, list_node, direction=1):

        graph = nx.DiGraph()
        graph.add_edges_from([(list_node[4], list_node[0])])
        graph.add_edges_from([(list_node[1], list_node[5])])
        graph.add_edges_from([(list_node[6], list_node[2])])
        graph.add_edges_from([(list_node[3], list_node[7])])

        if direction == 1:
            graph.add_edges_from([(list_node[4], list_node[2])])
            graph.add_edges_from([(list_node[3], list_node[5])])
        elif direction == -1:
            graph.add_edges_from([(list_node[1], list_node[7])])
            graph.add_edges_from([(list_node[6], list_node[0])])

        return graph

#    @classmethod
#    def dict_to_object(cls, d):
#        if 'Cr' not in d.keys():
#            d['Cr'] = None
#        if 'C0r' not in d.keys():
#            d['C0r'] = None
#        obj = cls(d = d['d'], D = d['D'], B = d['B'], alpha = d['alpha'], i = d['i'], Z = d['Z'],
#                  Dw = d['Dw'], Cr = d['Cr'], C0r = d['C0r'],
#                  material = Material.dict_to_object(d['material']),
#                  contact_type = d['contact_type'],
#                  name=d['name'], mass=d['mass'])
#        return obj

#    def Copy(self):
#        if not hasattr(self, 'Cr'):
#            Cr = None
#        else:
#            Cr = self.Cr
#        if not hasattr(self, 'C0r'):
#            C0r = None
#        else:
#            C0r = self.C0r
#        obj = AngularBallBearing(d = self.d, D = self.D, B = self.B, alpha = self.alpha,
#                            i = self.i, Z = self.Z,
#                            Dw = self.Dw, Cr = Cr, C0r = C0r,
#                            oil = self.oil, material = self.material,
#                            contact_type = self.contact_type,
#                            name = self.name)
#        return obj

class SphericalBallBearing(RadialBearing):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _non_hash_attributes = ['E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _generic_eq = True

    symmetric = True
    taking_loads = 'both'
    generate_axial_load = False
    linkage = 'ball_joint'
    coeff_baselife = 3.
    class_name = 'SphericalBallBearing'
    cost_coefficient = 0.4
    cost_constant = 2

    def __init__(self, d: float, D: float, B: float, alpha: float = 0, i: int = 1, Z: int = None,
                 Dw: float = None, Cr: float = None, C0r: float = None,
                 material: Material = material_iso,
                 contact_type_point: bool = True, contact_type_linear: bool = False, contact_type_mixed: bool = False,
                 mass: float = None, width: float = 0.02,
                 oil_speed_limit: float = None,
                 grease_speed_limit: float = None,
                 name: str = '', infos: List[str] = None, speed_infos: List[float] = None):
        RadialBearing.__init__(self, d, D, B, alpha, i, Z, Dw, Cr, C0r,
                               material,
                               contact_type_point=contact_type_point, contact_type_linear=contact_type_linear, contact_type_mixed=contact_type_mixed,
                               mass=mass, width=width,
                               oil_speed_limit=oil_speed_limit,
                               grease_speed_limit=grease_speed_limit,
                               name=name, infos=infos, speed_infos=speed_infos)



    def equivalent_static_load(self, fr, fa=None):
        #Charge radiale statique équivalente
        if self.i == 1:
            X0 = 0.5
            Y0 = 0.22/math.tan(self.alpha)
        elif self.i == 2:
            X0 = 1
            Y0 = 0.44/math.tan(self.alpha)
        P0r = max(fr, X0*fr+Y0*fa)
        return P0r

    @classmethod
    def estimate_equivalent_dynamic_load(cls, fr):
        Pr = fr
        return Pr

    def equivalent_dynamic_load(self, fr, fa=0):
        alphap = self.alpha
        ksi = 1.5
        nu = 1
        e = ksi*math.tan(alphap)
        X1 = 1-0.4*ksi/nu
        X3 = 1
        Y1 = 0.4/nu*(1/math.tan(alphap))
        Y3 = 0.625/ksi*(1/math.tan(alphap))
        X2 = 1.625*X1
        Y2 = 1.625*Y1
        if self.i == 1:
            if abs(fa) <= e*fr:
                Pr = fr
            else:
                Pr = X1*fr+Y1*abs(fa)
        elif self.i == 2:
            if abs(fa) <= e*fr:
                Pr = X3*fr+Y3*abs(fa)
            else:
                Pr = X2*fr+Y2*abs(fa)
        return Pr

    def a_iso(self, kappa, ec, Cu, Pr):
        if kappa < 0.4:
            f = lambda coeff: (1-(2.5671-2.2649/(kappa**0.054381))**(0.83)*((coeff)**(1/3.)))
        elif kappa < 1:
            f = lambda coeff: (1-(2.5671-1.9987/(kappa**0.019087))**(0.83)*((coeff)**(1/3.)))
        else:
            f = lambda coeff: (1-(2.5671-1.9987/(kappa**0.071739))**(0.83)*((coeff)**(1/3.)))
        coeff0 = fsolve(f, ec*Cu/Pr)[0]
        coeff = min(coeff0, ec*Cu/Pr)
        a_iso = 0.1*(f(coeff)**(-9.3))
        return a_iso

    @classmethod
    def graph(cls, list_node, direction=1):

        graph = nx.DiGraph()
        graph.add_edges_from([(list_node[4], list_node[0])])
        graph.add_edges_from([(list_node[1], list_node[5])])
        graph.add_edges_from([(list_node[6], list_node[2])])
        graph.add_edges_from([(list_node[3], list_node[7])])

        graph.add_edges_from([(list_node[4], list_node[2])])
        graph.add_edges_from([(list_node[3], list_node[5])])
        graph.add_edges_from([(list_node[1], list_node[7])])
        graph.add_edges_from([(list_node[6], list_node[0])])

        return graph

#    @classmethod
#    def dict_to_object(cls, d):
#        if 'Cr' not in d.keys():
#            d['Cr'] = None
#        if 'C0r' not in d.keys():
#            d['C0r'] = None
#        obj = cls(d = d['d'], D = d['D'], B = d['B'], alpha = d['alpha'],
#                  i = d['i'], Z = d['Z'],
#                  Dw = d['Dw'], Cr = d['Cr'], C0r = d['C0r'],
#                  material = Material.dict_to_object(d['material']),
#                  contact_type = d['contact_type'],
#                  name=d['name'], mass=d['mass'])
#        return obj

#    def Copy(self):
#        if not hasattr(self, 'Cr'):
#            Cr = None
#        else:
#            Cr = self.Cr
#        if not hasattr(self, 'C0r'):
#            C0r = None
#        else:
#            C0r = self.C0r
#        obj = SphericalBallBearing(d = self.d, D = self.D, B = self.B, alpha = self.alpha,
#                            i = self.i, Z = self.Z,
#                            Dw = self.Dw, Cr = Cr, C0r = C0r,
#                            oil = self.oil, material = self.material,
#                            contact_type = self.contact_type,
#                            name = self.name)
#        return obj

# TODO remove alpha?
class RadialRollerBearing(RadialBearing):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['Dpw', 'Lw', 'h', 'E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _non_hash_attributes = ['Dpw', 'Lw', 'h', 'E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _generic_eq = True

    # symmetric = True
    linkage = 'cylindric'
    coeff_baselife = 10/3.
    cost_coefficient = 0.5
    cost_constant = 2.5

    def __init__(self, d: float, D: float, B: float, alpha: float, i: int = 1, Z: int = None,
                 Dw: float = None, Cr: float = None, C0r: float = None,
                 material: Material = material_iso,
                 contact_type_point: bool = True, contact_type_linear: bool = False,
                 contact_type_mixed: bool = False,
                 mass: float = None, width: float = 0.02,
                 oil_speed_limit: float = None,
                 grease_speed_limit: float = None,
                 name: str = '', infos: List[str] = None, speed_infos: List[float] = None):


        RadialBearing.__init__(self, d, D, B, alpha=alpha, i=1, Z=Z, Dw=Dw,
                               Cr=Cr, C0r=C0r,
                               material=material, contact_type_point=contact_type_point,
                               contact_type_linear=contact_type_linear, contact_type_mixed=contact_type_mixed,
                               mass=mass, width=width,
                               oil_speed_limit=oil_speed_limit,
                               grease_speed_limit=grease_speed_limit,
                               name=name, infos=infos, speed_infos=speed_infos)

#        self.typ = typ

        # estimation for the graph 2D description
        self.Dpw = (self.d + self.D)/2.
        self.Lw = 0.7*self.B
        self.h = self.B/2. - self.Lw/2. - 1e-4

    def equivalent_static_load(self, fr, fa=None):
        #Charge radiale statique équivalente
        if self.alpha != 0:
            x0 = 0.5*self.i
            y0 = 0.22*1/math.tan(self.alpha)*self.i
        else:
            x0 = 1
            y0 = 0
        P0r = max(fr, x0*fr+y0*fa)
        return P0r

    @classmethod
    def estimate_equivalent_dynamic_load(cls, fr):
        Pr = fr
        return Pr

    def equivalent_dynamic_load(self, fr, fa=0):

        ksi = 1.5 #param of the ISO 1281
        e = ksi*math.tan(self.alpha)
        nu = 1-0.15*math.sin(self.alpha)
        w = self.material.weibull_e*self.coeff_baselife

        if self.contact_type_point:
            if self.i == 1:
                Jr0p5 = 0.2288
#                Ja0p5 = 0.2782
                J10p5 = 0.5625
                J20p5 = 0.5875
            else: #analyse if the parameters are true for i>2
                Jr0p5 = 0.4577
#                Ja0p5 = 0.
                J10p5 = 0.6925
                J20p5 = 0.7233
        elif self.contact_type_linear:
            if self.i == 1:
                Jr0p5 = 0.2453
#                Ja0p5 = 0.3090
                J10p5 = 0.6495
                J20p5 = 0.6744
            else: #analyse if the parameters are true for i>2
                Jr0p5 = 0.4906
#                Ja0p5 = 0
                J10p5 = 0.7577
                J20p5 = 0.7867
        elif self.contact_type_mixed:
            if self.i == 1:
                Jr0p5 = 0.2369
                # TODO: check why next variable is unused
#                Ja0p5 = 0.2932
                J10p5 = 0.6044
                J20p5 = 0.6295
            else: #analyse if the parameters are true for i>2
                Jr0p5 = 0.4739
#                Ja0p5 = 0
                J10p5 = 0.7244
                J20p5 = 0.7543

        X1 = 1-Jr0p5/(J10p5*J20p5)**0.5*ksi/nu
        X2 = 2**(1-(1/w))*X1
        X3 = 1
        if self.alpha > 0:
            Y1 = Jr0p5*(1/math.tan(self.alpha))/(J10p5*J20p5)**0.5*1/nu
            Y2 = 2**(1-(1/w))*Y1
            Y3 = 1/ksi*(2**(1-(1/w))-1)*(1/math.tan(self.alpha))

        if self.alpha > 0:
            if self.i == 1:
                if abs(fa) <= e*fr:
                    Pr = fr
                else:
                    Pr = X1*fr+Y1*abs(fa)
            elif self.i == 2:
                if abs(fa) <= e*fr:
                    Pr = X3*fr+Y3*abs(fa)
                else:
                    Pr = X2*fr+Y2*abs(fa)
        else: # for alpha=0 the axial load is not include in the L10
            if self.i == 1:
                Pr = fr
            elif self.i == 2:
                Pr = X3*fr

        return Pr

    def a_iso(self, kappa, ec, Cu, Pr):
        if kappa < 0.4:
            f = lambda coeff: (1-(1.5859-1.3993/(kappa**0.054381))*((coeff)**0.4))
        elif kappa < 1:
            f = lambda coeff: (1-(1.5859-1.2348/(kappa**0.19087))*((coeff)**0.4))
        else:
            f = lambda coeff: (1-(1.5859-1.2348/(kappa**0.071739))*((coeff)**0.4))
        try:
            coeff0 = fsolve(f, ec*Cu/Pr)[0]
        except FloatingPointError:# TODO check this for a better solution
            coeff0 = ec*Cu/Pr
        coeff = min(coeff0, ec*Cu/Pr)
        if f(coeff) > 0:
            a_iso = 0.1*(f(coeff)**(-9.185))
        else:
            return 100.
        return a_iso


    def rolling_contour(self):

        p1 = vm.Point2D(-self.Lw/2., -self.Dw/2.)
        p2 = vm.Point2D(-self.Lw/2., self.Dw/2.)
        p3 = vm.Point2D(self.Lw/2., self.Dw/2.)
        p4 = vm.Point2D(self.Lw/2., -self.Dw/2.)
        rol = primitives2d.ClosedRoundedLineSegments2D([p1, p2, p3, p4],
                                                       {0: self.radius, 1: self.radius, 2: self.radius, 3: self.radius})
        return rol

    def rolling_contour_cad(self):
        p1 = vm.Point2D(-self.Lw/2., 0)
        p2 = vm.Point2D(-self.Lw/2., self.Dw/2.)
        p3 = vm.Point2D(self.Lw/2., self.Dw/2.)
        p4 = vm.Point2D(self.Lw/2., 0)
        rol = primitives2d.ClosedRoundedLineSegments2D([p1, p2, p3, p4], {1: self.radius, 2: self.radius}, True)
        return rol


    def plot_data(self, pos=0, pos_y=0,
                  quote=True, constructor=True, direction=1, stroke_width=1):

        plot_datas = []
        be_sup = self.external_ring_contour(direction=direction, sign_V=1)[1]
        be_sup1 = be_sup.translation((pos, pos_y), True)
        plot_datas.append(be_sup1.plot_data())

        be_inf = self.external_ring_contour(direction=direction, sign_V=-1)[1]
        be_inf1 = be_inf.translation((pos, pos_y), True)
        plot_datas.append(be_inf1.plot_data())

        bi_sup = self.internal_ring_contour(direction=direction, sign_V=1)[1]
        bi_sup1 = bi_sup.translation((pos, pos_y), True)
        plot_datas.append(bi_sup1.plot_data())

        bi_inf = self.internal_ring_contour(direction=direction, sign_V=-1)[1]
        bi_inf1 = bi_inf.translation((pos, pos_y), True)
        plot_datas.append(bi_inf1.plot_data())

        roller = self.rolling_contour()
        roller_sup = roller.translation((0, self.Dpw/2.), True)
        roller_sup1 = roller_sup.translation((pos, pos_y), True)
        plot_datas.append(roller_sup1.plot_data())

        roller_inf = roller.translation((0, -self.Dpw/2.), True)
        roller_inf1 = roller_inf.translation((pos, pos_y), True)
        plot_datas.append(roller_inf1.plot_data())

        return [plot_data.PrimitiveGroup(plot_datas)]

    def plot_contour(self, direction=1):

        be_sup = self.external_ring_contour(direction=direction, sign_V=1)
        be_inf = self.external_ring_contour(direction=direction, sign_V=-1)
        bi_sup = self.internal_ring_contour(direction=direction, sign_V=1)
        bi_inf = self.internal_ring_contour(direction=direction, sign_V=-1)
        roller = self.rolling_contour()
        roller_sup = roller.translation(vm.Vector2D(0, self.Dpw/2.), True)
        roller_inf = roller.translation(vm.Vector2D(0, -self.Dpw/2.), True)

        bg = vm.wires.Contour2D([be_sup, bi_sup, roller_sup, be_inf, bi_inf, roller_inf])
        return bg

    def Copy(self):
        if not hasattr(self, 'Cr'):
            Cr = None
        else:
            Cr = self.Cr
        if not hasattr(self, 'C0r'):
            C0r = None
        else:
            C0r = self.C0r
        if self.class_name == 'N':
            return N(d=self.d, D=self.D, B=self.B, i=self.i, Z=self.Z,
                     Dw=self.Dw, Cr=Cr, C0r=C0r,
                     oil=self.oil, material=self.material,
                     contact_type_point=self.contact_type_point, contact_type_linear=self.contact_type_linear, contact_type_mixed=self.contact_type_mixed,
                     name=self.name)
        elif self.class_name == 'NU':
            return NU(d=self.d, D=self.D, B=self.B, i=self.i, Z=self.Z,
                      Dw=self.Dw, Cr=Cr, C0r=C0r,
                      oil=self.oil, material=self.material,
                      contact_type_point=self.contact_type_point, contact_type_linear=self.contact_type_linear, contact_type_mixed=self.contact_type_mixed,
                      name=self.name)
        elif self.class_name == 'NF':
            return NF(d=self.d, D=self.D, B=self.B, i=self.i, Z=self.Z,
                      Dw=self.Dw, Cr=Cr, C0r=C0r,
                      oil=self.oil, material=self.material,
                      contact_type_point=self.contact_type_point, contact_type_linear=self.contact_type_linear, contact_type_mixed=self.contact_type_mixed,
                      name=self.name)
        elif self.class_name == 'NUP':
            return NUP(d=self.d, D=self.D, B=self.B, i=self.i, Z=self.Z,
                       Dw=self.Dw, Cr=Cr, C0r=C0r,
                       oil=self.oil, material=self.material,
                       contact_type_point=self.contact_type_point, contact_type_linear=self.contact_type_linear, contact_type_mixed=self.contact_type_mixed,
                       name=self.name)

class NUP(RadialRollerBearing):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['Dpw', 'Lw', 'h', 'E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _non_hash_attributes = ['Dpw', 'Lw', 'h', 'E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _generic_eq = True

    symmetric = True
    taking_loads = 'both'
    generate_axial_load = False
    class_name = 'NUP'
    cost_coefficient = 0.5
    cost_constant = 2.7

    # TODO: remove alpha?
    def __init__(self, d: float, D: float, B: float, i: int = 1, Z: int = None, Dw: float = None,
                 Cr: float = None, C0r: float = None,
                 material: Material = material_iso,
                 contact_type_point: bool = False, contact_type_linear: bool = True,
                 contact_type_mixed: bool = False,
                 mass: float = None,
                 oil_speed_limit: float = None,
                 grease_speed_limit: float = None,
                 name: str = '', infos: List[str] = None, speed_infos: List[float] = None):
        RadialRollerBearing.__init__(self, d, D, B, alpha=0, i=i, Z=Z, Dw=Dw, Cr=Cr,
                                     C0r=C0r,
                                     material=material,
                                     contact_type_point=contact_type_point, contact_type_linear=contact_type_linear, contact_type_mixed=contact_type_mixed,
                                     mass=mass,
                                     oil_speed_limit=oil_speed_limit,
                                     grease_speed_limit=grease_speed_limit,
                                     name=name, infos=infos, speed_infos=speed_infos)

    def internal_ring_contour(self, direction=1, sign_V=1):

        d1 = self.d1
        pbi2 = vm.Point2D(-direction*self.B/2., sign_V*d1/2.)
        pbi1 = vm.Point2D(-direction*(self.B/2. - self.h), sign_V*d1/2.)
        pbi0 = vm.Point2D(-direction*(self.B/2. - self.h), sign_V*(self.F/2.))
        pbi3 = vm.Point2D(-direction*self.B/2., sign_V*self.d/2.)
        pbi4 = vm.Point2D(direction*self.B/2., sign_V*self.d/2.)
        pbi5 = vm.Point2D(direction*self.B/2., sign_V*d1/2.)
        pbi6 = vm.Point2D(direction*(self.B/2. - self.h), sign_V*d1/2.)
        pbi7 = vm.Point2D(direction*(self.B/2. - self.h), sign_V*(self.F/2.))
        irc = primitives2d.ClosedRoundedLineSegments2D([pbi0, pbi1, pbi2, pbi3, pbi4, pbi5, pbi6, pbi7],
                                                       {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius,
                                                        5: self.radius, 6: self.radius}, adapt_radius=True)

        # cbi1 = vm.edges.Arc2D(pbi1, vm.Point2D(0, self.F/2), pbi7)



        return None, irc



    def external_ring_contour(self, direction=1, sign_V=1):

        D1 = self.D1
        pbe2 = vm.Point2D(-direction*self.B/2., sign_V*D1/2.)
        pbe1 = vm.Point2D(-direction*(self.B/2. - self.h), sign_V*D1/2.)
        pbe0 = vm.Point2D(-direction*(self.B/2. - self.h), sign_V*(self.E/2.))
        pbe3 = vm.Point2D(-direction*self.B/2., sign_V*self.D/2.)
        pbe4 = vm.Point2D(direction*self.B/2., sign_V*self.D/2.)
        pbe5 = vm.Point2D(direction*self.B/2., sign_V*D1/2.)
        pbe6 = vm.Point2D(direction*(self.B/2. - self.h), sign_V*D1/2.)
        pbe7 = vm.Point2D(direction*(self.B/2. - self.h), sign_V*(self.E/2.))
        be1 = primitives2d.ClosedRoundedLineSegments2D([pbe0, pbe1, pbe2, pbe3, pbe4, pbe5, pbe6, pbe7],
                                                       {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius,
                                                        5: self.radius, 6: self.radius}, adapt_radius=True)


        # erc = vm.Contour2D([be1])
        # cbe1 = vm.edges.Arc2D(pbe1, vm.Point2D(0, self.F/2), pbe7)
        return None, be1

    @classmethod
    def graph(cls, list_node, direction=1):

        graph = nx.DiGraph()
        graph.add_edges_from([(list_node[4], list_node[0])])
        graph.add_edges_from([(list_node[1], list_node[5])])
        graph.add_edges_from([(list_node[6], list_node[2])])
        graph.add_edges_from([(list_node[3], list_node[7])])

        graph.add_edges_from([(list_node[4], list_node[2])])
        graph.add_edges_from([(list_node[3], list_node[5])])
        graph.add_edges_from([(list_node[1], list_node[7])])
        graph.add_edges_from([(list_node[6], list_node[0])])

        return graph


class N(RadialRollerBearing):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['Dpw', 'Lw', 'h', 'E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _non_hash_attributes = ['Dpw', 'Lw', 'h', 'E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _generic_eq = True

    symmetric = True
    taking_loads = 'free'
    generate_axial_load = False
    class_name = 'N'
    cost_coefficient = 0.5
    cost_constant = 2.5

    def __init__(self, d: float, D: float, B: float, i: int = 1, Z: int = None, Dw: float = None,
                 Cr: float = None, C0r: float = None,
                 material: Material = material_iso,
                 contact_type_point: bool = False, contact_type_linear: bool = True, contact_type_mixed: bool = False,
                 mass: float = None,
                 oil_speed_limit: float = None,
                 grease_speed_limit: float = None,
                 name: str = '', infos: List[str] = None, speed_infos: List[float] = None):


        RadialRollerBearing.__init__(self, d, D, B, alpha=0, i=i, Z=Z, Dw=Dw, Cr=Cr,
                                     C0r=C0r,
                                     material=material,
                                     contact_type_point=contact_type_point, contact_type_linear=contact_type_linear, contact_type_mixed=contact_type_mixed,
                                     mass=mass,
                                     oil_speed_limit=oil_speed_limit,
                                     grease_speed_limit=grease_speed_limit,
                                     name=name, infos=infos, speed_infos=speed_infos)


    def internal_ring_contour(self, direction=1, sign_V=1):

        d1 = self.d1
        pbi2 = vm.Point2D(-direction*self.B/2., sign_V*d1/2.)
        pbi1 = vm.Point2D(-direction*(self.B/2. - self.h), sign_V*d1/2.)
        pbi0 = vm.Point2D(-direction*(self.B/2. - self.h), sign_V*(self.F/2.))
        pbi3 = vm.Point2D(-direction*self.B/2., sign_V*self.d/2.)
        pbi4 = vm.Point2D(direction*self.B/2., sign_V*self.d/2.)
        pbi5 = vm.Point2D(direction*self.B/2., sign_V*d1/2.)
        pbi6 = vm.Point2D(direction*(self.B/2. - self.h), sign_V*d1/2.)
        pbi7 = vm.Point2D(direction*(self.B/2. - self.h), sign_V*(self.F/2.))
        irc = primitives2d.ClosedRoundedLineSegments2D([pbi0, pbi1, pbi2, pbi3, pbi4, pbi5, pbi6, pbi7],
                                                       {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius,
                                                        5: self.radius, 6: self.radius}, adapt_radius=True)

        # cbi1 = vm.edges.Arc2D(pbi1, vm.Point2D(0, self.F/2), pbi7)

        return None, irc



    def external_ring_contour(self, direction=1, sign_V=1):

        D1 = self.E + 0.1*(self.D - self.E)
        pbe2 = vm.Point2D(-direction*self.B/2., sign_V*D1/2.)
        pbe1 = vm.Point2D(-direction*(self.B/2. - self.h), sign_V*(self.E/2.))
        pbe3 = vm.Point2D(-direction*self.B/2., sign_V*self.D/2.)
        pbe4 = vm.Point2D(direction*self.B/2., sign_V*self.D/2.)
        pbe5 = vm.Point2D(direction*self.B/2., sign_V*D1/2.)
        pbe6 = vm.Point2D(direction*(self.B/2. - self.h), sign_V*(self.E/2.))
        be1 = primitives2d.ClosedRoundedLineSegments2D([pbe1, pbe2, pbe3, pbe4, pbe5, pbe6],
                                                       {1: self.radius,
                                                        2: self.radius, 3: self.radius, 4: self.radius}, adapt_radius=True)

        # cbe1 = vm.edges.Arc2D(pbe1, vm.Point2D(0, self.F/2), pbe6)


        # erc = vm.Contour2D([be1])
        return None, be1

    @classmethod
    def graph(cls, list_node, direction=1):

        graph = nx.DiGraph()
        graph.add_edges_from([(list_node[4], list_node[0])])
        graph.add_edges_from([(list_node[1], list_node[5])])
        graph.add_edges_from([(list_node[6], list_node[2])])
        graph.add_edges_from([(list_node[3], list_node[7])])

        return graph


class NF(RadialRollerBearing):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['Dpw', 'Lw', 'h', 'E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _non_hash_attributes = ['Dpw', 'Lw', 'h', 'E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _generic_eq = True

    symmetric = True
    taking_loads = 'right'
    generate_axial_load = False
    class_name = 'NF'
    cost_coefficient = 0.5
    cost_constant = 2.5


    def __init__(self, d: float, D: float, B: float, i: int = 1, Z: int = None,
                 Dw: float = None, Cr: float = None, C0r: float = None,
                 material: Material = material_iso,
                 contact_type_point: bool = False, contact_type_linear: bool = True, contact_type_mixed: bool = False,
                 mass: float = None,
                 oil_speed_limit: float = None,
                 grease_speed_limit: float = None,
                 name: str = '', infos: List[str] = None, speed_infos: List[float] = None):


        RadialRollerBearing.__init__(self, d, D, B, alpha=0, i=i, Z=Z, Dw=Dw, Cr=Cr,
                                     C0r=C0r,
                                     material=material,
                                     contact_type_point=contact_type_point, contact_type_linear=contact_type_linear, contact_type_mixed=contact_type_mixed,
                                     mass=mass,
                                     oil_speed_limit=oil_speed_limit,
                                     grease_speed_limit=grease_speed_limit,
                                     name=name, infos=infos, speed_infos=speed_infos)


    def internal_ring_contour(self, direction=1, sign_V=1):

        d1 = self.d1
        pbi2 = vm.Point2D(-direction*self.B/2., sign_V*d1/2.)
        pbi1 = vm.Point2D(-direction*(self.B/2. - self.h), sign_V*d1/2.)
        pbi0 = vm.Point2D(-direction*(self.B/2. - self.h), sign_V*(self.F/2.))
        pbi3 = vm.Point2D(-direction*self.B/2., sign_V*self.d/2.)
        pbi4 = vm.Point2D(direction*self.B/2., sign_V*self.d/2.)
        pbi5 = vm.Point2D(direction*self.B/2., sign_V*d1/2.)
        pbi6 = vm.Point2D(direction*(self.B/2. - self.h), sign_V*d1/2.)
        pbi7 = vm.Point2D(direction*(self.B/2. - self.h), sign_V*(self.F/2.))
        irc = primitives2d.ClosedRoundedLineSegments2D([pbi0, pbi1, pbi2, pbi3, pbi4, pbi5, pbi6, pbi7],
                                                       {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius,
                                                        5: self.radius, 6: self.radius}, adapt_radius=True)

        # cbi1 = vm.edges.Arc2D(pbi1, vm.Point2D(0, self.F/2), pbi7)

        return None, irc


    def external_ring_contour(self, direction=1, sign_V=1):

        D1 = self.D1
        D2 = self.E + 0.1*(self.D - self.E)
        pbe2 = vm.Point2D(-direction*self.B/2., sign_V*D1/2.)
        pbe1 = vm.Point2D(-direction*(self.B/2. - self.h), sign_V*D1/2.)
        pbe0 = vm.Point2D(-direction*(self.B/2. - self.h), sign_V*(self.E/2.))
        pbe3 = vm.Point2D(-direction*self.B/2., sign_V*self.D/2.)
        pbe4 = vm.Point2D(direction*self.B/2., sign_V*self.D/2.)
        pbe5 = vm.Point2D(direction*self.B/2., sign_V*D2/2.)
        pbe6 = vm.Point2D(direction*(self.B/2. - self.h), sign_V*(self.E/2.))
        be1 = primitives2d.ClosedRoundedLineSegments2D([pbe0, pbe1, pbe2, pbe3, pbe4, pbe5, pbe6],
                                                       {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius,
                                                        5: self.radius}, adapt_radius=True)

        # cbe1 = vm.edges.Arc2D(pbe1, vm.Point2D(0, self.F/2), pbe6)


        # erc = vm.Contour2D([be1])

        return None, be1

    @classmethod
    def graph(cls, list_node, direction=1):

        graph = nx.DiGraph()
        graph.add_edges_from([(list_node[4], list_node[0])])
        graph.add_edges_from([(list_node[1], list_node[5])])
        graph.add_edges_from([(list_node[6], list_node[2])])
        graph.add_edges_from([(list_node[3], list_node[7])])

        if direction == 1:
            graph.add_edges_from([(list_node[4], list_node[2])])
            graph.add_edges_from([(list_node[3], list_node[5])])
        elif direction == -1:
            graph.add_edges_from([(list_node[1], list_node[7])])
            graph.add_edges_from([(list_node[6], list_node[0])])

        return graph

class NU(RadialRollerBearing):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['Dpw', 'Lw', 'h', 'E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _non_hash_attributes = ['Dpw', 'Lw', 'h', 'E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _generic_eq = True

    symmetric = True
    taking_loads = 'free'
    generate_axial_load = False
    class_name = 'NU'
    cost_coefficient = 0.5
    cost_constant = 2.5

    def __init__(self, d: float, D: float, B: float, i: int = 1, Z: int = None, Dw: float = None,
                 Cr: float = None, C0r: float = None,
                 material: Material = material_iso,
                 contact_type_point: bool = False, contact_type_linear: bool = True,
                 contact_type_mixed: bool = False,
                 mass: float = None,
                 oil_speed_limit: float = None,
                 grease_speed_limit: float = None,
                 name: str = '', infos: List[str] = None, speed_infos: List[float] = None):


        RadialRollerBearing.__init__(self, d, D, B, alpha=0, i=i, Z=Z, Dw=Dw, Cr=Cr,
                                     C0r=C0r,
                                     material=material,
                                     contact_type_point=contact_type_point, contact_type_linear=contact_type_linear, contact_type_mixed=contact_type_mixed,
                                     mass=mass,
                                     oil_speed_limit=oil_speed_limit,
                                     grease_speed_limit=grease_speed_limit,
                                     name=name, infos=infos, speed_infos=speed_infos)


    def internal_ring_contour(self, direction=1, sign_V=1):
        d1 = self.F - 0.1*(self.F - self.d)
        pbi2 = vm.Point2D(-direction*self.B/2., sign_V*d1/2.)
        pbi1 = vm.Point2D(-direction*(self.B/2. - self.h), sign_V*(self.F/2.))
        pbi3 = vm.Point2D(-direction*self.B/2., sign_V*self.d/2.)
        pbi4 = vm.Point2D(direction*self.B/2., sign_V*self.d/2.)
        pbi5 = vm.Point2D(direction*self.B/2., sign_V*d1/2.)
        pbi6 = vm.Point2D(direction*(self.B/2. - self.h), sign_V*(self.F/2.))
        irc = primitives2d.ClosedRoundedLineSegments2D([pbi1, pbi2, pbi3, pbi4, pbi5, pbi6],
                                                       {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius},
                                                       adapt_radius=True)

        # cbi1 = vm.edges.Arc2D(pbi1, vm.Point2D(0, self.F/2), pbi6)

        return None, irc


    def external_ring_contour(self, direction=1, sign_V=1):
        D1 = self.D1
        pbe2 = vm.Point2D(-direction*self.B/2., sign_V*D1/2.)
        pbe1 = vm.Point2D(-direction*(self.B/2. - self.h), sign_V*D1/2.)
        pbe0 = vm.Point2D(-direction*(self.B/2. - self.h), sign_V*(self.E/2.))
        pbe3 = vm.Point2D(-direction*self.B/2., sign_V*self.D/2.)
        pbe4 = vm.Point2D(direction*self.B/2., sign_V*self.D/2.)
        pbe5 = vm.Point2D(direction*self.B/2., sign_V*D1/2.)
        pbe6 = vm.Point2D(direction*(self.B/2. - self.h), sign_V*D1/2.)
        pbe7 = vm.Point2D(direction*(self.B/2. - self.h), sign_V*(self.E/2.))
        be1 = primitives2d.ClosedRoundedLineSegments2D([pbe0, pbe1, pbe2, pbe3, pbe4, pbe5, pbe6, pbe7],
                                                       {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius,
                                                        5: self.radius, 6: self.radius}, adapt_radius=True)
        # erc = vm.Contour2D([be1])
        # cbe1 = vm.edges.Arc2D(pbe1, vm.Point2D(0, self.F/2), pbe7)

        return None, be1

    @classmethod
    def graph(cls, list_node, direction=1):

        graph = nx.DiGraph()
        graph.add_edges_from([(list_node[4], list_node[0])])
        graph.add_edges_from([(list_node[1], list_node[5])])
        graph.add_edges_from([(list_node[6], list_node[2])])
        graph.add_edges_from([(list_node[3], list_node[7])])

        return graph

class TaperedRollerBearing(RadialRollerBearing, AngularBallBearing):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['beta', 'Dpw', 'Lw', 'h', 'E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _non_hash_attributes = ['beta', 'Dpw', 'Lw', 'h', 'E', 'F', 'd1', 'D1', 'radius', 'slack', 'mass', 'cost', 'name']
    _generic_eq = True

    symmetric = False
    taking_loads = 'right'
    linkage = 'cylindric_joint'
    generate_axial_load = True
    coeff_baselife = 10/3.
    class_name = 'TaperedRollerBearing'
    cost_coefficient = 0.4
    cost_constant = 2


    def __init__(self, d: float, D: float, B: float, alpha: float, i: int = 1, Z: int = None,
                 Dw: float = None, Cr: float = None, C0r: float = None,
                 material: Material = material_iso,
                 contact_type_point: bool = False, contact_type_linear: bool = True, contact_type_mixed: bool = False,
                 mass: float = None,
                 oil_speed_limit: float = None,
                 grease_speed_limit: float = None,
                 name: str = '', infos: List[str] = None, speed_infos: List[float] = None):


        if Dw is None:
            self.Dw = (D - d)/7.*math.cos(alpha)

        RadialRollerBearing.__init__(self, d, D, B, alpha=alpha, i=i, Z=Z,
                                     Dw=Dw, Cr=Cr, C0r=C0r,
                                     material=material,
                                     contact_type_point=contact_type_point, contact_type_linear=contact_type_linear, contact_type_mixed=contact_type_mixed,
                                     mass=mass,
                                     oil_speed_limit=oil_speed_limit,
                                     grease_speed_limit=grease_speed_limit,
                                     name=name, infos=infos, speed_infos=speed_infos)

        # estimation for the graph 2D description
        self.Dpw = (self.d + self.D)/2.
        self.Lw = 0.7*self.B
        self.beta = math.atan(self.Dw/self.Dpw*math.sin(self.alpha))

    def internal_ring_contour(self, direction=1, sign_V=1):

        shift_bi = 5e-4
#        shift_be = 1e-3

        p0 = vm.Point2D(0, sign_V*self.Dpw/2.)
        p0v = p0.to_vector()
        p1v = p0v.translation((math.cos(self.alpha), -direction*sign_V*math.sin(self.alpha)), True)
        p1 = p1v.to_point()
        l1 = vm.edges.Line2D(p0, p1)
        l1.rotation(p0, direction*sign_V*self.beta)
        l1.translation((-0.8*direction*self.Dw/2.*math.sin(self.alpha), -sign_V*0.8*self.Dw/2.*math.cos(self.alpha)))
        l2 = l1.translation((-0.2*direction*self.Dw/2.*math.sin(self.alpha), -sign_V*0.2*self.Dw/2.*math.cos(self.alpha)), True)
        pbi3 = vm.Point2D(direction*(self.B/2. - shift_bi), sign_V*self.d/2.)
        pbi3v = pbi3.to_vector()
        pbi3Tv = pbi3v.translation((0, 1))
        pbi3T = pbi3Tv.to_point()
        pbi4 = vm.Point2D(-direction*(self.B/2.), sign_V*self.d/2.)
        pbi4v = pbi4.to_vector()
        pbi4Tv = pbi4v.translation((0, 1))
        pbi4T = pbi4Tv.to_point()
        l3 = vm.edges.Line2D(pbi3, pbi3T)
        l4 = vm.edges.Line2D(pbi4, pbi4T)
        pbi2 = vm.Point2D.line_intersection(l1, l3)
        pbi5 = vm.Point2D.line_intersection(l1, l4)
        l5 = vm.edges.Line2D(vm.Point2D(direction*self.Lw/2., 0), vm.Point2D(direction*self.Lw/2., 1))
        l5.rotation(vm.Point2D(0, sign_V*self.Dpw/2.), -sign_V*direction*self.alpha)
        l6 = vm.edges.Line2D(vm.Point2D(-direction*self.Lw/2., 0), vm.Point2D(-direction*self.Lw/2., 1))
        l6.rotation(vm.Point2D(0, sign_V*self.Dpw/2.), -sign_V*direction*self.alpha)
        pbi1 = vm.Point2D.line_intersection(l1, l5)
        pbi0 = vm.Point2D.line_intersection(l2, l5)
        pbi6 = vm.Point2D.line_intersection(l1, l6)
        pbi7 = vm.Point2D.line_intersection(l2, l6)

        bi1 = primitives2d.ClosedRoundedLineSegments2D([pbi0, pbi1, pbi2, pbi3, pbi4, pbi5, pbi6, pbi7],
                                                       {1: self.radius, 2: self.radius, 3: self.radius,
                                                        4: self.radius, 5: self.radius,
                                                        6: self.radius}, adapt_radius=False)

        # irc = vm.Contour2D([bi1])
#
        return bi1

    def external_ring_contour(self, direction=1, sign_V=1):

#        shift_bi = 5e-4
        shift_be = 1e-3

        p0 = vm.Point2D(0, sign_V*self.Dpw/2.)
        p1 = p0.to_vector().translation(vm.Vector2D(math.cos(self.alpha), -direction*sign_V*math.sin(self.alpha)), True).to_point()
        l0 = vm.edges.Line2D(p0, p1)
        l0.rotation(p0, -direction*sign_V*self.beta)
        l0.translation(vm.Vector2D(direction*self.Dw/2.*math.sin(self.alpha), sign_V*self.Dw/2.*math.cos(self.alpha)))
        pbe3 = vm.Point2D(direction*self.B/2., sign_V*self.D/2.)
        pbe3T = pbe3.to_vector().translation(vm.Vector2D(0, 1)).to_point()
        pbe4 = vm.Point2D(-direction*(self.B/2. - shift_be), sign_V*self.D/2.)
        pbe4T = pbe4.to_vector().translation(vm.Vector2D(0, 1)).to_point()
        l3 = vm.edges.Line2D(pbe3, pbe3T)
        l4 = vm.edges.Line2D(pbe4, pbe4T)
        pbe2 = vm.Point2D.line_intersection(l0, l3)
        pbe5 = vm.Point2D.line_intersection(l0, l4)
        be1 = primitives2d.ClosedRoundedLineSegments2D([pbe2, pbe3, pbe4, pbe5],
                                                       {0: self.radius,
                                                        1: self.radius,
                                                        2: self.radius,
                                                        3: self.radius},
                                                       adapt_radius=True)

        return be1

    def volmdlr_primitives(self, center=vm.O3D, axis=vm.X3D):
        axis.normalize()

        y = axis.random_unit_normal_vector()
#        y.vector = npy.round(y.vector,3)
#        y.vector = y.vector/y.Norm()

#        z=vm.Vector3D(npy.cross(x.vector,y.vector))
#        z = axis.Cross(y)

        #Internal Ring
        IRC = vm.wires.Contour2D(self.internal_ring_contour().primitives)

        irc = primitives3d.RevolvedProfile(center, axis, y, IRC, center,
                                           axis, angle=2*math.pi, name='Internal Ring')
        #External Ring

        ERC = vm.wires.Contour2D(self.external_ring_contour().primitives)
        erc = primitives3d.RevolvedProfile(center, axis, y, ERC, center,
                                           axis, angle=2*math.pi, name='External Ring')

        volumes = [irc, erc]
        return volumes



    def rolling_contour(self, direction=1, sign_V=1):

        r1 = vm.Point2D(--direction*self.Lw/2., self.Dw/2. - self.Lw/2.*math.tan(self.beta))
        r2 = vm.Point2D(-direction*self.Lw/2., self.Dw/2. + self.Lw/2.*math.tan(self.beta))
        r3 = vm.Point2D(-direction*self.Lw/2., -self.Dw/2. - self.Lw/2.*math.tan(self.beta))
        r4 = vm.Point2D(--direction*self.Lw/2., -self.Dw/2. + self.Lw/2.*math.tan(self.beta))
        rol = primitives2d.ClosedRoundedLineSegments2D([r1, r2, r3, r4],
                                                       {0: self.radius, 1: self.radius, 2: self.radius, 3: self.radius})

#        bg = vm.Contour2D([rol])
        return vm.wires.Contour2D(rol.primitives)

    def plot_data(self, pos=0, pos_y=0, direction=1, quote=True, constructor=True, stroke_width=1):

        plot_datas = []
        be_sup = self.external_ring_contour(direction=direction, sign_V=1)
        be_sup1 = be_sup.translation((pos, pos_y), True)
        plot_datas.append(be_sup1.plot_data())

        be_inf = self.external_ring_contour(direction=direction, sign_V=-1)
        be_inf1 = be_inf.translation((pos, pos_y), True)
        plot_datas.append(be_inf1.plot_data())

        bi_sup = self.internal_ring_contour(direction=direction, sign_V=1)
        bi_sup1 = bi_sup.translation((pos, pos_y), True)
        plot_datas.append(bi_sup1.plot_data())

        bi_inf = self.internal_ring_contour(direction=direction, sign_V=-1)
        bi_inf1 = bi_inf.translation((pos, pos_y), True)
        plot_datas.append(bi_inf1.plot_data())
#
        roller_sup = self.rolling_contour(direction=direction, sign_V=1)
        roller_sup = roller_sup.rotation(vm.Point2D(0, 0), -direction*self.alpha, True)
        roller_sup = roller_sup.translation((0, self.Dpw/2.), True)
        roller_sup1 = roller_sup.translation((pos, pos_y), True)
        plot_datas.append(roller_sup1.plot_data())

        roller_inf = self.rolling_contour(direction=direction, sign_V=-1)
        roller_inf = roller_inf.rotation(vm.Point2D(0, 0), direction*self.alpha, True)
        roller_inf = roller_inf.translation((0, -self.Dpw/2.), True)
        roller_inf1 = roller_inf.translation((pos, pos_y), True)
        plot_datas.append(roller_inf1.plot_data())

        return [plot_data.PrimitiveGroup(plot_datas)]

    def plot_contour(self, direction=1):

        be_sup = self.external_ring_contour(direction=direction, sign_V=1)
        be_inf = self.external_ring_contour(direction=direction, sign_V=-1)
        bi_sup = self.internal_ring_contour(direction=direction, sign_V=1)
        bi_inf = self.internal_ring_contour(direction=direction, sign_V=-1)
        roller_sup = self.rolling_contour(direction=direction, sign_V=1)
        roller_sup = roller_sup.rotation(vm.Point2D(0, 0), -direction*self.alpha, True)
        roller_sup = roller_sup.translation(vm.Vector2D(0, self.Dpw/2.), True)
        roller_inf = self.rolling_contour(direction=direction, sign_V=-1)
        roller_inf = roller_inf.rotation(vm.Point2D(0, 0), direction*self.alpha, True)
        roller_inf = roller_inf.translation(vm.Vector2D(0, -self.Dpw/2.), True)

        bg = vm.wires.Contour2D([be_sup, bi_sup, roller_sup, be_inf, bi_inf, roller_inf])
        return bg

#    @classmethod
#    def dict_to_object(cls, d):
#        if 'Cr' not in d.keys():
#            d['Cr'] = None
#        if 'C0r' not in d.keys():
#            d['C0r'] = None
#        obj = cls(d = d['d'], D = d['D'], B = d['B'], alpha = d['alpha'], i = d['i'], Z = d['Z'],
#                  Dw = d['Dw'], Cr = d['Cr'], C0r = d['C0r'],
#                  material = Material.dict_to_object(d['material']),
#                  contact_type = d['contact_type'],
#                  name=d['name'], mass=d['mass'])
#        return obj

    def Copy(self):
        if not hasattr(self, 'Cr'):
            Cr = None
        else:
            Cr = self.Cr
        if not hasattr(self, 'C0r'):
            C0r = None
        else:
            C0r = self.C0r
        obj = TaperedRollerBearing(d=self.d, D=self.D, B=self.B, alpha=self.alpha,
                                   i=self.i, Z=self.Z,
                                   Dw=self.Dw, Cr=Cr, C0r=C0r,
                                   oil=self.oil, material=self.material,
                                   contact_type_point=self.contact_type_point, contact_type_linear=self.contact_type_linear, contact_type_mixed=self.contact_type_mixed,
                                   name=self.name)
        return obj

class BearingCatalog(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = ['bearings_by_types', 'bearings_by_dict']
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, bearings: List[RadialBearing], name: str = ''):
        self.bearings = bearings
        self.bearings_by_types = {}
        for bearing in bearings:
            if bearing.__class__ in self.bearings_by_types:
                self.bearings_by_types[bearing.__class__].append(bearing)
            else:
                self.bearings_by_types[bearing.__class__] = [bearing]

        bearings_by_types = {}
        for bearing in bearings:
            d = str(round(bearing.d, 3))
            D = str(round(bearing.D, 3))

            if bearing.__class__ in bearings_by_types:
                if d not in bearings_by_types.keys():
                    bearings_by_types[bearing.__class__][d] = {}
                    bearings_by_types[bearing.__class__][d][D] = [bearing]
                else:
                    if D not in bearings_by_types[d].keys():
                        bearings_by_types[bearing.__class__][d][D] = []
                    else:
                        bearings_by_types[bearing.__class__][d][D].append(bearing)
            else:
                bearings_by_types[bearing.__class__] = {}
                bearings_by_types[bearing.__class__][d] = {}
                bearings_by_types[bearing.__class__][d][D] = [bearing]

        self.bearings_by_dict = {}
        for classe, dict_classe in bearings_by_types.items():
            self.bearings_by_dict[classe] = {}
            for d, dict_d in dict_classe.items():
                self.bearings_by_dict[classe][d] = {}
                for D, dict_D in dict_d.items():
                    list_Cr = []
                    for bearing in dict_D:
                        list_Cr.append(bearing.Cr)
                    self.bearings_by_dict[classe][d][D] = [dict_D[i] for i in npy.argsort(list_Cr)]

        DessiaObject.__init__(self, name=name)

    def find_duplicates(self):
        duplicates = []
        for bearing in self.bearings:
            if self.bearings.count(bearing) > 1:
                duplicates.append(bearing)
        return duplicates



    @classmethod
    def load_from_dataframe(cls, dataframe, catalog_name):
        bearings = []
        for index in dataframe.index:
            typ_rlt = dataframe.loc[index, 'typ_bearing']
            d = round(dataframe.loc[index, 'd'], 3)
            D = round(dataframe.loc[index, 'D'], 3)
            B = round(dataframe.loc[index, 'B'], 3)
            i = dataframe.loc[index, 'i']
            Z = dataframe.loc[index, 'Z']
            Dw = dataframe.loc[index, 'Dw']
            mass = None
            alp = dataframe.loc[index, 'alpha']
            if str(alp) == 'nan':
                alpha = 0
            else:
                alpha = alp
            Cr = dataframe.loc[index, 'Cr']*1e3
            C0r = dataframe.loc[index, 'C0r']

            if typ_rlt == 'radial_roller_bearing':
                typ = dataframe.loc[index, 'mounting']
                if typ == 'NUP':
                    bearings.append(NUP(d, D, B, i,
                                        Z, Dw, Cr, C0r,
                                        mass=mass))
                elif typ == 'N':
                    bearings.append(N(d, D, B, i,
                                      Z, Dw, Cr, C0r,
                                      mass=mass))
                elif typ == 'NF':
                    bearings.append(NF(d, D, B, i,
                                       Z, Dw, Cr, C0r,
                                       mass=mass,))
                elif typ == 'NU':
                    bearings.append(NU(d, D, B, i,
                                       Z, Dw, Cr, C0r,
                                       mass=mass,))

            elif typ_rlt == 'radial_ball_bearing':
                bearings.append(RadialBallBearing(d, D, B, i, Z,
                                                  Dw, Cr, C0r, mass=mass))

            elif typ_rlt == 'angular_ball_bearing':
                bearings.append(AngularBallBearing(d, D, B, i, Z, Dw, alpha,
                                                   Cr, C0r, mass=mass))

            elif typ_rlt == 'spherical_ball_bearing':
                bearings.append(SphericalBallBearing(d, D, B, i, Z, Dw, alpha,
                                                     Cr, C0r, mass=mass))

            elif typ_rlt == 'tapered_roller_bearing':
                bearings.append(TaperedRollerBearing(d, D, B, i, Z, Dw,
                                                     Cr, C0r, mass=mass))
        return cls(bearings, catalog_name)

    def search_bearing_catalog(self, bearing_class, d, D):

        if bearing_class in self.bearings_by_types:

            bearings = self.bearings_by_types[bearing_class]
            list_bearings = []
            list_sort = []
            for bearing in bearings:

                if (bearing.d >= d) and (bearing.D <= D):
                    list_bearings.append(bearing)
                    list_sort.append(bearing.Cr)
            arg_list_sort = npy.argsort(list_sort)
            return [list_bearings[i] for i in arg_list_sort]
        else:
            return []

    def next_bearing_catalog(self, bearing_class, d, D):
        try:
            next_bearings = self.bearings_by_dict[bearing_class][str(d)][str(D)]
            return next_bearings
        except KeyError:
            return False

    def check(self):
        for bearing_class, bearings in self.bearings.items():
            for bearing in bearings:
                if not bearing.Check():
                        return False
        return True

    def invalid_bearings(self):
        invalid_bearings = []
        for bearing in self.bearings:
            if not bearing.Check():
                invalid_bearings.append(bearing)
        return invalid_bearings

    def plot(self):
        d = [b.d for b in self.bearings]
        D = [b.D for b in self.bearings]
        B = [b.B for b in self.bearings]
        Cr = [b.Cr for b in self.bearings]
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
        ax1.plot(d, D, 'o')
        ax1.set_xlabel('d')
        ax1.set_ylabel('D')

        ax2.plot(d, B, 'o')
        ax2.set_xlabel('d')
        ax2.set_ylabel('B')

        ax3.plot(d, Cr, 'o')
        ax3.set_xlabel('d')
        ax3.set_ylabel('Cr')

        ax4.plot(D, Cr, 'o')
        ax4.set_xlabel('D')
        ax4.set_ylabel('Cr')


    def _display_angular(self):

        filters = [{'attribute' : 'd'},
                   {'attribute' : 'D'},
                   {'attribute' : 'B'},
                   {'attribute' : 'Cr'},
                   {'attribute' : 'C0r'}]

        bearings_index = {b:i for i, b in enumerate(self.bearings)}

        datasets_values = []
        datasets_names = []

        values = []
        for bearing in self.bearings:
            value = {'d': bearing.d,
                     'D': bearing.D,
                     'B': bearing.B,
                     'Cr': bearing.Cr,
                     'C0r': bearing.C0r}
            values.append(value)

        for class_, bearings in self.bearings_by_types.items():
            datasets_names.append(class_.__name__)
            dataset_values_indices = []
            for bearing in bearings:
                dataset_values_indices.append(bearings_index[bearing])
            datasets_values.append(dataset_values_indices)
        nds = len(datasets_values)
        datasets = []
        for ids, (dataset_values, name) in enumerate(zip(datasets_values, datasets_names)):
            datasets.append({'label' : name,
                             'color' : matplotlib.colors.to_hex(matplotlib.colors.hsv_to_rgb((ids/(nds-1), 0.8, 0.7))),
                             'values' : dataset_values,
                             })

        displays = [{'angular_component': 'results',
                     'filters': filters,
                     'references_attribute': 'bearings',
                     'values': values,
                     'datasets': datasets}]
        return displays

bearing_classes_ = [RadialBallBearing, AngularBallBearing,
                    NUP, N,
                    NF, NU,
                    TaperedRollerBearing]

dict_bearing_classes = {str(RadialBallBearing): RadialBallBearing,
                        str(AngularBallBearing): AngularBallBearing,
                        str(NUP): NUP,
                        str(N): N,
                        str(NF): NF,
                        str(NU): NU,
                        str(TaperedRollerBearing): TaperedRollerBearing}

strength_bearing_classes = {str(RadialBallBearing): 1,
                            str(AngularBallBearing): 1,
                            str(NUP): 3,
                            str(N): 3,
                            str(NF): 3,
                            str(NU): 3,
                            str(TaperedRollerBearing): 2}

class ConceptualBearingCombination(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = ['bearing_classes']
    _generic_eq = True

    def __init__(self, bearing_classes: List[RadialBearing], directions: List[int], mounting: Mounting,
                 name: str = ''):

        self.bearing_classes = bearing_classes
        self.directions = directions
        self.mounting = mounting

        DessiaObject.__init__(self, name=name)

#    def __eq__(self, other_eb):
#
#        equal = True
#        for bg, other_bg in zip(self.bearing_classes, other_eb.bearing_classes):
#            equal = equal and bg == other_bg
#        equal = equal and self.directions == other_eb.directions
#        equal = equal and self.mounting == other_eb.mounting
#        return equal
#
#    def __hash__(self):
#        h = 3*len(self.mounting) + 127*sum(self.directions)
#        for class_ in self.bearing_classes:
#            h += len(class_.__name__)
#        return h

    def bearing_combination(self, bearings):

        if self.mounting.left and self.mounting.right:
            connection_bi = Mounting(left=True, right=True)
            connection_be = Mounting(left=True, right=True)


        elif not self.mounting.left and not self.mounting.right:
            connection_bi = Mounting(left=True, right=True)
            connection_be = Mounting(left=False, right=False)
        elif self.mounting.left:
            connection_bi = Mounting(left=False, right=True)
            connection_be = Mounting(left=True, right=False)
        elif self.mounting.right:
            connection_bi = Mounting(left=True, right=False)
            connection_be = Mounting(left=False, right=True)
#        if self.mounting == 'both':
#            connection_bi=['left', 'right']
#            connection_be=['left', 'right']
#        elif self.mounting == 'free':
#            connection_bi=['left', 'right']
#            connection_be=[]
#        elif self.mounting == 'right':
#            connection_bi=['left']
#            connection_be=['right']
#        elif self.mounting == 'left':
#            connection_bi=['right']
#            connection_be=['left']
        return BearingCombination(bearings, self.directions, radial_load_linkage=[True]*len(bearings),
                                  connection_bi=connection_bi, connection_be=connection_be,
                                  behavior_link=self.mounting)

    def check_kinematic(self):

        list_node_output = list(range(8))
        list_node_bearings = []
        for i, bearing_classe in enumerate(self.bearing_classes):
            list_node_bearings.append(npy.arange(8*(i + 1), 8*(i + 2)))

        nx_graph = self.graph(list_node_bearings, list_node_output)

        check = self.check_viability_angular_bearing(nx_graph, list_node_bearings,
                                                     list_node_output)
        return check

    def graph(self, list_node_bearings, list_node_output):

        list_left_bearing = list_node_bearings[0]
        list_right_bearing = list_node_bearings[-1]
        nx_graph = self.bearing_classes[0].graph(list_left_bearing, self.directions[0])

        bg0_load = self.bearing_classes[0].taking_loads

        if self.directions[0] == -1:
            if self.bearing_classes[0].taking_loads == 'left':
                bg0_load = 'right'
            elif self.bearing_classes[0].taking_loads == 'right':
                bg0_load = 'left'
        # if (self.mounting.left) or (self.mounting == 'both'):
        if self.mounting.left:
            if bg0_load in ['left', 'both']:
                nx_graph.add_edges_from([(list_node_output[1], list_left_bearing[1])])
                nx_graph.add_edges_from([(list_left_bearing[0], list_node_output[0])])
        # if (self.mounting == 'right') or (self.mounting == 'both'):
        if self.mounting.right:
            if bg0_load in ['right', 'both']:
                nx_graph.add_edges_from([(list_node_output[3], list_left_bearing[3])])
                nx_graph.add_edges_from([(list_left_bearing[2], list_node_output[2])])

        for num_bg, (dir1, dir2, bg1, bg2) in enumerate(zip(self.directions[0: -1], self.directions[1:], \
                    self.bearing_classes[0: -1], self.bearing_classes[1:])):
            list_nd1 = list_node_bearings[num_bg]
            list_nd2 = list_node_bearings[num_bg + 1]
            nx_graph = nx.compose(bg2.graph(list_nd2, dir2), nx_graph)

            bg1_load = bg1.taking_loads
            bg2_load = bg2.taking_loads
            if dir1 == -1:
                if bg1.taking_loads == 'left':
                    bg1_load = 'right'
                elif bg1.taking_loads == 'right':
                    bg1_load = 'left'
            if dir2 == -1:
                if bg2.taking_loads == 'left':
                    bg2_load = 'right'
                elif bg2.taking_loads == 'right':
                    bg2_load = 'left'

            if (bg1_load == 'right' and bg2_load != 'right') or \
                (bg1_load != 'left' and bg2_load == 'left'):
                nx_graph.add_edges_from([(list_nd2[0], list_nd1[4])])
                nx_graph.add_edges_from([(list_nd1[5], list_nd2[1])])
            elif (bg1_load == 'left' and bg2_load != 'left') or \
                (bg1_load != 'right' and bg2_load == 'right'):
                nx_graph.add_edges_from([(list_nd2[2], list_nd1[6])])
                nx_graph.add_edges_from([(list_nd1[7], list_nd2[3])])
            else:
                nx_graph.add_edges_from([(list_nd2[0], list_nd1[4])])
                nx_graph.add_edges_from([(list_nd1[5], list_nd2[1])])
                nx_graph.add_edges_from([(list_nd2[2], list_nd1[6])])
                nx_graph.add_edges_from([(list_nd1[7], list_nd2[3])])

        bg_end_load = self.bearing_classes[-1].taking_loads

        if self.directions[-1] == -1:
            if self.bearing_classes[-1].taking_loads == 'left':
                bg_end_load = 'right'
            elif self.bearing_classes[-1].taking_loads == 'right':
                bg_end_load = 'left'
        # if (self.mounting == 'left') or (self.mounting == 'both'):
        if self.mounting.left:
            if bg_end_load in ['left', 'both']:
                nx_graph.add_edges_from([(list_node_output[6], list_right_bearing[6])])
                nx_graph.add_edges_from([(list_right_bearing[7], list_node_output[7])])
        # if (self.mounting == 'right') or (self.mounting == 'both'):
        if self.mounting.right:
            if bg_end_load in ['right', 'both']:
                nx_graph.add_edges_from([(list_node_output[4], list_right_bearing[4])])
                nx_graph.add_edges_from([(list_right_bearing[5], list_node_output[5])])

        return nx_graph

    def check_viability_angular_bearing(self, nx_graph, list_node_bearings, li_node_output):
        # axial load generate by angular_bearing
        node_axial_ring = []
        # if self.mounting in ['left', 'both']:

        if self.mounting.left:
            node_axial_ring.append(li_node_output[0])
            node_axial_ring.append(li_node_output[7])
        # if self.mounting in ['right', 'both']:
        if self.mounting.right:
            node_axial_ring.append(li_node_output[5])
            node_axial_ring.append(li_node_output[2])

        node_axial_input = []
        for i, bearing_classe in enumerate(self.bearing_classes):
            if bearing_classe.generate_axial_load:
                if ((self.directions[i] == 1) and bearing_classe.taking_loads == 'right') or \
                    ((self.directions[i] == -1) and bearing_classe.taking_loads == 'left'):
                    node_axial_input.extend([list_node_bearings[i][j] for j in [3, 4]])
                else:
                    node_axial_input.extend([list_node_bearings[i][j] for j in [1, 6]])

        #check angular and tapered bearing
        valid = True
        for nd_axial_input in node_axial_input:
            valid_input = False
            for nd_axial_ring in node_axial_ring:
                try:
                    list(nx.all_shortest_paths(nx_graph, source=nd_axial_input,
                                               target=nd_axial_ring))
                    valid_input = True
                except nx.NetworkXNoPath:
                    pass
            if valid_input == False:
                valid = False

        #check axial mounting load
        # if self.mounting in ['right', 'both']:
        if self.mounting.right:
            node_input = li_node_output[3]
            node_output = li_node_output[5]

            try:
                if (node_input in nx_graph) and (node_output in nx_graph):
                    list(nx.all_shortest_paths(nx_graph, source=node_input,
                                               target=node_output))
                else:

                    valid = False
            except nx.NetworkXNoPath:
                valid = False
                pass
        # if self.mounting in ['left', 'both']:
        if self.mounting.left:
            node_input = li_node_output[6]
            node_output = li_node_output[0]

            try:
                if (node_input in nx_graph) and (node_output in nx_graph):
                    list(nx.all_shortest_paths(nx_graph, source=node_input,
                                               target=node_output))
                else:
                    valid = False
            except nx.NetworkXNoPath:
                valid = False
                pass


        def AnalyseConnection(node_input, node_output):
            valid = True
            try:
                if (node_input in nx_graph) and (node_output in nx_graph):
                    list(nx.all_shortest_paths(nx_graph, source=node_input,
                                               target=node_output))
                else:
                    valid = False
            except nx.NetworkXNoPath:
                valid = False
                pass
            return valid

#        if self.mounting == 'free':
#            valid_free1 = AnalyseConnection(list_node_bearings[0][3], list_node_bearings[-1][5])
#            valid_free2 = AnalyseConnection(list_node_bearings[-1][6], list_node_bearings[0][0])
#            if not valid_free1 and not valid_free2:
#                valid = True
#            else:
#                valid = False

        # if self.mounting == 'right':
        if self.mounting.right:
            if valid:
                valid = AnalyseConnection(list_node_bearings[0][4], li_node_output[2])
            if valid:
                valid = AnalyseConnection(list_node_bearings[-1][3], li_node_output[5])

        # if self.mounting == 'left':
        if self.mounting.left:
            if valid:
                valid = AnalyseConnection(list_node_bearings[0][6], li_node_output[0])
            if valid:
                valid = AnalyseConnection(list_node_bearings[-1][1], li_node_output[7])

        return valid

class BearingCombination(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, bearings: List[RadialBearing], directions: List[int],
                 radial_load_linkage: List[float], internal_pre_load: bool=False,
                 connection_bi: Mounting = None,
                 connection_be: Mounting = None,
                 behavior_link: Mounting = None,
                 axial_positions: List[float] = None,
                 internal_diameters: List[float] = None,
                 external_diameters: List[float] = None,
                 length: float = None,
                 name: str = ''):
        self.bearings = bearings
        self.radial_load_linkage = radial_load_linkage
        self.internal_pre_load = internal_pre_load
        self.connection_be = connection_be
        self.connection_bi = connection_bi
        if connection_bi is None:
            self.connection_bi = Mounting(left=True, right=True)
        else:
            self.connection_bi = connection_bi
        if connection_be is None:
            self.connection_be = Mounting(left=True, right=True)
        else:
            self.connection_be = connection_be
        if behavior_link is None:
            self.behavior_link = Mounting(left=True, right=True)
        else:
            self.behavior_link = behavior_link
        self.mass = 0
        self.cost = 0
        self.directions = directions
        self.axial_positions = axial_positions
        self.internal_diameters = internal_diameters
        self.external_diameters = external_diameters
        self.length = length

        for bg in bearings:
            if bg.mass is not None:
                self.mass += bg.mass
            if bg.cost is not None:
                self.cost += bg.cost
        self.B = 0
        for bg in bearings:
            self.B += bg.B
        self.D = 0
        for bg in bearings:
            self.D = max(self.D, bg.D)
        self.d = math.inf
        for bg in bearings:
            self.d = min(self.d, bg.d)
        self.number_bearing = len(self.bearings)

        DessiaObject.__init__(self, name=name)


#    def __eq__(self, other_eb):
#
#        equal = True
#        for bg, other_bg in zip(self.bearings, other_eb.bearings):
#            equal = equal and bg == other_bg
#        equal = (equal and self.directions == other_eb.directions
#                       and self.radial_load_linkage == other_eb.radial_load_linkage
#                       and self.internal_pre_load == other_eb.internal_pre_load
#                       and self.connection_be == other_eb.connection_be
#                       and self.connection_bi == other_eb.connection_bi
#                       and self.behavior_link == other_eb.behavior_link)
#        if hasattr(self, 'axial_positions') and hasattr(other_eb, 'axial_positions'):
#            equal = (equal and self.axial_positions == other_eb.axial_positions
#                           and self.internal_diameters == other_eb.internal_diameters
#                           and self.external_diameters == other_eb.external_diameters
#                           and self.length == other_eb.length)
#        elif hasattr(self, 'axial_positions') or hasattr(other_eb, 'axial_positions'):
#            equal = False
#        return equal
#
#    def __hash__(self):
#        h = 0
#        for bg in self.bearings:
#            h += hash(bg)
#        h += sum(self.directions)
#        return h

    # def plot_graph(self):

    #     for gp in self.graph:
    #         list_graph = []
    #         pos_m = -self.B/2.
    #         for bg in gp:
    #             graph = bg.PlotGraph()
    #             graph = graph.translation((pos_m + bg.B/2., 0), True)
    #             pos_m += bg.B
    #             list_graph.append(graph)
    #         list_graph = vm.wires.Contour2D(list_graph)
    # #        list_graph = list_graph.Translation((pos, 0), True)
    #         f,a = list_graph.MPLPlot(color='b',arrow= True)

    def update(self, axial_positions, internal_diameters, external_diameters, length):
        # TODO Why axial position is not in init?
        self.axial_positions = axial_positions# TODO: move this in bearing assembly, for plots pass a center parameter
        self.internal_diameters = internal_diameters
        self.external_diameters = external_diameters
        self.length = length

    def external_bearing(self, sign=1):
        B = self.B
        d = self.d
        D = self.D
        ep = min(0.1*D, 0.1*d)
        De = D + ep
        Dg, Dd = D, D
        if self.connection_be.left:
            Dg = Dg - ep
        if self.connection_be.right:
            Dd = Dd -ep
        be = primitives2d.ClosedRoundedLineSegments2D([vm.Point2D(-B/2., sign*D/2.), vm.Point2D(-B/2., sign*Dg/2.),
                                                       vm.Point2D(-B/2. - ep, sign*Dg/2.), vm.Point2D(-B/2. - ep, sign*De/2.),
                                                       vm.Point2D(B/2. + ep, sign*De/2.), vm.Point2D(B/2. + ep, sign*Dd/2.),
                                                       vm.Point2D(B/2., sign*Dd/2.), vm.Point2D(B/2., sign*D/2.),
                                                       vm.Point2D(-B/2., sign*D/2.)], {})
        return vm.wires.Contour2D(be.primitives)

    def internal_bearing(self, sign=1):
        B = self.B
        d = self.d
        D = self.D
        ep = min(0.1*D, 0.1*d)
        di = d - ep
        dg, dd = d, d
        if self.connection_bi.left:
            dg = dg + ep
        if self.connection_bi.right:
            dd = dd + ep
        bi = primitives2d.ClosedRoundedLineSegments2D([vm.Point2D(-B/2., sign*d/2.), vm.Point2D(-B/2., sign*dg/2.),
                                                       vm.Point2D(-B/2. - ep, sign*dg/2.), vm.Point2D(-B/2. - ep, sign*di/2.),
                                                       vm.Point2D(B/2. + ep, sign*di/2.), vm.Point2D(B/2. + ep, sign*dd/2.),
                                                       vm.Point2D(B/2., sign*dd/2.), vm.Point2D(B/2., sign*d/2.),
                                                       vm.Point2D(-B/2., sign*d/2.)], {})
        return vm.wires.Contour2D(bi.primitives)

    def bearing_box(self, sign=1):

        box = primitives2d.ClosedRoundedLineSegments2D([vm.Point2D(self.axial_positions, sign*self.internal_diameters/2.),
                                                        vm.Point2D(self.axial_positions, sign*self.external_diameters/2.),
                                                        vm.Point2D(self.axial_positions + self.length, sign*self.external_diameters/2.),
                                                        vm.Point2D(self.axial_positions + self.length, sign*self.internal_diameters/2.),
                                                        vm.Point2D(self.axial_positions, sign*self.internal_diameters/2.)], {})
        return vm.wires.Contour2D(box.primitives)

    def plot_data(self, pos=0, pos_y=0, box=False, typ=None, bearing_combination_result=None, quote=False, constructor=True):

        export_data = []

        if box:
            box_sup = self.bearing_box(1).translation(vm.Vector2D(pos, pos_y), True)
            export_data.append(box_sup.plot_data())
            box_inf = self.bearing_box(-1).translation(vm.Vector2D(pos, pos_y), True)
            export_data.append(box_inf.plot_data())


        be_sup = self.external_bearing(sign=1).translation(vm.Vector2D(pos, pos_y), True)
        export_data.append(be_sup.plot_data())
        be_inf = self.external_bearing(sign=-1).translation(vm.Vector2D(pos, pos_y), True)
        export_data.append(be_inf.plot_data())
        bi_sup = self.internal_bearing(sign=1).translation(vm.Vector2D(pos, pos_y), True)
        export_data.append(bi_sup.plot_data())
        bi_inf = self.internal_bearing(sign=-1).translation(vm.Vector2D(pos, pos_y), True)
        export_data.append(bi_inf.plot_data())


        pos_m = -self.B/2.

        for bg, di in zip(self.bearings, self.directions):
            try:
                cont = bg.plot_data(pos=pos_m + bg.B/2. + pos, pos_y=pos_y, direction=di)
            except TypeError:
                cont = bg.plot_data(pos=pos_m + bg.B/2. + pos, pos_y=pos_y)
            pos_m += bg.B
            export_data.extend(cont[0].primitives)

        if typ == 'Load':
            pos_m = -self.B/2.
            for bg_ref, bg_simu in zip(self.bearings, bearing_combination_result):
                export_data.extend(bg_simu.PlotDataLoad(pos=pos + pos_m + bg_ref.B/2., d=bg_ref.d, D=bg_ref.D,
                                                        B=bg_ref.B, d1=bg_ref.d1, D1=bg_ref.D1))
                pos_m += bg_ref.B


        return [plot_data.PrimitiveGroup(export_data)]

    def plot_contour2D(self, pos=0, a=None, box=True, typ='Graph'):
        be_sup = self.external_bearing(sign=1)
        be_inf = self.external_bearing(sign=-1)
        bi_sup = self.internal_bearing(sign=1)
        bi_inf = self.internal_bearing(sign=-1)
        contour = [be_sup, be_inf, bi_sup, bi_inf]
        linkage_area = vm.wires.Contour2D(contour)
        linkage_area = linkage_area.translation(vm.Vector2D(pos, 0), True)

        contour = []
        pos_m = -self.B/2.
        for bg, di in zip(self.bearings, self.directions):
            cont = bg.plot_contour(direction=di)
            cont = cont.translation(vm.Vector2D(pos_m + bg.B/2., 0), True)
            pos_m += bg.B
            contour.append(cont)
        assembly_bg = vm.wires.Contour2D(contour)
        assembly_bg = assembly_bg.translation(vm.Vector2D(pos, 0), True)

        return linkage_area, assembly_bg

    def solve_axial_load(self):
        shaft = unidimensional.Body(0, -1, name='Shaft')
        ground = unidimensional.Body(0, 0.5, name='Ground')

        p_shaft = unidimensional.Load(shaft, 500)
        id_ground = unidimensional.ImposedDisplacement(ground, 0.)
        imposed_displacements = [id_ground]
        loads = [p_shaft]


        component, nonlinear_linkages = self.elementary_axial_load(ground, shaft, 0)
        bodies = [ground, shaft]
        for bir, bor in component:
            bodies.append(bir)
            bodies.append(bor)

        sm = unidimensional.UnidimensionalModel(bodies, [], nonlinear_linkages, loads,
                                                imposed_displacements)
        result = sm.Solve(500)
        result.Plot(intensity_factor=1e-5)

    def elementary_axial_load(self, ground, shaft, pos, radial_load, bearing_result, axial_load=None):
        nb_bg_radial = sum([1 if p  else 0 for p in self.radial_load_linkage])

        component = []
        nonlinear_linkages = []
        axial_bearings = []
        loads = []
        k1 = 1e4
        j1 = 0
        posx = pos
        check_axial_load = True

        global_axial_load = 0



        for num_bg, (bg, bg_result) in enumerate(zip(self.bearings, bearing_result)):
            component_item = []
            component_item.append(unidimensional.Body(posx, bg.d/2., name='Inner ring{}'.format(num_bg)))
            component_item.append(unidimensional.Body(posx, bg.D/2., name='Outer ring{}'.format(num_bg)))
            bir = component_item[0]
            bor = component_item[1]

            if bg.taking_loads == 'both':
                pos1 = bir.initial_position
                pos2 = bor.initial_position
                link1 = unidimensional.CompressionSpring(bir, bor, k1, -j1, 'bearing {}'.format(num_bg))
                link2 = unidimensional.CompressionSpring(bor, bir, k1, -j1, 'bearing {}'.format(num_bg))
                nonlinear_linkages.append(link1)
                nonlinear_linkages.append(link2)

                axial_bearings.append([link1, link2])
#            elif bg.taking_loads == 'free':
#                link1 = unidimensional.CompressionSpring(bir, bor, 10, -j1, 'bearing {}'.format(num_bg))
##                link2 = unidimensional.CompressionSpring(bor, bir, 100, -j1, 'bearing {}'.format(num_bg))
#                nonlinear_linkages.append(link1)
##                nonlinear_linkages.append(link2)
#                axial_bearings.append([link1])
            elif bg.generate_axial_load:

                Fp = radial_load/nb_bg_radial*math.tan(bg.alpha)
                if self.directions[num_bg] == -1:
                    global_axial_load += Fp
                    link = unidimensional.CompressionSpring(bor, bir, k1, -j1, 'bearing {}'.format(num_bg))

                    nonlinear_linkages.append(link)


                    axial_bearings.append([link])

                    loads.append(unidimensional.Load(bor, Fp))
                    loads.append(unidimensional.Load(bir, -Fp))
                elif self.directions[num_bg] == 1:
                    global_axial_load += -Fp
                    link = unidimensional.CompressionSpring(bir, bor, k1, -j1, 'bearing {}'.format(num_bg))

                    nonlinear_linkages.append(link)

                    axial_bearings.append([link])
                    loads.append(unidimensional.Load(bor, -Fp))
                    loads.append(unidimensional.Load(bir, Fp))
            elif bg.taking_loads != 'free':
                if self.directions[num_bg] == -1:
                    link = unidimensional.CompressionSpring(bor, bir, k1, -j1, 'bearing {}'.format(num_bg))
                    nonlinear_linkages.append(link)

                    axial_bearings.append([link])
                elif self.directions[num_bg] == 1:

                    link = unidimensional.CompressionSpring(bir, bor, k1, -j1, 'bearing {}'.format(num_bg))
                    nonlinear_linkages.append(link)

                    axial_bearings.append([link])


            # check_radial_linkage = self.radial_load_linkage[num_bg]
            # if check_radial_linkage:
            #     if nb_bg_radial:
            #         bg_result.radial_load.append(radial_load/nb_bg_radial)
            component.append(component_item)
            posx += bg.B

            if axial_load is not None:
                global_axial_load += axial_load
                if (self.behavior_link == 'right') and (global_axial_load <= 0):
                    check_axial_load = False
                    break
                if (self.behavior_link == 'left') and (global_axial_load >= 0):
                    check_axial_load = False
                    break

        if len(component) > 1:
            for bg1, bg2 in zip(component[0:-1], component[1:]):
                pos1 = bg1[0].initial_position
                pos2 = bg2[0].initial_position
                nonlinear_linkages.append(unidimensional.UnilateralContact(bg1[0], bg2[0], pos2 - pos1, name='Inner rings'))
                nonlinear_linkages.append(unidimensional.UnilateralContact(bg1[1], bg2[1], pos2 - pos1, name='Outer rings'))
        if self.connection_be.left:

            bor = component[0][1]

            nonlinear_linkages.append(unidimensional.UnilateralContact(ground, bor, 0.0001, name='Outer rings'))
        if self.connection_be.right:

            bor = component[-1][1]

            nonlinear_linkages.append(unidimensional.UnilateralContact(bor, ground, 0.0001, name='Outer rings'))
        if self.connection_bi.left:

            bir = component[0][0]

            nonlinear_linkages.append(unidimensional.UnilateralContact(shaft, bir, 0.0001, name='Inner rings'))
        if self.connection_bi.right:

            bir = component[-1][0]

            nonlinear_linkages.append(unidimensional.UnilateralContact(bir, shaft, 0.0001, name='Inner rings'))

        return component, nonlinear_linkages, loads, axial_bearings, check_axial_load

    @classmethod
    def estimate_base_life_time(cls, L10s):
        sum_L10_inv = 0
        for L10 in L10s:
            sum_L10_inv += (1/L10)**1.5
        if sum_L10_inv == 0:
            return 0
        return sum_L10_inv**(-1/1.5)

    def base_life_time(self, bearing_combination_simulation_result, S=0.9):

        for bearing_result in bearing_combination_simulation_result.bearing_simulation_results:
            bearing_result.radial_load = []
            bearing_result.axial_load = []

        nb_bg_radial = sum([1 if (p is True) else 0 for p in self.radial_load_linkage])
        result_bgs = bearing_combination_simulation_result.bearing_simulation_results

        for radial_load, axial_load in zip(bearing_combination_simulation_result.radial_loads,
                                           bearing_combination_simulation_result.axial_loads):

            if (not self.behavior_link.free) and (abs(axial_load) >= 1e-4):


                check_axial_load = self.axial_load(axial_load, radial_load, bearing_combination_simulation_result)

                if check_axial_load == False:
                    return False
            else:

                for num_bg, bearing_result in enumerate(bearing_combination_simulation_result.bearing_simulation_results):
                    check_radial_linkage = self.radial_load_linkage[num_bg]

                    if check_radial_linkage:
                        bearing_result.radial_load.append(radial_load/nb_bg_radial)
                    bearing_result.axial_load.append(0)


        time = bearing_combination_simulation_result.operating_times
        speed = bearing_combination_simulation_result.speeds

        for bg, bg_result in zip(self.bearings,
                                 bearing_combination_simulation_result.bearing_simulation_results):
            try:
                L10 = bg.adjusted_life_time(Fr=bg_result.radial_load, Fa=bg_result.axial_load,
                                            N=speed, t=time, Cr=bg.Cr, S=S)
            except FloatingPointError:
                L10 = False

            except BearingL10Error:
                L10 = False
            if (str(L10) != 'nan') and (L10 != False):
                bg_result.L10 = L10
            else:
                bg_result.L10 = False

        sum_L10_inv = 0
        valid_L10 = True
        for bg_result in result_bgs:
            if (bg_result.L10 is not False) and (bg_result.L10 != 0):
                sum_L10_inv += (1/bg_result.L10)**1.5
            else:
                valid_L10 = False
        if valid_L10:
            bearing_combination_simulation_result.L10 = sum_L10_inv**(-1/1.5)
        else:
            bearing_combination_simulation_result.L10 = False

    def axial_load(self, axial_load, radial_load, bearing_combination_simulation_result):

        result_bgs = bearing_combination_simulation_result.bearing_simulation_results

        shaft = unidimensional.Body(0, 0, name='Shaft')
        ground = unidimensional.Body(0, 0.05, name='Ground')

        p_shaft = unidimensional.Load(shaft, axial_load)
        id_ground = unidimensional.ImposedDisplacement(ground, 0.)
        imposed_displacements = [id_ground]
        loads = [p_shaft]

        bodies = [ground, shaft]
        nonlinear_linkages = []

        component, nonlinear_linkages_iter, loads_iter, axial_bearings, check_axial_load \
            = self.elementary_axial_load(ground, shaft, 0, radial_load, result_bgs, axial_load)



        loads = loads + loads_iter

        for bir, bor in component:
            bodies.append(bir)
            bodies.append(bor)



        nonlinear_linkages.extend(nonlinear_linkages_iter)



        sm = unidimensional.UnidimensionalModel(bodies, [], nonlinear_linkages, loads,
                                                imposed_displacements)




        if check_axial_load:
            try:
                result_sm = sm.Solve(10)
    #            bearing_combination_simulation_result.axial_load_model = result_sm

                for num_bg, (axial_linkage, (bir, bor)) in enumerate(zip(axial_bearings, component)):
                    for link in axial_linkage:
                        if link in result_sm.activated_nonlinear_linkages:
                            positions = (result_sm.positions[bir], result_sm.positions[bor])
                            result_bgs[num_bg].axial_load.append(link.Strains(positions))
                return True
            except unidimensional.ModelConvergenceError:
                print('Convergence Error')

                return False
                pass
            except IndexError:

                print('Convergence Error')

                return False
                pass
        else:
            return False

    def plot(self, pos=0, a=None, box=True, typ='Load', ind_load_case=0):
        """
        Generate a Plot

        :param pos: axial position of bearing assembly center
        :param a: object define the append graph (default value is None)
        :param box: draw the box parameter of the bearing assembly
        :param typ: define the aditionnal draw (default is None), 'Graph' draw the graph connection between bearing, 'Load' define the load
        """

        plot_data.plot_canvas(self.plot_data()[0])
        # linkage_area, assembly_bg = self.plot_contour2D(pos, a, box, typ)
        # edge_style=plot_data.EdgeStyle(color_stroke=plot_data.colors.GREEN)
        # list_plot_data=[]
        # list_plot_data.append(linkage_area.plot_data(edge_style = edge_style))


        # edge_style=plot_data.EdgeStyle(color_stroke=plot_data.colors.BLACK)

        # list_plot_data.append(assembly_bg.plot_data(edge_style = edge_style))

        # if typ == 'Graph':
        #     list_graph = []
        #     pos_m = -self.B/2.
        #     # for bg_ref, bg_simu in zip(self.bearings, bearing_combination_result):
        #     #     graph = bg_simu.PlotGraph(d = bg_ref.d, D = bg_ref.D,
        #     #                           B = bg_ref.B, d1 = bg_ref.d1, D1 = bg_ref.D1)
        #     #     graph = graph.Translation(vm.Vector2D((pos_m + bg_ref.B/2., 0)), True)
        #     #     pos_m += bg_ref.B
        #     #     list_graph.append(graph)
        #     list_graph = vm.wires.Contour2D(list_graph)
        #     list_graph = list_graph.translation(vm.Vector2D(pos, 0), True)
        #     edge_style = plot_data.EdgeStyle(color_stroke=plot_data.colors.BLUE)
        #     list_plot_data.append(list_graph.plot_data(edge_style = edge_style))


        # elif typ == 'Load':
        #     pos_m = -self.B/2.

        #     max_load = 0
        #     for bg in self.bearings:
        #         for nd in bg.load_bearing_results[ind_load_case].list_node:
        #             if nd.load is not None:
        #                 max_load = max(nd.load, max_load)

        #     for bg in self.bearings:
        #         bg.PlotLoad(a, pos = pos + pos_m + bg.B/2., d = bg.d, D = bg.D,
        #                     B = bg.B, d1 = bg.d1, D1 = bg.D1,
        #                     ind_load_case = ind_load_case,
        #                     max_load = max_load)
        #         pos_m += bg.B

        # if box:
        #     box_sup = self.bearing_box(1)
        #     box_inf = self.bearing_box(-1)
        #     cont_box = [box_sup, box_inf]
        #     contour_box = vm.wires.Contour2D(cont_box)
        #     contour_box = contour_box.translation(vm.Vector2D(pos, 0), True)
        #     contour_box.MPLPlot(a, color='r')




    def volume_model(self, center=vm.Point3D(0, 0, 0), axis=vm.Vector3D(1, 0, 0)):
        groups = []
#        position = self.axial_positions
        center_bearing = center+0.5*(self.bearings[0].B - self.B)*axis
        for bearing, di in zip(self.bearings, self.directions):
            try:

                groups.extend(bearing.volmdlr_primitives(center=center_bearing, direction=di))
            except TypeError:
                groups.extend(bearing.volmdlr_primitives(center=center_bearing))
            center_bearing += bearing.B*axis


        model = vm.core.VolumeModel(groups, self.name)
        return model

    def volmdlr_volume_model(self):
        model = self.volume_model()
        return model

    def to_markdown(self):
        infos = ''
        for i, bearing in enumerate(self.bearings):
            infos += 'Bearing ' + str(i) + '\n\n'
            infos += bearing.to_markdown()

        return infos

#    def Dict(self, subobjects_id={}, stringify_keys=True):
#        """
#        Export dictionary
#        """
#        d={}
#        d['directions'] = self.directions
#
#        d['radial_load_linkage'] = self.radial_load_linkage
#        d['internal_pre_load'] = self.internal_pre_load
#        d['connection_bi'] = self.connection_bi
#        d['connection_be'] = self.connection_be
#        d['behavior_link'] = self.behavior_link
#
#        bearings = []
#        for bearing in self.bearings:
#            if bearing in subobjects_id:
#                bearings.append(subobjects_id[bearing])
#            else:
#                bearings.append(bearing.to_dict())
#        d['bearings'] = bearings
#
#        if stringify_keys:
#            return StringifyDictKeys(d)
#        return d
#
#    @classmethod
#    def dict_to_object(cls, d):
#        bearings = []
#        for bearing_s in d['bearings']:
#            object_class = bearing_s['object_class']
#            module = object_class.rsplit('.', 1)[0]
#            exec('import ' + module)
#            class_ = eval(object_class)
#            bearings.append(class_.dict_to_object(bearing_s))
#        obj = cls(bearings = bearings, directions = d['directions'], radial_load_linkage = d['radial_load_linkage'],
#                  internal_pre_load = 0, connection_bi = d['connection_bi'],
#                  connection_be = d['connection_be'], behavior_link = d['behavior_link'])
#        return obj

class BearingAssembly(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, bearing_combinations: List[BearingCombination],
                 pre_load: float = 0,
                 axial_positions: List[float] = None,
                 internal_diameters: List[float] = None, axial_pos: List[float] = None,
                 external_diameters: List[float] = None, length: float = None,
                 overall_length: float = None, mass: float = None,
                 cost: float = None,
                 name: str = ''):

        self.bearing_combinations = bearing_combinations
        self.cr_equ = self.cr_equ()
        self.pre_load = pre_load
        self.load_bearing_assembly_results = None
        self.axial_positions = axial_positions
        self.internal_diameters = internal_diameters
        self.axial_pos = axial_pos
        self.external_diameters = external_diameters
        self.length = length
        self.B = 0
        for bc in bearing_combinations:
            self.B += bc.B
        self.D = 0
        for bc in bearing_combinations:
            self.D = max(bc.D, self.D)
        self.d = 0
        for bc in bearing_combinations:
            self.d = max(bc.d, self.d)
        if overall_length is None:
            self.overall_length = self.B
            if axial_positions is not None:
                self.overall_length += axial_positions[1] - axial_positions[0]
        else:
            self.overall_length = overall_length

        self.number_bearing = sum([bc.number_bearing for bc in self.bearing_combinations])
        self.number_bearing_first_bc = self.bearing_combinations[0].number_bearing
        self.number_bearing_second_bc = self.bearing_combinations[1].number_bearing

        DessiaObject.__init__(self, name=name)

#    def __eq__(self, other_eb):
#        equal = True
#        for bc, other_bc in zip(self.bearing_combinations, other_eb.bearing_combinations):
#            equal = equal and bc == other_bc
#        equal = equal and self.pre_load == other_eb.pre_load
#        if hasattr(self, 'axial_positions') and hasattr(other_eb, 'axial_positions'):
#            equal = equal and list(self.axial_positions) == list(other_eb.axial_positions)
#        elif hasattr(self, 'axial_positions') or hasattr(other_eb, 'axial_positions'):
#            equal = False
#        return equal
#
#    def __hash__(self):
#        h = 0
#        for bc in self.bearing_combinations:
#            h += hash(bc)
#        return h

    def update(self, axial_positions, internal_diameters, axial_pos,
               external_diameters, length):
#        if axial_position is not None:
#            axial_position = list(axial_position)
        self.axial_positions = axial_positions
        self.internal_diameters = internal_diameters
        self.axial_pos = axial_pos
        self.external_diameters = external_diameters
        self.length = length
        for num_linkage, assembly_bg in enumerate(self.bearing_combinations):
            pos = self.axial_pos[num_linkage] - self.axial_positions[num_linkage]
            assembly_bg.update(pos, self.internal_diameters[num_linkage], self.external_diameters[num_linkage],
                               self.length[num_linkage])
        self.overall_length = self.B + axial_positions[1] - axial_positions[0]

    def cr_equ(self):
        Cr_equ = 0
        for li_bg in self.bearing_combinations:
            for bg in li_bg.bearings:
                Cr_equ += (bg.Cr)
        return (Cr_equ)

    @property
    def mass(self):
        mass = 0
        for li_bg in self.bearing_combinations:
            mass += li_bg.mass
        if self.axial_positions is not None:
            poly2d = self.shaft()
            mass += poly2d.area()
        return mass

    @property
    def cost(self):
        cost = 0
        for bc in self.bearing_combinations:
            cost += bc.cost
        return cost

    def shaft(self):
        d1 = self.bearing_combinations[0].d
        d2 = self.bearing_combinations[1].d
        B1 = self.bearing_combinations[0].B
        B2 = self.bearing_combinations[1].B
        pos_mid = (self.axial_positions[0] + self.axial_positions[1])/2.
        shaft = primitives2d.ClosedRoundedLineSegments2D([vm.Point2D(self.axial_positions[0] - B1/2. - 5e-3, -d1/2.),
                                                          vm.Point2D(self.axial_positions[0] - B1/2. - 5e-3, d1/2.),
                                                          vm.Point2D(pos_mid, d1/2.),
                                                          vm.Point2D(self.axial_positions[1] + B2/2. + 5e-3, d2/2.),
                                                          vm.Point2D(self.axial_positions[1] + B2/2. + 5e-3, -d2/2.),
                                                          vm.Point2D(pos_mid, -d2/2.)
                                                          ],
                                                         {}, adapt_radius=False)
        # shaft = vm.Polygon2D([vm.Point2D((self.axial_positions[0] - B1/2. - 5e-3, -d1/2.)),
        #               vm.Point2D((self.axial_positions[0] - B1/2. - 5e-3, d1/2.)),
        #               vm.Point2D((pos_mid, d1/2.)),
        #               vm.Point2D((pos_mid, d2/2.)),
        #               vm.Point2D((self.axial_positions[1] + B2/2. + 5e-3, d2/2.)),
        #               vm.Point2D((self.axial_positions[1] + B2/2. + 5e-3, -d2/2.)),
        #               vm.Point2D((pos_mid, -d2/2.)),
        #               vm.Point2D((pos_mid, -d1/2.)),
        #               vm.Point2D((self.axial_positions[0] - B1/2. - 5e-3, -d1/2.))])
        return vm.wires.Contour2D(shaft.primitives)

    def cad_shaft(self, center=vm.O3D, axis=vm.X3D):
        # TODO: mutualization of this in parent class?
        axis.normalize()

        y = axis.random_unit_normal_vector()
        z = axis.cross(y)

        #Internal Ring
        d1 = self.bearing_combinations[0].d
        d2 = self.bearing_combinations[1].d
        B1 = self.bearing_combinations[0].B
        B2 = self.bearing_combinations[1].B
        pos_mid = (self.axial_positions[0] + self.axial_positions[1])/2.
        if d1 == d2:
            p1 = vm.Point2D(self.axial_positions[0] - B1/2. - 5e-3, 0)
            p2 = vm.Point2D(self.axial_positions[0] - B1/2. - 5e-3, d1/2.)
            p3 = vm.Point2D(pos_mid, d1/2.)
            p5 = vm.Point2D(self.axial_positions[1] + B2/2. + 5e-3, d2/2.)
            p6 = vm.Point2D(self.axial_positions[1] + B2/2. + 5e-3, 0)
            shaft = primitives2d.OpenedRoundedLineSegments2D([p6, p5, p3, p2, p1],
                                                             {},
                                                             adapt_radius=True)
        else:
            p1 = vm.Point2D(self.axial_positions[0] - B1/2. - 5e-3, 0)
            p2 = vm.Point2D(self.axial_positions[0] - B1/2. - 5e-3, d1/2.)
            p3 = vm.Point2D(pos_mid, d1/2.)
            p4 = vm.Point2D(pos_mid, d2/2.)
            p5 = vm.Point2D(self.axial_positions[1] + B2/2. + 5e-3, d2/2.)
            p6 = vm.Point2D(self.axial_positions[1] + B2/2. + 5e-3, 0)
            shaft = primitives2d.OpenedRoundedLineSegments2D([p6, p5, p4, p3, p2, p1],
                                                             {},
                                                             adapt_radius=True)
        l1 = vm.edges.LineSegment2D(p1, p6)
        shaft_cont = vm.wires.Contour2D(shaft.primitives + [l1])

        irc = primitives3d.RevolvedProfile(center, axis, z, shaft_cont, center, axis, angle=2*math.pi, color=[204/255, 12/255, 12/255], name='Shaft')

        return [irc]

    def plot_data(self, box=True, typ=None, constructor=False):

        shaft = self.shaft()
        # contour_shaft = vm.wires.Contour2D(shaft.primitives)
        export_data = [shaft.plot_data()]


        for assembly_bg, pos in zip(self.bearing_combinations, self.axial_positions):
            export_data.extend(assembly_bg.plot_data(pos, box, quote=False, constructor=constructor)[0].primitives)
        return [plot_data.PrimitiveGroup(primitives=export_data)]

    # def plot(self, box=True, typ=None, ind_load_case=0):

    #     shaft = self.shaft()
    #     # contour_shaft = vm.Contour2D([shaft])
    #     f, a = shaft.MPLPlot()

    #     for assembly_bg, pos in zip(self.bearing_combinations, self.axial_positions):
    #         assembly_bg.plot(pos, a, box, typ, ind_load_case)

    def graph(self):
        for li_bg in self.bearing_combinations:
            G, positions, li_axial_link, nd_axial_load = li_bg.Graph()
#            plt.figure()
#            nx.draw_networkx(G, pos = positions)

    def check_load(self, bearing_assembly_simulation_result):
        loads = bearing_assembly_simulation_result.loads
        valid_axial_load = True
        for ind_load_case, load_cases in enumerate(loads):

            axial_load = 0
            for (pos, ld, tq) in load_cases:
                axial_load += ld[0]
            valid = False
            for bc in self.bearing_combinations:
                if axial_load > 0:
                    if 'right' in bc.transfert_load:
                        valid = True
                elif axial_load < 0:
                    if 'left' in bc.transfert_load:
                        valid = True
            if not valid:
                valid_axial_load = False
        return valid_axial_load

    @classmethod
    def quick_shaft_load(cls, positions, loads):

        (pos1, pos2) = positions
        ground = genmechanics.Part('ground')
        shaft1 = genmechanics.Part('shaft1')
        p1 = npy.array([pos1, 0, 0])
        p2 = npy.array([pos2, 0, 0])
        bearing1 = linkages.FrictionlessBallLinkage(ground, shaft1, p1, [0, 0, 0], 'bearing1')
        bearing2 = linkages.FrictionlessLinearAnnularLinkage(ground, shaft1, p2, [0, 0, 0], 'bearing2')

        bearing_combination_results = []


        for ind_load_case, load_cases in enumerate(loads):
            load1 = []

            for pos, ld, tq in load_cases:

                load1.append(gm_loads.KnownLoad(shaft1, pos, [0, 0, 0], ld, tq, 'input'))
            load2 = gm_loads.SimpleUnknownLoad(shaft1, [(pos1 + pos2)/2, 0, 0], [0, 0, 0], [], [0], 'output torque')

            imposed_speeds = [(bearing1, 0, 100)]
        
            mech = genmechanics.Mechanism([bearing1, bearing2], ground, imposed_speeds, load1, [load2])

            axial_load = 0
            for pos, ld, tq in load_cases:
                axial_load += ld[0]
            bearing_combination_result = []
            for bg in [bearing1, bearing2]:
                tensor = mech.GlobalLinkageForces(bg, 1)
                fr = (tensor[1]**2 + tensor[2]**2)**(0.5)
                bearing_combination_result.append(fr)
            bearing_combination_result.append(axial_load)
            bearing_combination_results.append(bearing_combination_result)
        return bearing_combination_results

    def shaft_load(self, positions, bearing_assembly_simulation_result, fr_bearings = None, verbose=False):

        result_bcs = bearing_assembly_simulation_result.bearing_combination_simulation_results
        result_bgs = [bg.bearing_simulation_results for bg in result_bcs]
        for bearing_results, bearing_combination_result in zip(result_bgs, result_bcs):
            bearing_combination_result.radial_loads = []
            bearing_combination_result.axial_loads = []
            for bearing_result in bearing_results:
                bearing_result.radial_load = []
                bearing_result.axial_load = []

        loads = bearing_assembly_simulation_result.loads
        (pos1, pos2) = positions
        if not fr_bearings:
            ground = genmechanics.Part('ground')
            shaft1 = genmechanics.Part('shaft1')
            p1 = npy.array([pos1, 0, 0])
            p2 = npy.array([pos2, 0, 0])
            bearing1 = linkages.FrictionlessBallLinkage(ground, shaft1, p1, [0, 0, 0], 'bearing1')
            bearing2 = linkages.FrictionlessLinearAnnularLinkage(ground, shaft1, p2, [0, 0, 0], 'bearing2')
    
        for ind_load_case, load_cases in enumerate(loads):
            axial_load = 0
            for pos, ld, tq in load_cases:
                axial_load += ld[0]
            if not fr_bearings:
                load1 = []
                for pos, ld, tq in load_cases:
                    load1.append(gm_loads.KnownLoad(shaft1, pos, [0, 0, 0], ld, tq, 'input'))
                load2 = gm_loads.SimpleUnknownLoad(shaft1, [(pos1 + pos2)/2, 0, 0], [0, 0, 0], [], [0], 'output torque')
                imposed_speeds = [(bearing1, 0, 100)]
    
                mech = genmechanics.Mechanism([bearing1, bearing2], ground, imposed_speeds, load1, [load2])

               
                for li_bg, bg, bearing_combination_result in zip(self.bearing_combinations,
                                                                 [bearing1, bearing2],
                                                                 result_bcs):
                    tensor = mech.LocalLinkageForces(bg, 1)
                    fr = (tensor[1]**2 + tensor[2]**2)**(0.5)
                    bearing_combination_result.radial_loads.append(fr)
                    
            else:
    
                for indice_bearing, bearing_combination_result in enumerate(result_bcs):

                    fr = fr_bearings[ind_load_case][indice_bearing]
                    bearing_combination_result.radial_loads.append(fr)
                    
            for num_linkage, (bearing_combination, load_bearing_combination_result) in enumerate(zip(self.bearing_combinations,
                                                                                                     result_bcs)):
                radial_load = load_bearing_combination_result.radial_loads[-1]
                bearing_result = result_bgs[num_linkage]
                nb_bg_radial = sum([1 if p  else 0 for p in bearing_combination.radial_load_linkage])
                for num_bg, (bg, bg_result) in enumerate(zip(bearing_combination.bearings, bearing_result)):
                    
                    check_radial_linkage = bearing_combination.radial_load_linkage[num_bg]
                    if check_radial_linkage:
                        if nb_bg_radial:
                            bg_result.radial_load.append(radial_load/nb_bg_radial)

            
            bearing_assembly_simulation_result = self.axial_load(positions, axial_load, bearing_assembly_simulation_result)



        try:

            
            self.base_life_time(bearing_assembly_simulation_result, verbose)
        except BearingL10Error:
            raise BearingL10Error()

    @classmethod
    def estimate_base_life_time(cls, L10s):
        sum_L10_inv = 0


        for L10 in L10s:
            sum_L10_inv += (1/L10)**1.5
        if sum_L10_inv == 0:

            return 0
        return sum_L10_inv**(-1/1.5)

    def base_life_time(self, bearing_assembly_simulation_result, verbose=False):

        result_bcs = bearing_assembly_simulation_result.bearing_combination_simulation_results
        result_bgs = [bg.bearing_simulation_results for bg in result_bcs]
       
        for bearing_combination, bearing_result in zip(self.bearing_combinations, result_bgs):
           
            for bg, bg_result in zip(bearing_combination.bearings, bearing_result):
           
                time = bearing_assembly_simulation_result.operating_times
                speed = bearing_assembly_simulation_result.speeds
                try:
                    # print(2879797)
                    # print(bearing_combination.bearings[0])

                    L10 = bg.base_life_time(Fr=bg_result.radial_load, Fa=bg_result.axial_load,
                                            N=speed, t=time, Cr=bg.Cr, verbose=verbose)


#                if (str(L10) != 'nan') and (L10 != False):
                    bg_result.L10 = L10
                except BearingL10Error:
                    raise BearingL10Error()

        sum_L10_inv = 0
        valid_L10 = True
        L10_min = 1000000000
        for bearing_result in result_bgs:
            for bg_result in bearing_result:
#                if bg_result.L10 != False:
                if bg_result.L10 < L10_min:
                    L10_min = bg_result.L10
                if (bg_result.L10 is not False) and (bg_result.L10 != 0):
                    sum_L10_inv += (1/bg_result.L10)**1.5
                else:
                    valid_L10 = False

        bearing_assembly_simulation_result.L10 = L10_min

    def axial_load(self, positions, axial_load, bearing_assembly_simulation_result):
        result_bcs = bearing_assembly_simulation_result.bearing_combination_simulation_results
        result_bgs = [bg.bearing_simulation_results for bg in result_bcs]

        shaft = unidimensional.Body(0, 0, name='Shaft')
        ground = unidimensional.Body(0, 0.05, name='Ground')

        p_shaft = unidimensional.Load(shaft, axial_load)
        id_ground = unidimensional.ImposedDisplacement(ground, 0.)
        imposed_displacements = [id_ground]
        loads = [p_shaft]
        bc_axial_bearings = []

        bodies = [ground, shaft]
        nonlinear_linkages = []
        components = []
        Fp = []
        generate_axial_load = False
        for bearing_combination in self.bearing_combinations:
            for bearing in bearing_combination.bearings:
                if bearing.generate_axial_load:
                    generate_axial_load = True

        if abs(axial_load) > 1e-8 or generate_axial_load:
            for num_linkage, (bearing_combination, load_bearing_combination_result) in enumerate(zip(self.bearing_combinations,
                                                                                                     result_bcs)):
    
                pos = positions[num_linkage]
                radial_load = load_bearing_combination_result.radial_loads[-1]
                bearing_result = result_bgs[num_linkage]
    
                component, nonlinear_linkages_iter, loads_iter, axial_bearings, __\
                    = bearing_combination.elementary_axial_load(ground, shaft, pos, \
                                                              radial_load, bearing_result)
    
                nb_bg_radial = sum([1 if p  else 0 for p in bearing_combination.radial_load_linkage])
                Fp.append(radial_load/nb_bg_radial*math.tan(bearing_combination.bearings[0].alpha))
                if (not bearing_combination.behavior_link.free) and (axial_load != 0):
                    bc_axial_bearings.append(axial_bearings)
                    loads = loads + loads_iter
    
                    components.append(component)
                    for bir, bor in component:
                        bodies.append(bir)
                        bodies.append(bor)
                    nonlinear_linkages.extend(nonlinear_linkages_iter)
                else:
                    bc_axial_bearings.append([])
                    components.append([])
    
            sm = unidimensional.UnidimensionalModel(bodies, [], nonlinear_linkages, loads,
                                                    imposed_displacements)
    
    
    
            try:
                result_sm = sm.Solve(500)
    
    
    #            bearing_assembly_simulation_result.axial_load_model = result_sm
    
                for num_bc, (axial_linkages, component) in enumerate(zip(bc_axial_bearings, components)):
    
                    if len(axial_linkages) == 0:
                        for result_bg in result_bgs[num_bc]:
                            result_bg.axial_load.append(0)
                    for num_bg, (axial_linkage, (bir, bor)) in enumerate(zip(axial_linkages, component)):
                        link_axial_load = False
                        for link in axial_linkage:
    
                            if link in result_sm.activated_nonlinear_linkages:
                                positions = (result_sm.positions[bir], result_sm.positions[bor])
                                # if len(Fp)==2:
                                #     if axial_load<0:
                                #         bearing_position_1=1
                                #         bearing_position_2=0
                                #     else:
                                #         bearing_position_1=0
                                #         bearing_position_2=1
                                #     if abs(axial_load)+abs(Fp[bearing_position_1])>abs(Fp[bearing_position_2]):
    
                                #         result_bgs[num_bc][num_bg].axial_load.append(abs(axial_load)+abs(Fp[bearing_position_1]))
                                #     else:
                                #         result_bgs[num_bc-1][num_bg].axial_load.append(-abs(axial_load)+abs(Fp[bearing_position_2]))
    
                                # else:
                                result_bgs[num_bc][num_bg].axial_load.append(link.Strains(positions))
    
    
                                link_axial_load = True
                                break
                        if not link_axial_load:
                                result_bgs[num_bc][num_bg].axial_load.append(0.0)
    
                return bearing_assembly_simulation_result

            except unidimensional.ModelConvergenceError:
                print('Convergence Error')
    
                pass
        else:
        
            for bc in result_bgs:
                for bg in bc:
                    
                    bg.axial_load.append(0.0)
    
            return bearing_assembly_simulation_result
    def volmdlr_primitives(self):
        primitives = []
        cad_shaft = self.cad_shaft()
        primitives.extend(cad_shaft)

        y_pos = 0
        z_pos = 0

        for assembly_bg, x_pos in zip(self.bearing_combinations, self.axial_positions):
            center = vm.Point3D(x_pos, y_pos, z_pos)
            primitives.extend(assembly_bg.volume_model(center=center).primitives)

        return primitives

    def volume_model(self):
        groups = self.volmdlr_primitives()

        model = vm.core.VolumeModel(groups)

        return model


class BearingSimulationResult(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, axial_load: float = None, radial_load: float = None,
                 L10: float = None, name: str = ''):
        if axial_load is None:
            self.axial_load = []
        else:
            self.axial_load = axial_load
        if radial_load is None:
            self.radial_load = []
        else:
            self.radial_load = radial_load
        self.L10 = L10
        DessiaObject.__init__(self, name=name)



class BearingCombinationSimulationResult(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, bearing_simulation_results: List[BearingSimulationResult],
                 axial_loads: List[float] = None, radial_loads: List[float] = None,
                 speeds: List[float] = None, operating_times: List[float] = None,
                 max_axial_load: float = None, max_radial_load: float = None,
                 name: str = ''):
        if axial_loads is None:
            self.axial_loads = []
        else:
            self.axial_loads = axial_loads
        if radial_loads is None:
            self.radial_loads = []
        else:
            self.radial_loads = radial_loads
        self.bearing_simulation_results = bearing_simulation_results

        if speeds is not None:
            self.speeds = speeds
        if operating_times is not None:
            self.operating_times = operating_times
#        if axial_load_model is not None:
#            self.axial_load_model = axial_load_model

        DessiaObject.__init__(self, name=name)

    @property
    def max_axial_load(self):
        max_axial_load = 0
        for bearing_simulation_result in self.bearing_simulation_results:
            max_axial_load = max(max_axial_load, max(bearing_simulation_result.axial_load))
        return max_axial_load

    @property
    def max_radial_load(self):
        max_radial_load = 0
        for bearing_simulation_result in self.bearing_simulation_results:
            max_radial_load = max(max_radial_load, max(bearing_simulation_result.radial_load))
        return max_radial_load




class BearingAssemblySimulationResult(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, bearing_combination_simulation_results: List[BearingCombinationSimulationResult],
                 loads: List[float], speeds, operating_times,
                 bearing_combination_first: BearingCombinationSimulationResult = None,
                 bearing_combination_second: BearingCombinationSimulationResult = None,
                 L10: float = None, name: str = ''):
        self.loads = loads
        self.speeds = speeds
        self.operating_times = operating_times
#        self.axial_load_model = axial_load_model
        self.L10 = L10
        self.bearing_combination_simulation_results = bearing_combination_simulation_results
        self.bearing_combination_first = self.bearing_combination_simulation_results[0]
        self.bearing_combination_second = self.bearing_combination_simulation_results[1]

        DessiaObject.__init__(self, name=name)




class BearingAssemblySimulation(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, bearing_assembly: BearingAssembly,
                 bearing_assembly_simulation_result: BearingAssemblySimulationResult,
                 name: str = ''):
        self.bearing_assembly = bearing_assembly
        self.bearing_assembly_simulation_result = bearing_assembly_simulation_result

        DessiaObject.__init__(self, name=name)

        self.overall_length = self.bearing_assembly.overall_length
        self.mass = self.bearing_assembly.mass
        self.cost = self.bearing_assembly.cost
        self.number_bearing = self.bearing_assembly.number_bearing
        self.number_bearing_first_bc = self.bearing_assembly.number_bearing_first_bc
        self.number_bearing_second_bc = self.bearing_assembly.number_bearing_second_bc
        self.first_max_axial_load = self.bearing_assembly_simulation_result.bearing_combination_first.max_axial_load
        self.first_max_radial_load = self.bearing_assembly_simulation_result.bearing_combination_first.max_radial_load
        self.second_max_axial_load = self.bearing_assembly_simulation_result.bearing_combination_second.max_axial_load
        self.second_max_radial_load = self.bearing_assembly_simulation_result.bearing_combination_second.max_radial_load

    def volmdlr_volume_model(self):
        model = self.bearing_assembly.volume_model()
        return model
    def volmdlr_primitives(self):
        volmdlr_primitives = self.bearing_assembly.volmdlr_primitives()
        return volmdlr_primitives

    def plot_data(self, box=True, typ=None, constructor=False):
        plot_data = self.bearing_assembly.plot_data()
        return plot_data



class BearingCombinationSimulation(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, bearing_combination: BearingCombination,
                 bearing_combination_simulation_result: BearingCombinationSimulationResult,
                 name: str = ''):
        self.bearing_combination = bearing_combination
        self.bearing_combination_simulation_result = bearing_combination_simulation_result

        DessiaObject.__init__(self, name=name)

        self.mass = self.bearing_combination.mass
        self.cost = self.bearing_combination.cost
        self.number_bearing = self.bearing_combination.number_bearing
        self.max_axial_load = self.bearing_combination_simulation_result.max_axial_load
        self.max_radial_load = self.bearing_combination_simulation_result.max_radial_load



    def volmdlr_volume_model(self):
        model = self.bearing_combination.volume_model()
        return model

    def plot_data(self, box=True, typ=None, constructor=False):
        plot_data = self.bearing_combination.plot_data(box=box)
        return plot_data




class BearingL10Error(Exception):
    def __init__(self):
        super().__init__('L10 simulation failed due to dynamic equivalent load equal 0')

class CatalogSearchError(Exception):
    def __init__(self):
        super().__init__('No available bearings with d and D defined')
