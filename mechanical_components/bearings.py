import numpy as npy
npy.seterr(divide='raise', over='ignore', under='ignore')
#import math as mt
from scipy import interpolate
#import os
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math
from mechanical_components import shafts_assembly 

#import copy
import json

from scipy.optimize import fsolve
import networkx as nx
import matplotlib.pyplot as plt

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

#oil_kinematic_viscosity
iso_vg_1500={'data':[[47.21238870380181,922.5481847223729],
                     [76.41592953982855,191.5471560481642],
                     [110.70796589064605,54.90426918109079]],'x':'Linear','y':'Log'}
iso_vg_1000={'data':[[41.68141577143845,877.2173704075102],
                     [62.477877333256444,261.78400435754804],
                     [100.53097486106454,57.74149074755608]],'x':'Linear','y':'Log'}
iso_vg_680={'data':[[38.80530942959304,777.3038394524846],
                    [57.168142067138206,251.441902731903],
                    [89.46902691125547,61.96159004047747]],'x':'Linear','y':'Log'}
iso_vg_460={'data':[[36.15044283907511,580.3315115122488],
                    [59.159291488756054,159.77215436382392],
                    [85.48672598293739,53.80881018274548]],'x':'Linear','y':'Log'}
iso_vg_320={'data':[[32.16814191075703,551.8160309554283],
                    [57.38937973338331,131.93199565920736],
                    [81.06194763704671,48.65076991453211]],'x':'Linear','y':'Log'}
iso_vg_220={'data':[[29.95575273781169,407.8526478060212],
                    [56.725664649565616,96.53454045940936],
                    [83.05309705866458,35.24081769455843]],'x':'Linear','y':'Log'}
iso_vg_150={'data':[[27.964601231111455,307.58489032135725],
                    [50.97345196587479,89.9597273365993],
                    [87.25663773831015,23.313898800373792]],'x':'Linear','y':'Log'}
iso_vg_100={'data':[[33.05309674590221,148.89034266049572],
                    [60.7079655778837,41.82579554586569],
                    [91.23893866662823,15.115797660575524]],'x':'Linear','y':'Log'}
iso_vg_68={'data':[[29.95575273781169,113.42390186278217],
                   [56.94690231581072,34.19139868782362],
                   [89.91150432882806,11.749567781125915]],'x':'Linear','y':'Log'}
iso_vg_46={'data':[[29.070795817584123,76.56429373059628],
                   [60.48672582655621,21.946066333596434],
                   [99.4247781895095,7.316992362396092]],'x':'Linear','y':'Log'}
iso_vg_32={'data':[[27.52212381353886,56.0220352459023],
                   [58.27433665361086,17.058761946001017],
                   [82.16814222351938,8.510952177980519]],'x':'Linear','y':'Log'}
iso_vg_22={'data':[[30.619469906711767,32.840621976693456],
                   [57.38937973338331,12.481883087082235],
                   [90.79646124905564,4.939173694948054]],'x':'Linear','y':'Log'}
iso_vg_15={'data':[[23.982300928318086,28.519522512213047],
                   [44.115044695711276,13.126893028298408],
                   [77.07964670872863,4.889651598606255]],'x':'Linear','y':'Log'}
iso_vg_10={'data':[[25.088495514790754,17.231530421142708],
                   [46.548673619984115,8.092752773048398],
                   [66.01769875891955,4.649391098015179]],'x':'Linear','y':'Log'}


#Ordre de rangement du coefficient de contamination de l'huile: Dpw mini/ Dpw maxi/ grade/ coeff de contamination
dict_oil_contamination={0:{0.1:{1:1,2:0.7,3:0.55,4:0.4,5:0.2,6:0.05,7:0}},
                        0.1:{math.inf:{1:1,2:0.85,3:0.7,4:0.5,5:0.3,6:0.05,7:0}}}


class Oil:
    
    def __init__(self, oil_data):
        self.oil_data = oil_data
        self.oil_kinematic_viscosity_curve = self.KinematicViscosity(oil_data)
    
    def FunCoeff(self,x,data,type_x='Linear',type_y='Linear'):
        if type_x == 'Log': 
            x = math.log10(x)
        f = interpolate.interp1d(list(data[:,0]),list(data[:,1]), fill_value='extrapolate')
        sol = float(f(x))
        if type_y == 'Log':
            sol = 10**sol
        return sol
    
    def KinematicViscosity(self,oil_kinematic_viscosity):
        oil_kinematic_viscosity_curve = {}
        for (key,val) in oil_kinematic_viscosity.items():
            val_np = npy.array(val)
            if key not in ['x','y']:
                oil_kinematic_viscosity_curve[key] = {}
                A = (math.log10(math.log10(0.6+val_np[0,1]))-math.log10(math.log10(0.6+val_np[-1,1])))/(math.log10(val_np[0,0])-math.log10(val_np[-1,0]))
                B = math.log10(math.log10(0.6+val_np[0,1]))-A*math.log10(val_np[0,0])
                oil_kinematic_viscosity_curve[key]['A'] = A
                oil_kinematic_viscosity_curve[key]['B'] = B
        return oil_kinematic_viscosity_curve['data']
    
    def OilParameterContamination(self,Dpw,grade):
        for k,v in dict_oil_contamination.items():
            if (Dpw>=k) and (Dpw<list(v.keys())[0]):
                return list(v.values())[0][grade]
            
    def Dict(self):
        """
        Export dictionary
        """
        d={}
        for k,v in self.__dict__.items():
            tv=type(v)
            if tv == npy.int64:
                d[k]=int(v)
            elif tv==npy.float64:
                d[k]=float(v)
            else:
                d[k]=v
            
        return d
    
    @classmethod
    def DictToObject(cls, d):
        obj = cls(oil_data = d['oil_data'])
        return obj
    
oil_iso_vg_1500=Oil(iso_vg_1500)
oil_iso_vg_1000=Oil(iso_vg_1000)
oil_iso_vg_680=Oil(iso_vg_680)
oil_iso_vg_460=Oil(iso_vg_460)
oil_iso_vg_320=Oil(iso_vg_320)
oil_iso_vg_220=Oil(iso_vg_220)
oil_iso_vg_150=Oil(iso_vg_150)
oil_iso_vg_100=Oil(iso_vg_100)
oil_iso_vg_68=Oil(iso_vg_68)
oil_iso_vg_46=Oil(iso_vg_46)
oil_iso_vg_32=Oil(iso_vg_32)
oil_iso_vg_22=Oil(iso_vg_22)
oil_iso_vg_15=Oil(iso_vg_15)
oil_iso_vg_10=Oil(iso_vg_10)
        
class Material:
    def __init__(self, weibull_e=9/8., weibull_c=31/3., weibull_h=7/3.,
                 B1=551.13373/0.483, mu_delta=0.83, c_gamma=0.05):

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


    def Dict(self):
        """Export dictionary
        """
        d={}
        for k,v in self.__dict__.items():
            tv=type(v)
            if tv==npy.int64:
                d[k]=int(v)
            elif tv==npy.float64:
                d[k]=float(v)
            else:
                d[k]=v
            
        return d
    
    @classmethod
    def DictToObject(cls, d):
        obj = cls(weibull_e = d['weibull_e'], weibull_c = d['weibull_c'], 
                  weibull_h = d['weibull_h'],
                  B1 = d['B1'], mu_delta = d['mu_delta'], c_gamma = d['c_gamma'])
        return obj

material_iso=Material()

#class LoadBearing:
#    def __init__(self, class_name, typ=None, direction=None):
#        self.class_name = class_name
#        self.typ = typ
#        self.direction = direction
#                
#    def PlotGraph(self, d, D, B, d1, D1, ind_load_case=0):
#        
#        list_node = self.load_bearing_results[ind_load_case].list_node
#        delta_1 = (D - D1)/10.
#        delta_2 = (d1 - d)/10.
#        positions = {}
#        positions[list_node[0]] = vm.Point2D((-B/2., D/2. - delta_1))
#        positions[list_node[1]] = vm.Point2D((-B/2., D/2. - 2*delta_1))
#        positions[list_node[2]] = vm.Point2D((-B/2., d/2. + 2*delta_2))
#        positions[list_node[3]] = vm.Point2D((-B/2., d/2. + delta_2))
#        positions[list_node[4]] = vm.Point2D((B/2., D/2. - delta_1))
#        positions[list_node[5]] = vm.Point2D((B/2., D/2. - 2*delta_1))
#        positions[list_node[6]] = vm.Point2D((B/2., d/2. + 2*delta_2))
#        positions[list_node[7]] = vm.Point2D((B/2., d/2. + delta_2))
#        
#        list_line = []
#        check_line = []
#        for nd in list_node:
#            np = self.load_bearing_results[ind_load_case].next[nd]
#            for n in np:
#                if (n in list_node) and [nd, n] not in check_line:
#                    n1 = positions[nd]
#                    n2 = positions[n]
#                    list_line.append(vm.LineSegment2D(n1, n2))
#                    check_line.append([nd, n])
#                    
#        load_arrow = vm.Contour2D(list_line)
#        return load_arrow
#                
#    def PlotLoad(self, a, pos, d, D, B, d1, D1, ind_load_case=0, max_load=None):
#        
#        list_node = self.load_bearing_results[ind_load_case].list_node
#        transversale_load = self.load_bearing_results[ind_load_case].transversale_load
#        internal_ring_load = self.load_bearing_results[ind_load_case].internal_ring_load
#        external_ring_load = self.load_bearing_results[ind_load_case].external_ring_load
#        
#        delta_1 = (D - D1)/10.
#        delta_2 = (d1 - d)/10.
#        positions = {}
#        positions[list_node[0]] = vm.Point2D((-B/2., D/2. - delta_1))
#        positions[list_node[1]] = vm.Point2D((-B/2., D/2. - 3*delta_1))
#        positions[list_node[2]] = vm.Point2D((-B/2., d/2. + 3*delta_2))
#        positions[list_node[3]] = vm.Point2D((-B/2., d/2. + delta_2))
#        positions[list_node[4]] = vm.Point2D((B/2., D/2. - delta_1))
#        positions[list_node[5]] = vm.Point2D((B/2., D/2. - 3*delta_1))
#        positions[list_node[6]] = vm.Point2D((B/2., d/2. + 3*delta_2))
#        positions[list_node[7]] = vm.Point2D((B/2., d/2. + delta_2))
#        
#        if max_load is None:
#            max_load = 0
#            for nd in list_node:
#                if nd.load is not None:
#                    max_load = max(nd.load, max_load)
#
#        if (self.direction == 1) or (self.load_bearing_results[ind_load_case].sub_direction == 1):
#            if transversale_load != 0:
#                list_line = [vm.LineSegment2D(positions[list_node[3]], positions[list_node[5]])]
#                list_line.append(vm.LineSegment2D(positions[list_node[4]], positions[list_node[2]]))
#                load_arrow = vm.Contour2D(list_line)
#                load_arrow = load_arrow.Translation(vm.Vector2D((pos, 0)), True)
#                load_arrow.MPLPlot(a,'-b', True, width = B/10./max_load*transversale_load)
#        elif (self.direction == -1) or (self.load_bearing_results[ind_load_case].sub_direction == -1):
#            if transversale_load != 0:
#                list_line = [vm.LineSegment2D(positions[list_node[6]], positions[list_node[0]])]
#                list_line.append(vm.LineSegment2D(positions[list_node[1]], positions[list_node[7]]))
#                load_arrow = vm.Contour2D(list_line)
#                load_arrow = load_arrow.Translation(vm.Vector2D((pos, 0)), True)
#                load_arrow.MPLPlot(a,'-b', True, width = B/10./max_load*transversale_load)
#        if internal_ring_load != 0:
#            list_line = [vm.LineSegment2D(positions[list_node[3]], positions[list_node[7]])]
#            list_line.append(vm.LineSegment2D(positions[list_node[6]], positions[list_node[2]]))
#            load_arrow = vm.Contour2D(list_line)
#            load_arrow = load_arrow.Translation(vm.Vector2D((pos, 0)), True)
#            load_arrow.MPLPlot(a,'-b', True, width = B/10./max_load*internal_ring_load)
#        if external_ring_load != 0:
#            list_line = [vm.LineSegment2D(positions[list_node[1]], positions[list_node[5]])]
#            list_line.append(vm.LineSegment2D(positions[list_node[4]], positions[list_node[0]]))
#            load_arrow = vm.Contour2D(list_line)
#            load_arrow = load_arrow.Translation(vm.Vector2D((pos, 0)), True)
#            load_arrow.MPLPlot(a,'-b', True, width = B/10./max_load*external_ring_load)
#            

            
class RadialBearing:
    symmetric = None
    taking_loads = None
    generate_axial_load = None
    linkage = None
    
    _jsonschema = {
        "definitions": {},
        "$schema": "http://json-schema.org/draft-07/schema#",
        "type": "object",
        "title": "mechanical_components.bearings.RadialBearing Base Schema",
        "required": [
            "d",
            "D",
            "B",
            "alpha",
            "i",
            "Z",
            "Dw"
          ],
          "properties": {
            "d": {"type": "number",  "examples": [0.02],
                  "physical_quantity": 'distance', "unit": "m",
                  "editable": "true", "description": "Internal diameter"},
            "D": {"type": "number",  "examples": [0.02],
                  "physical_quantity": 'distance', "unit": "m",
                  "editable": "true", "description": "External diameter"},
            "B": {"type": "number",  "examples": [0.04],
                  "physical_quantity": 'distance', "unit": "m",
                  "editable": "true", "description": "Width"},
            "alpha": {"type": "number",  "examples": [0.],
                      "physical_quantity": 'angle', "unit": "rad",
                      "editable": "true", "description": "Contact angle"},
            "i": {"type": "number", "multipleOf": 1.0, "examples": [1],
                  "editable": "true", "description": "Row number"},
            "Z": {"type": "number", "multipleOf": 1.0, "examples": [13],
                  "editable": "true", "description": "Rolling number"},
            "Dw": {"type": "number",  "examples": [0.01],
                   "physical_quantity": 'distance', "unit": "m",
                   "editable": "true", "description": "Rolling diameter"},
            "Cr": {"type": "number",  "examples": [29500],
                   "physical_quantity": 'force', "unit": "N",
                   "editable": "true", "description": "Based dynamic load"},
            "C0r": {"type": "number",  "examples": [25000],
                    "physical_quantity": 'force', "unit": "N",
                    "editable": "true", "description": "Based static load"},
            "contact_type": {"examples": ["linear_contact"],
                     "physical_quantity": 'contact_property',
                     "anyOf": [{"type": "string", "enum": [ "linear_contact", "point_contact", "mixed_contact" ]}, {"type": "null"}],
                     "editable": "true", "description": "Contact type"},
            "mass": {"type": "number",  "examples": [0.04],
                     "physical_quantity": 'mass', "unit": "kg",
                     "editable": "true", "description": "Mass"},
            "name": {"type": "string",  "examples": 'SNR 6062',
                     "editable": "true", "description": "Name"},
            "cost": {"type": "number",  "examples": [4.0],
                     "physical_quantity": 'cost', "unit": "Euro",
                     "editable": "false", "description": "Cost"},
            "E": {"type": "number",  "examples": [0.04],
                  "physical_quantity": 'distance', "unit": "m",
                  "editable": "false", "description": "External rolling diameter"},
            "F": {"type": "number",  "examples": [0.04],
                  "physical_quantity": 'distance', "unit": "m",
                  "editable": "false", "description": "Internal rolling diameter"},
            "d1": {"type": "number",  "examples": [0.04],
                  "physical_quantity": 'distance', "unit": "m",
                  "editable": "false", "description": "External diameter of internal ring"},
            "D1": {"type": "number",  "examples": [0.04],
                  "physical_quantity": 'distance', "unit": "m",
                  "editable": "false", "description": "Internal diameter of external ring"},
            }
          }

    
    # TODO: remove oil
    def __init__(self, d, D, B, alpha, i, Z, Dw, Cr=None, C0r=None, 
                 material=material_iso, contact_type=None, mass=None, name='', metadata={}):

        self.d = d
        self.D = D
        self.B = B
        self.Dpw = (d + D)/2.
        self.i = i
        
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
        self.contact_type = contact_type
#        self.sub_direction = 0
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
        self.metadata = metadata
        if mass is None:
            self.mass = self.Mass()
        else:
            self.mass = mass
        self.cost = self.mass*self.cost_coefficient + self.cost_constant
            
        
    def __eq__(self, other_bearing):
        if self.class_name != other_bearing.class_name:
            return False
        
        for k,v in self.__dict__.items():
            if k in ['d', 'D', 'B', 'alpha', 'i', 'Z', 'Dw', 'Cr', 'C0r', 'mass']:
                v2 = getattr(other_bearing, k)
                if v != v2:
                    return False
        return True
    
    def __hash__(self):
        h = int(self.d*4e3) + int(self.D*12e3) + int(self.B*1e3)+self.i+ int(1000*self.mass)
        h += len(self.__class__.__name__)
        return h
        
    def Check(self):
        if self.d <= 0.:
            return False
        if self.d >= self.D:
            return False
        return True
    
    @classmethod
    def EstimateBaseLifeTime(cls, Fr, N, t, Cr):
        total_cycles = 0.
        Pr = 0.
        for fr, ni, ti in zip(Fr, N, t):
            C = fr**(cls.coeff_baselife)
            if C != 0.:
                cycles = ni * ti * 2 * math.pi
                Pr += cycles * C
                total_cycles += cycles
        if total_cycles == 0.:
            return [math.inf]
        else:
            Pr = (Pr / total_cycles) ** (1/cls.coeff_baselife)
            L10 = (Cr/Pr)**(cls.coeff_baselife)
            return L10
        
    def BaseLifeTime(self,Fr, Fa, N, t, Cr):
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
            C = self.EquivalentDynamicLoad(fr, fa)**(self.coeff_baselife)
            if C != 0.:
                cycles = ni * ti * 2 * math.pi
                Pr += cycles * C
                total_cycles += cycles

        if total_cycles != 0:
            Pr = (Pr / total_cycles) ** (1/self.coeff_baselife)
            L10 = (Cr/Pr)**(self.coeff_baselife)
            return L10
        else:
            raise BearingL10Error()
        
    def AdjustedLifeTime(self, Fr, Fa, N, t, T, Cr=None, C0r=None, S=0.9):
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
        
        for fr, fa, n, ti, Ti  in zip(Fr, Fa, N, t, T):
            if (((fr != 0.) or (fa != 0.)) and (n > 0.)):
                cycles = n * ti * 2 * math.pi
                total_cycles += cycles
                
                a1 = ((1-self.material.c_gamma)
                     * (math.log(1/S)/math.log(100/90.))**(1/self.material.weibull_e)
                     + self.material.c_gamma)
                L10 = self.BaseLifeTime([fr], [fa], [n], [ti], self.Cr)
                Pr = self.EquivalentDynamicLoad(fr, fa)
                # viscosité cinématique de référence
                if n < (1000*2*math.pi/60.):
                    nu1 = 45000*(n*60/(2*math.pi))**(-0.83)*(self.Dpw*1e3)**(-0.5)
                else:
                    nu1 = 4500*(n*60/(2*math.pi))**(-0.5)*(self.Dpw*1e3)**(-0.5)
                
                coeff_oil = self.oil.oil_kinematic_viscosity_curve
                nu = 10**(10**(coeff_oil['A']*math.log10(Ti)+coeff_oil['B']))-0.6
                kappa = nu/nu1
                # Oil Contamination
                ec = self.oil.OilParameterContamination(self.Dpw,3)
                # Wear limit load
                if self.Dpw<0.1:
                    Cu = self.C0r/8.2
                else:
                    Cu = self.C0r/8.2*(100/(self.Dpw*1e3))**0.3
                    
                kappa = min(kappa, 4)                
                a_iso = self.AIso(kappa, ec, Cu, Pr)       
                a_iso = min(50., a_iso)
    
                Lpi = a1*a_iso*L10 # Corrected lifetime
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
#                check_rules_iter, val_rules_iter = rules_snr.RuleAxialLoad(fr, fa, n, level_axial_load='constant_load')
#                val_rules = min(val_rules_iter, val_rules)
#                if check_rules_iter == False:
#                    check_rules = False
#        return check_rules, val_rules
    
    def Mass(self):
        # TODO: enhance this but without querying CAD volumes!
        return 7800 * math.pi*self.B*(self.D-self.d) * (self.d+self.D)

    def CADVolumes(self, center = vm.o3D, axis = vm.x3D):
        # TODO: mutualization of this in parent class?
        axis.Normalize()
        
        y = axis.RandomUnitNormalVector()
        z = axis.Cross(y)
        
        #Internal Ring
        IRC = self.InternalRingContour()        
        irc = primitives3D.RevolvedProfile(center, axis, z, [IRC], center,
                                         axis, angle=2*math.pi, name='Internal Ring')
        #External Ring
        ERC=self.ExternalRingContour()
        erc=primitives3D.RevolvedProfile(center, axis, z, [ERC], center,
                                         axis, angle=2*math.pi,name='External Ring')
        #roller
        ROL=self.RollingContourCAD()
        
        radius=self.F/2.+self.slack+self.Dw/2.
        rollers=[]
        theta=2*math.pi/self.Z
        
        for zi in range(int(self.Z)):
            center_roller = center + radius*math.cos(zi*theta) * y + radius*math.sin(zi*theta) * z
            rollers.append(primitives3D.RevolvedProfile(center_roller, axis, z, [ROL],
                                                    center_roller, axis,
                                                    angle=2*math.pi,name='Roller {}'.format(zi+1)))
        
        volumes = [irc, erc] + rollers
        return volumes   
    
    def FreeCADExport(self, fcstd_filepath, python_path='python', 
                      freecad_lib_path='/usr/lib/freecad/lib', export_types=['fcstd']):
        model = self.VolumeModel()
#        tolerance = self.D/130.
        model.FreeCADExport(fcstd_filepath, python_path=python_path,
                            freecad_lib_path=freecad_lib_path,
                            export_types=export_types)
        
    def PlotDataQuote(self, pos=0):
        delta_quote = 0.05*self.B
        plot_data = []
        #internal diameter
        quote_x = 1.1*self.B/2.
        line1 = vm.LineSegment2D(vm.Point2D((0, self.d/2.)), vm.Point2D((quote_x + delta_quote, self.d/2.)))
        line1.Translation(vm.Vector2D((pos, 0)))
        li_data = [line1.PlotData(color = (0,0,0), stroke_width = 0.1, dash = True)]
        line2 = vm.LineSegment2D(vm.Point2D((0, -self.d/2.)), vm.Point2D((quote_x + delta_quote, -self.d/2.)))
        line2.Translation(vm.Vector2D((pos, 0)))
        li_data.append(line2.PlotData(color = (0,0,0), stroke_width = 0.1, dash = True))
        line3 = vm.LineSegment2D(vm.Point2D((quote_x, self.d/2.)), vm.Point2D((quote_x, -self.d/2.)))
        line3.Translation(vm.Vector2D((pos, 0)))
        li_data.append(line3.PlotData(color = (0,0,0), stroke_width = 0.1, dash = False, marker = 'triangle_quote'))
        
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
        line1 = vm.LineSegment2D(vm.Point2D((0, self.D/2.)), vm.Point2D((quote_x - delta_quote, self.D/2.)))
        line1.Translation(vm.Vector2D((pos, 0)))
        li_data = [line1.PlotData(color = (0,0,0), stroke_width = 0.1, dash = True)]
        line2 = vm.LineSegment2D(vm.Point2D((0, -self.D/2.)), vm.Point2D((quote_x - delta_quote, -self.D/2.)))
        line2.Translation(vm.Vector2D((pos, 0)))
        li_data.append(line2.PlotData(color = (0,0,0), stroke_width = 0.1, dash = True))
        line3 = vm.LineSegment2D(vm.Point2D((quote_x, self.D/2.)), vm.Point2D((quote_x, -self.D/2.)))
        line3.Translation(vm.Vector2D((pos, 0)))
        li_data.append(line3.PlotData(color = (0,0,0), stroke_width = 0.1, dash = False, marker = 'triangle_quote'))
        
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
        line1 = vm.LineSegment2D(vm.Point2D((-self.B/2., -self.D/2.)), vm.Point2D((-self.B/2., -self.D/2. - quote_x - delta_quote)))
        line1.Translation(vm.Vector2D((pos, 0)))
        li_data = [line1.PlotData(color = (0,0,0), stroke_width = 0.1, dash = True)]
        line2 = vm.LineSegment2D(vm.Point2D((self.B/2., -self.D/2.)), vm.Point2D((self.B/2., -self.D/2. - quote_x - delta_quote)))
        line2.Translation(vm.Vector2D((pos, 0)))
        li_data.append(line2.PlotData(color = (0,0,0), stroke_width = 0.1, dash = True))
        line3 = vm.LineSegment2D(vm.Point2D((self.B/2., -self.D/2. - quote_x)), vm.Point2D((-self.B/2., -self.D/2. - quote_x)))
        line3.Translation(vm.Vector2D((pos, 0)))
        li_data.append(line3.PlotData(color = (0,0,0), stroke_width = 0.1, dash = False, marker = 'triangle_quote'))
        
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
        
    def Plot(self, direction=1, a=None, typ=None):
        bg = self.PlotContour(direction)
        if a is None:
            f, a = bg.MPLPlot(style = '-k')
        else:
            bg.MPLPlot(a,'-k')
            
        if typ == 'Graph':
            graph = self.PlotGraph()
            graph.MPLPlot(a, '--b', True)
            
        elif typ == 'Load':
            self.PlotLoad(a)

    def VolumeModel(self, center = vm.o3D, axis = vm.x3D):
        model=vm.VolumeModel([(self.name, self.CADVolumes(center, axis))])
        return model   
    
#    mass = property(Mass)
    
    def Dict(self, stringify_keys=True):
        """Export dictionary
        """
        d = {}
        d['d'] = self.d
        d['D'] = self.D
        d['B'] = self.B
        d['alpha'] = self.alpha
        d['i'] = self.i
        d['Z'] = self.Z
        d['Dw'] = self.Dw
        if hasattr(self, 'Cr'):
            d['Cr'] = self.Cr
        if hasattr(self, 'C0r'):
            d['C0r'] = self.C0r
        
        d['name'] = self.name
        d['metadata'] = self.metadata
        d['contact_type'] = self.contact_type
        d['class_name'] = self.class_name
        d['material'] = self.material.Dict()
        d['mass'] = self.mass
        
        if 'load_bearing_results' in d:
            load_bearing_results = []
            for bg_result in self.load_bearing_results:
                load_bearing_results.append(bg_result.Dict())
            d['load_bearing_results'] = load_bearing_results
        
        if stringify_keys:
            return StringifyDictKeys(d)
        
        return d
    
    @classmethod
    def DictToObject(cls, d):
        if d['class_name'] == 'RadialBallBearing':
            return RadialBallBearing.DictToObject(d)
        elif d['class_name'] == 'AngularBallBearing':
            return AngularBallBearing.DictToObject(d)
        elif d['class_name'] == 'SphericalBallBearing':
            return SphericalBallBearing.DictToObject(d)
        elif d['class_name'] == 'N':
            return N.DictToObject(d)
        elif d['class_name'] == 'NU':
            return NU.DictToObject(d)
        elif d['class_name'] == 'NUP':
            return NUP.DictToObject(d)
        elif d['class_name'] == 'NF':
            return NF.DictToObject(d)
        elif d['class_name'] == 'TaperedRollerBearing':
            return TaperedRollerBearing.DictToObject(d)
        else:
            print(d['class_name'])
            raise ValueError
        
    def to_shaft(self):
        return shafts_assembly.Shaft(self.PlotContour(), name=self.name)
    
    
# =============================================================
# Object générique roulement cylindrique
# =============================================================
#class ThrustBallBearings(persistent.Persistent):
#    #Butée axiale à bille
#    
#class ThrustRollerBearings(persistent.Persistent):
#    #Butée axiale à rouleaux
#    
#class ThrustNeedleRollerBearings(ThrustRollerBearings,persistent.Persistent):
#    #Butée axiale à rouleaux avec cage intégrée
#    
class RadialBallBearing(RadialBearing):
    symmetric = True
    taking_loads = 'both'
    generate_axial_load = False
    linkage = 'ball'
    coeff_baselife = 3.
    class_name = 'RadialBallBearing'
    cost_coefficient = 0.2
    cost_constant = 1
    
    _jsonschema = {
        "definitions": {
                "RadialBearing": RadialBearing._jsonschema},
        "allOf": [
            { "$ref": "#/definitions/RadialBearing" },
            {'type': 'object',
             "properties": {
                'alpha': {'const': 0},
              },
#            "required": ['alpha']
            }
          ]
      }
    
    def __init__(self, d, D, B, i=1, Z=None, Dw=None, Cr=None, C0r=None,
                 material=material_iso, contact_type=None, mass=None, name='', metadata={}):
        RadialBearing.__init__(self, d, D, B, alpha=0, i=i, Z=Z, Dw=Dw, Cr=Cr,
                               C0r=C0r, material=material,
                               contact_type=contact_type, mass=mass,
                               name=name, metadata=metadata)
        
        # estimation for the graph 2D description
        h1 = self.Dw/2. - (self.E - self.D1)/2.
        self.h = self.B/2. - self.Dw/2.*math.sin(math.acos(h1/(self.Dw/2.))) - 1e-4
        
    def EquivalentStaticLoad(self, fr, fa=None):
        #Charge radiale statique équivalente
        X0 = 0.6
        Y0 = 0.5
        P0r = max(fr,X0*fr+Y0*fa)
        return P0r
    
    @classmethod
    def EstimateEquivalentDynamicLoad(cls, fr):
        Pr = fr
        return Pr
    
    def EquivalentDynamicLoad(self, fr, fa=0):
        alphap = fsolve((lambda alphap:math.cos(5/180.*math.pi)/math.cos(alphap) \
                        -(1.+0.012534*(fa/(self.i*self.Z*((self.Dw*1e3)**2) \
                        *math.sin(alphap)))**(2/3.))),self.alpha + 1.)[0]
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
    
    def AIso(self, kappa, ec, Cu, Pr):
        if kappa < 0.4:
            f = lambda coeff:(1-(2.5671-2.2649/(kappa**0.054381))**(0.83)*((coeff)**(1/3.)))
        elif kappa < 1:
            f = lambda coeff:(1-(2.5671-1.9987/(kappa**0.019087))**(0.83)*((coeff)**(1/3.)))
        else:
            f = lambda coeff:(1-(2.5671-1.9987/(kappa**0.071739))**(0.83)*((coeff)**(1/3.)))
        coeff0 = fsolve(f,ec*Cu/Pr)[0]
        coeff = min(coeff0, ec*Cu/Pr)
        a_iso = 0.1*(f(coeff)**(-9.3))
        return a_iso
    
    def InternalRingContour(self):
        
        pbi2 = vm.Point2D((-self.B/2., self.d1/2.))
        pbi1 = pbi2.Translation(vm.Vector2D((self.h, 0)))
        pbi3 = vm.Point2D((-self.B/2., self.d/2.))
        pbi4 = vm.Point2D((self.B/2., self.d/2.))
        pbi5 = vm.Point2D((self.B/2., self.d1/2.))
        pbi6 = pbi5.Translation(vm.Vector2D((-self.h, 0)))
        bi1 = primitives2D.RoundedLineSegments2D([pbi6, pbi5, pbi4, pbi3, pbi2, pbi1],
                                                 {1: self.radius, 
                                                  2: self.radius,
                                                  3: self.radius,
                                                  4: self.radius},
                                                  False,
                                                  adapt_radius=True)
        cbi1 = vm.Arc2D(pbi1, vm.Point2D((0, self.F/2)), pbi6)
        return vm.Contour2D([bi1, cbi1])
        
    def ExternalRingContour(self):
        
        pbe2 = vm.Point2D((-self.B/2., self.D1/2.))
        pbe1 = pbe2.Translation(vm.Vector2D((self.h, 0)))
        pbe3 = vm.Point2D((-self.B/2., self.D/2.))
        pbe4 = vm.Point2D((self.B/2., self.D/2.))
        pbe5 = vm.Point2D((self.B/2., self.D1/2.))
        pbe6 = pbe5.Translation(vm.Vector2D((-self.h, 0)))
#        betest = primitives2D.RoundedLineSegments2D([pbe1, pbe2, pbe3, pbe4, pbe5, pbe6],{}, False)
#        betest.MPLPlot()

        be1 = primitives2D.RoundedLineSegments2D([pbe1, pbe2, pbe3, pbe4, pbe5, pbe6],
                                                 {1: self.radius, 
                                                  2: self.radius,
                                                  3: self.radius,
                                                  4: self.radius},
                                                  False,
                                                  adapt_radius=True)
        cbe1 = vm.Arc2D(pbe6, vm.Point2D((0, self.E/2)), pbe1)
        return vm.Contour2D([be1, cbe1])
    
    def RollingContour(self):
        
        p0 = vm.Point2D((0, 0))
        c1 = vm.Circle2D(p0, self.Dw/2.) 
        return vm.Contour2D([c1])
    
    def RollingContourCAD(self):
        p0 = vm.Point2D((-self.Dw/2., 0))
        p1 = vm.Point2D((0, self.Dw/2.))
        p2 = vm.Point2D((self.Dw/2., 0))
        a1 = vm.Arc2D(p0, p1, p2)
        l1 = vm.LineSegment2D(p2,p0)
#        c1 = vm.Circle2D(p0, self.Dw/2.) 
        return vm.Contour2D([a1, l1])
    
    def PlotContour(self, direction=1):
        
        be_sup = self.ExternalRingContour()
        bi_sup = self.InternalRingContour()
        ball_sup = self.RollingContour()
        ball_sup.Translation(vm.Vector2D((0, self.Dpw/2.)))
        
        bearing_sup = vm.Contour2D([be_sup, bi_sup, ball_sup])
        bearing_inf = bearing_sup.Rotation(vm.Point2D((0, 0)), math.pi, True)
        
        bg = vm.Contour2D([bearing_sup, bearing_inf])
        return bg
    
    def PlotData(self, pos=0, quote=True, constructor=True, direction=1):
        
        be_sup = self.ExternalRingContour()
        be_sup1 = be_sup.Translation(vm.Vector2D((pos, 0)), True)
        export_D3 = [be_sup1.PlotData('be_sup', fill = 'url(#diagonal-stripe-1)')]
        bi_sup = self.InternalRingContour()
        bi_sup1 = bi_sup.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(bi_sup1.PlotData('bi_sup', fill = 'url(#diagonal-stripe-1)'))
        ball_sup = self.RollingContour()
        ball_sup.Translation(vm.Vector2D((0, self.Dpw/2.)))
        ball_sup1 = ball_sup.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(ball_sup1.PlotData('ball_sup', fill = None))
        
        be_inf = be_sup.Rotation(vm.Point2D((0, 0)), math.pi, True)
        be_inf1 = be_inf.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(be_inf1.PlotData('be_inf', fill = 'url(#diagonal-stripe-1)'))
        bi_inf = bi_sup.Rotation(vm.Point2D((0, 0)), math.pi, True)
        bi_inf1 = bi_inf.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(bi_inf1.PlotData('bi_inf', fill = 'url(#diagonal-stripe-1)'))
        ball_inf = ball_sup.Rotation(vm.Point2D((0, 0)), math.pi, True)
        ball_inf1 = ball_inf.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(ball_inf1.PlotData('ball_inf', fill = None))
        
        if constructor:
            line1 = vm.LineSegment2D(vm.Point2D((-self.B/2., self.d/2.)), vm.Point2D((-self.B/2., -self.d/2.)))
            line1.Translation(vm.Vector2D((pos, 0)))
            li_data = [line1.PlotData(color = (0,0,0), stroke_width = 0.05, dash = False, marker = None)]
            line2 = vm.LineSegment2D(vm.Point2D((self.B/2., self.d/2.)), vm.Point2D((self.B/2., -self.d/2.)))
            line2.Translation(vm.Vector2D((pos, 0)))
            li_data.append(line2.PlotData(color = (0,0,0), stroke_width = 0.05, dash = False, marker = None))
            pt_data = {}
            pt_data['name'] = 'constructor line'
            pt_data['type'] = 'line'
            pt_data['plot_data'] = li_data
            export_D3.append(pt_data)
        
        if quote:
            export_D3.extend(self.PlotDataQuote(pos))
        
        return export_D3
    
    @classmethod
    def Graph(cls, list_node, direction=1):
        
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
            
    @classmethod
    def DictToObject(cls, d):
        if 'Cr' not in d.keys():
            d['Cr'] = None
        if 'C0r' not in d.keys():
            d['C0r'] = None
        obj = cls(d = d['d'], D = d['D'], B = d['B'], i = d['i'], Z = d['Z'], 
                  Dw = d['Dw'], Cr = d['Cr'], C0r = d['C0r'],
                  material = Material.DictToObject(d['material']),  
                  contact_type = d['contact_type'],
                  name=d['name'], metadata=d['metadata'], mass=d['mass'])
        return obj
    
    def Copy(self):
        if not hasattr(self, 'Cr'):
            Cr = None
        else:
            Cr = self.Cr
        if not hasattr(self, 'C0r'):
            C0r = None
        else:
            C0r = self.C0r
        obj = RadialBallBearing(d = self.d, D = self.D, B = self.B, i = self.i, Z = self.Z, 
                  Dw = self.Dw, Cr = Cr, C0r = C0r,
                  oil = self.oil, material = self.material,  
                  contact_type = self.contact_type,
                  name = self.name, metadata = self.metadata)
        return obj
        
class AngularBallBearing(RadialBearing):
    symmetric = False
    taking_loads = 'right'
    generate_axial_load = True
    linkage = 'ball_joint'
    coeff_baselife = 3.
    class_name = 'AngularBallBearing'
    cost_coefficient = 0.4
    cost_constant = 1.5
    
    _jsonschema = {
        "definitions": {
                "RadialBearing": RadialBearing._jsonschema},
        "allOf": [
            { "$ref": "#/definitions/RadialBearing" },
            {'type': 'object',
             "properties": {
                'i': {'const': 1},
              },
#            "required": ['i']
            }
          ]
      }
    
    def __init__(self, d, D, B, alpha, i=1, Z=None, Dw=None, Cr=None, C0r=None ,
                 material=material_iso, contact_type=None, mass=None, name='', metadata={}):
        RadialBearing.__init__(self, d, D, B, alpha=alpha, i=1, Z=Z, Dw=Dw, Cr=Cr,
                               C0r=C0r, material=material, 
                               contact_type=contact_type, mass=mass, name=name, metadata=metadata)

        
        # estimation for the graph 2D description
        h1 = self.Dw/2. - (self.E - self.D1)/2.
        self.h1 = self.B/2. - self.Dw/2.*math.sin(math.acos(h1/(self.Dw/2.))) - 1e-4
        h2 = 0.95*self.Dw/2.
        self.h2 = self.B/2. - self.Dw/2.*math.sin(math.acos(h2/(self.Dw/2.))) - 1e-4
        self.D2 = 0.6*(self.D - self.E) + self.E
        self.d2 = 0.6*(self.F - self.d) + self.d
        
    def EquivalentStaticLoad(self, fr, fa=None):
        #Charge radiale statique équivalente
        if self.i == 1:
            X0 = 0.5
            evol_Y0 = [0.52, 0.5, 0.46, 0.42, 0.38, 0.33, 0.29, 0.26, 0.22]
            evol_alpha = npy.array([5, 10, 15, 20, 25, 30, 35, 40, 45])/180*math.pi
            f = interpolate.interp1d(list(evol_alpha),evol_Y0, fill_value='extrapolate')
            Y0 = float(f(self.alpha))
        elif self.i == 2:
            X0 = 1
            evol_Y0 = [1.04, 1, 0.92, 0.84, 0.76, 0.66, 0.58, 0.52, 0.44]
            evol_alpha = npy.array([5, 10, 15, 20, 25, 30, 35, 40, 45])/180*math.pi
            f = interpolate.interp1d(list(evol_alpha),evol_Y0, fill_value='extrapolate')
            Y0 = float(f(self.alpha))
        P0r = max(fr,X0*fr+Y0*fa)
        return P0r
    
    @classmethod
    def EstimateEquivalentDynamicLoad(cls, fr):
        Pr = fr
        return Pr
    
    def EquivalentDynamicLoad(self, fr, fa = 0):
        alphap = fsolve((lambda alphap:math.cos(self.alpha)/math.cos(alphap) \
                    -(1+0.012534*(fa/(self.i*self.Z*((self.Dw*1e3)**2)*math.sin(alphap)))**(2/3.))),self.alpha + 1)[0]
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
            Y1 = fsolve((lambda Y1:Y1-(0.4/math.tan(alpha1))/(1-(1/3.)*math.sin(alpha1))),1)[0]
        if self.alpha > 15/180.*math.pi:
            e = (1-X1)/Y1
            Y3 = 0.625/e
        elif self.alpha <= 15/180.*math.pi:
            Y3 = 0.625/ksi*(1/math.tan(alphap))
        X2 = 1.625*X1
        Y2 = 1.625*Y1
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
    
    def AIso(self, kappa, ec, Cu, Pr):
        if kappa < 0.4:
            f = lambda coeff:(1-(2.5671-2.2649/(kappa**0.054381))**(0.83)*((coeff)**(1/3.)))
        elif kappa < 1:
            f = lambda coeff:(1-(2.5671-1.9987/(kappa**0.019087))**(0.83)*((coeff)**(1/3.)))
        else:
            f = lambda coeff:(1-(2.5671-1.9987/(kappa**0.071739))**(0.83)*((coeff)**(1/3.)))
        coeff0 = fsolve(f,ec*Cu/Pr)[0]
        coeff = min(coeff0, ec*Cu/Pr)
        a_iso = 0.1*(f(coeff)**(-9.3))
        return a_iso
    
    def InternalRingContour(self, direction=1, sign_V=1):
        
        pbi2 = vm.Point2D((direction*self.B/2., sign_V*self.d2/2.))
        pbi1 = vm.Point2D((direction*(self.B/2. - self.h2), sign_V*(self.Dpw/2. - self.Dw/2.*0.95)))
        pbi3 = vm.Point2D((direction*self.B/2., sign_V*self.d/2.))
        pbi4 = vm.Point2D((-direction*self.B/2., sign_V*self.d/2.))
        pbi5 = vm.Point2D((-direction*self.B/2., sign_V*self.d1/2.))
        pbi6 = pbi5.Translation(vm.Vector2D((direction*self.h1, 0)))
        bi1 = primitives2D.RoundedLineSegments2D([pbi6, pbi5, pbi4, pbi3, pbi2, pbi1], {1: self.radius, 
                                             2: self.radius, 3: self.radius, 4: self.radius}, False, adapt_radius = True)
        cbi1 = vm.Arc2D(pbi1, vm.Point2D((0, sign_V*self.F/2)), pbi6)
        irc = vm.Contour2D([bi1, cbi1])

        return irc
        
    def ExternalRingContour(self, direction=1, sign_V=1):
        
        pbe2 = vm.Point2D((direction*self.B/2., sign_V*self.D1/2.))
        pbe1 = pbe2.Translation(vm.Vector2D((-direction*self.h1, 0)))
        pbe3 = vm.Point2D((direction*self.B/2., sign_V*self.D/2.))
        pbe4 = vm.Point2D((-direction*self.B/2., sign_V*self.D/2.))
        pbe5 = vm.Point2D((-direction*self.B/2., sign_V*self.D2/2.))
        pbe6 = vm.Point2D((-direction*(self.B/2. - self.h2), sign_V*(self.Dpw/2. + self.Dw/2.*0.95)))
        be1 = primitives2D.RoundedLineSegments2D([pbe1, pbe2, pbe3, pbe4, pbe5, pbe6], {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius}, False, adapt_radius = True)
        cbe1 = vm.Arc2D(pbe6, vm.Point2D((0, sign_V*self.E/2)), pbe1)
        erc = vm.Contour2D([be1, cbe1])
        return erc

    
    def RollingContour(self):
        
        p0 = vm.Point2D((0, 0))
        c1 = vm.Circle2D(p0, self.Dw/2.) 
        return vm.Contour2D([c1])
    
    def RollingContourCAD(self):
        
        p0 = vm.Point2D((-self.Dw/2., 0))
        p1 = vm.Point2D((0, self.Dw/2.))
        p2 = vm.Point2D((self.Dw/2., 0))
        a1 = vm.Arc2D(p0, p1, p2)
        l1 = vm.LineSegment2D(p2,p0)
#        c1 = vm.Circle2D(p0, self.Dw/2.) 
        return vm.Contour2D([a1, l1])
    
    def PlotContour(self, direction=1):
        
        be_sup = self.ExternalRingContour(direction = direction, sign_V = 1)
        be_inf = self.ExternalRingContour(direction = direction, sign_V = -1)
        bi_sup = self.InternalRingContour(direction = direction, sign_V = 1)
        bi_inf = self.InternalRingContour(direction = direction, sign_V = -1)
        ball = self.RollingContour()
        ball_sup = ball.Translation(vm.Vector2D((0, self.Dpw/2.)), True)
        ball_inf = ball.Translation(vm.Vector2D((0, -self.Dpw/2.)), True)
        bg = vm.Contour2D([be_sup, bi_sup, ball_sup, be_inf, bi_inf, ball_inf])
        return bg
    
    def PlotData(self, pos=0, quote=True, constructor=True, direction=1):
        
        be_sup = self.ExternalRingContour(direction = direction, sign_V = 1)
        be_sup1 = be_sup.Translation(vm.Vector2D((pos, 0)), True)
        export_D3 = [be_sup1.PlotData('be_sup', fill = 'url(#diagonal-stripe-1)')]
        bi_sup = self.InternalRingContour(direction = direction, sign_V = 1)
        bi_sup1 = bi_sup.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(bi_sup1.PlotData('bi_sup', fill = 'url(#diagonal-stripe-1)'))
        ball = self.RollingContour()
        ball_sup = ball.Translation(vm.Vector2D((0, self.Dpw/2.)), True)
        ball_sup1 = ball_sup.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(ball_sup1.PlotData('ball_sup', fill = None))
        
        be_inf = self.ExternalRingContour(direction = direction, sign_V = -1)
        be_inf1 = be_inf.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(be_inf1.PlotData('be_inf', fill = 'url(#diagonal-stripe-1)'))
        bi_inf = self.InternalRingContour(direction = direction, sign_V = -1)
        bi_inf1 = bi_inf.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(bi_inf1.PlotData('bi_inf', fill = 'url(#diagonal-stripe-1)'))
        ball_inf = ball.Translation(vm.Vector2D((0, -self.Dpw/2.)), True)
        ball_inf1 = ball_inf.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(ball_inf1.PlotData('ball_inf', fill = None))
        
        if constructor:
            line1 = vm.LineSegment2D(vm.Point2D((-self.B/2., self.d/2.)), vm.Point2D((-self.B/2., -self.d/2.)))
            line1.Translation(vm.Vector2D((pos, 0)))
            li_data = [line1.PlotData(color = (0,0,0), stroke_width = 0.05, dash = False, marker = None)]
            line2 = vm.LineSegment2D(vm.Point2D((self.B/2., self.d/2.)), vm.Point2D((self.B/2., -self.d/2.)))
            line2.Translation(vm.Vector2D((pos, 0)))
            li_data.append(line2.PlotData(color = (0,0,0), stroke_width = 0.05, dash = False, marker = None))
            pt_data = {}
            pt_data['name'] = 'constructor line'
            pt_data['type'] = 'line'
            pt_data['plot_data'] = li_data
            export_D3.append(pt_data)
        
        if quote:
            export_D3.extend(self.PlotDataQuote(pos))
            
        return export_D3
    
    @classmethod
    def Graph(cls, list_node, direction=1):
        
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
            
    @classmethod
    def DictToObject(cls, d):
        if 'Cr' not in d.keys():
            d['Cr'] = None
        if 'C0r' not in d.keys():
            d['C0r'] = None
        obj = cls(d = d['d'], D = d['D'], B = d['B'], alpha = d['alpha'], i = d['i'], Z = d['Z'], 
                  Dw = d['Dw'], Cr = d['Cr'], C0r = d['C0r'],
                  material = Material.DictToObject(d['material']),  
                  contact_type = d['contact_type'],
                  name=d['name'], metadata=d['metadata'], mass=d['mass'])
        return obj
    
    def Copy(self):
        if not hasattr(self, 'Cr'):
            Cr = None
        else:
            Cr = self.Cr
        if not hasattr(self, 'C0r'):
            C0r = None
        else:
            C0r = self.C0r
        obj = AngularBallBearing(d = self.d, D = self.D, B = self.B, alpha = self.alpha, 
                            i = self.i, Z = self.Z, 
                            Dw = self.Dw, Cr = Cr, C0r = C0r,
                            oil = self.oil, material = self.material,  
                            contact_type = self.contact_type,
                            name = self.name, metadata = self.metadata)
        return obj
            
class SphericalBallBearing(RadialBearing):
    symmetric = True
    taking_loads = 'both'
    generate_axial_load = False
    linkage = 'ball_joint'
    coeff_baselife = 3.
    class_name = 'SphericalBallBearing'
    cost_coefficient = 0.4
    cost_constant = 2
    
    _jsonschema = {
        "definitions": {
                "RadialBearing": RadialBearing._jsonschema},
        "allOf": [
            { "$ref": "#/definitions/RadialBearing" },
            {'type': 'object',
             "properties": {
                'alpha': {'const': 0.},
              },
#            "required": ['alpha']
            }
          ]
      }
    
    def __init__(self, d, D, B, alpha=0, i=1, Z=None, Dw=None, Cr=None, C0r=None,
                 material=material_iso, contact_type=None, mass=None, name='',
                 metadata={}):
        RadialBearing.__init__(self, d, D, B, alpha, i, Z, Dw, Cr, C0r,
                               material, contact_type, mass, name, metadata=metadata)
        
        
    def EquivalentStaticLoad(self, fr, fa=None):
        #Charge radiale statique équivalente
        if self.i == 1:
            X0 = 0.5
            Y0 = 0.22/math.tan(self.alpha)
        elif self.i == 2:
            X0 = 1
            Y0 = 0.44/math.tan(self.alpha)
        P0r = max(fr,X0*fr+Y0*fa)
        return P0r
    
    @classmethod
    def EstimateEquivalentDynamicLoad(cls, fr):
        Pr = fr
        return Pr
    
    def EquivalentDynamicLoad(self, fr, fa = 0):
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
    
    def AIso(self, kappa, ec, Cu, Pr):
        if kappa < 0.4:
            f = lambda coeff:(1-(2.5671-2.2649/(kappa**0.054381))**(0.83)*((coeff)**(1/3.)))
        elif kappa < 1:
            f = lambda coeff:(1-(2.5671-1.9987/(kappa**0.019087))**(0.83)*((coeff)**(1/3.)))
        else:
            f = lambda coeff:(1-(2.5671-1.9987/(kappa**0.071739))**(0.83)*((coeff)**(1/3.)))
        coeff0 = fsolve(f,ec*Cu/Pr)[0]
        coeff = min(coeff0, ec*Cu/Pr)
        a_iso = 0.1*(f(coeff)**(-9.3))
        return a_iso
    
    @classmethod
    def Graph(cls, list_node, direction=1):
        
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
    
    @classmethod
    def DictToObject(cls, d):
        if 'Cr' not in d.keys():
            d['Cr'] = None
        if 'C0r' not in d.keys():
            d['C0r'] = None
        obj = cls(d = d['d'], D = d['D'], B = d['B'], alpha = d['alpha'], 
                  i = d['i'], Z = d['Z'], 
                  Dw = d['Dw'], Cr = d['Cr'], C0r = d['C0r'],
                  material = Material.DictToObject(d['material']),  
                  contact_type = d['contact_type'],
                  name=d['name'], metadata=d['metadata'], mass=d['mass'])
        return obj
    
    def Copy(self):
        if not hasattr(self, 'Cr'):
            Cr = None
        else:
            Cr = self.Cr
        if not hasattr(self, 'C0r'):
            C0r = None
        else:
            C0r = self.C0r
        obj = SphericalBallBearing(d = self.d, D = self.D, B = self.B, alpha = self.alpha, 
                            i = self.i, Z = self.Z, 
                            Dw = self.Dw, Cr = Cr, C0r = C0r,
                            oil = self.oil, material = self.material,  
                            contact_type = self.contact_type,
                            name = self.name, metadata = self.metadata)
        return obj
            
# TODO remove alpha?
class RadialRollerBearing(RadialBearing):
    """
    Abstract Class
    """
    symmetric = True
    linkage = 'cylindric'
    coeff_baselife = 10/3.
    cost_coefficient = 0.5
    cost_constant = 2.5
    
    _jsonschema = {
        "definitions": {
                "RadialBearing": RadialBearing._jsonschema},
        "allOf": [
            { "$ref": "#/definitions/RadialBearing" },
            {'type': 'object',
             "properties": {
                'i': {'const': 1},
                'contact_type': {'const': 'linear_contact'},
              },
            "required": ['contact_type']
            }
          ]
      }
    
    def __init__(self, d, D, B, alpha, i=1, Z=None, Dw=None, Cr=None, C0r=None,
                 material=material_iso,
                 contact_type='linear_contact',
                 mass=None, name='', metadata={}):
        RadialBearing.__init__(self, d, D, B, alpha=alpha, i=1, Z=Z, Dw=Dw, 
                               Cr=Cr, C0r=C0r,
                               material=material, contact_type=contact_type, 
                               mass=mass, name=name, metadata=metadata)
#        self.typ = typ
        
        # estimation for the graph 2D description
        self.Dpw = (self.d + self.D)/2.
        self.Lw = 0.7*self.B
        self.h = self.B/2. - self.Lw/2. - 1e-4
        
    def EquivalentStaticLoad(self, fr, fa=None):
        #Charge radiale statique équivalente
        if self.alpha != 0:
            x0 = 0.5*self.i
            y0 = 0.22*1/math.tan(self.alpha)*self.i
        else:
            x0 = 1
            y0 = 0
        P0r = max(fr,x0*fr+y0*fa)
        return P0r
    
    @classmethod
    def EstimateEquivalentDynamicLoad(cls, fr):
        Pr = fr
        return Pr
    
    def EquivalentDynamicLoad(self, fr, fa = 0):
        
        ksi = 1.5 #param of the ISO 1281
        e = ksi*math.tan(self.alpha)
        nu = 1-0.15*math.sin(self.alpha)
        w = self.material.weibull_e*self.coeff_baselife
  
        if self.contact_type == 'point_contact':
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
        elif self.contact_type == 'linear_contact':
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
        elif self.contact_type == 'mixed_contact':
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
                if fa <= e*fr:
                    Pr = fr
                else:
                    Pr = X1*fr+Y1*fa
            elif self.i == 2:
                if fa <= e*fr:
                    Pr = X3*fr+Y3*fa
                else:
                    Pr = X2*fr+Y2*fa
        else: # for alpha=0 the axial load is not include in the L10
            if self.i == 1:
                Pr = fr
            elif self.i == 2:
                Pr = X3*fr
        return Pr
    
    def AIso(self, kappa, ec, Cu, Pr):
        if kappa < 0.4:
            f = lambda coeff:(1-(1.5859-1.3993/(kappa**0.054381))*((coeff)**0.4))
        elif kappa < 1:
            f = lambda coeff:(1-(1.5859-1.2348/(kappa**0.19087))*((coeff)**0.4))
        else:
            f = lambda coeff:(1-(1.5859-1.2348/(kappa**0.071739))*((coeff)**0.4))
        try:
            coeff0 = fsolve(f,ec*Cu/Pr)[0]
        except FloatingPointError:# TODO check this for a better solution
            coeff0 = ec*Cu/Pr
        coeff = min(coeff0, ec*Cu/Pr)
        if f(coeff) > 0:
            a_iso = 0.1*(f(coeff)**(-9.185))
        else:
            return 100.
        return a_iso
    
    
    def RollingContour(self):
        
        p1 = vm.Point2D((-self.Lw/2.,-self.Dw/2.))
        p2 = vm.Point2D((-self.Lw/2.,self.Dw/2.))
        p3 = vm.Point2D((self.Lw/2.,self.Dw/2.))
        p4 = vm.Point2D((self.Lw/2.,-self.Dw/2.))
        rol = primitives2D.RoundedLineSegments2D([p1, p2, p3, p4], {0: self.radius, 
                                             1: self.radius, 2: self.radius, 3: self.radius}, True)
        return vm.Contour2D([rol])
    
    def RollingContourCAD(self):
        p1 = vm.Point2D((-self.Lw/2., 0))
        p2 = vm.Point2D((-self.Lw/2., self.Dw/2.))
        p3 = vm.Point2D((self.Lw/2., self.Dw/2.))
        p4 = vm.Point2D((self.Lw/2., 0))
        rol = primitives2D.RoundedLineSegments2D([p1, p2, p3, p4], {1: self.radius, 2: self.radius}, True)
        return vm.Contour2D([rol])
    
    def PlotData(self, pos=0, quote=True, constructor=True, direction=1):
        
        be_sup = self.ExternalRingContour(direction = direction, sign_V = 1)
        be_sup1 = be_sup.Translation(vm.Vector2D((pos, 0)), True)
        export_D3 = [be_sup1.PlotData('be_sup', fill = 'url(#diagonal-stripe-1)')]
                
        be_inf = self.ExternalRingContour(direction = direction, sign_V = -1)
        be_inf1 = be_inf.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(be_inf1.PlotData('be_inf', fill = 'url(#diagonal-stripe-1)'))
                
        bi_sup = self.InternalRingContour(direction = direction, sign_V = 1)
        bi_sup1 = bi_sup.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(bi_sup1.PlotData('bi_sup', fill = 'url(#diagonal-stripe-1)'))
                
        bi_inf = self.InternalRingContour(direction = direction, sign_V = -1)
        bi_inf1 = bi_inf.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(bi_inf1.PlotData('bi_inf', fill = 'url(#diagonal-stripe-1)'))
                
        roller = self.RollingContour()
        roller_sup = roller.Translation(vm.Vector2D((0, self.Dpw/2.)), True)
        roller_sup1 = roller_sup.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(roller_sup1.PlotData('roller_sup', fill = 'none'))
                
        roller_inf = roller.Translation(vm.Vector2D((0, -self.Dpw/2.)), True)
        roller_inf1 = roller_inf.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(roller_inf1.PlotData('roller_inf', fill = 'none'))
        
        if constructor:
            line1 = vm.LineSegment2D(vm.Point2D((-self.B/2., self.d/2.)), vm.Point2D((-self.B/2., -self.d/2.)))
            line1.Translation(vm.Vector2D((pos, 0)))
            li_data = [line1.PlotData(color = (0,0,0), stroke_width = 0.05, dash = False, marker = None)]
            line2 = vm.LineSegment2D(vm.Point2D((self.B/2., self.d/2.)), vm.Point2D((self.B/2., -self.d/2.)))
            line2.Translation(vm.Vector2D((pos, 0)))
            li_data.append(line2.PlotData(color = (0,0,0), stroke_width = 0.05, dash = False, marker = None))
            pt_data = {}
            pt_data['name'] = 'constructor line'
            pt_data['type'] = 'line'
            pt_data['plot_data'] = li_data
            export_D3.append(pt_data)
        
        if quote:
            export_D3.extend(self.PlotDataQuote(pos))
            
        return export_D3
    
    def PlotContour(self, direction=1):
        
        be_sup = self.ExternalRingContour(direction = direction, sign_V = 1)
        be_inf = self.ExternalRingContour(direction = direction, sign_V = -1)
        bi_sup = self.InternalRingContour(direction = direction, sign_V = 1)
        bi_inf = self.InternalRingContour(direction = direction, sign_V = -1)
        roller = self.RollingContour()
        roller_sup = roller.Translation(vm.Vector2D((0, self.Dpw/2.)), True)
        roller_inf = roller.Translation(vm.Vector2D((0, -self.Dpw/2.)), True)
        
        bg = vm.Contour2D([be_sup, bi_sup, roller_sup, be_inf, bi_inf, roller_inf])
        return bg
    
    @classmethod
    def DictToObject(cls, d):
        if 'Cr' not in d.keys():
            d['Cr'] = None
        if 'C0r' not in d.keys():
            d['C0r'] = None
        obj = cls(d = d['d'], D = d['D'], B = d['B'], i = d['i'], Z = d['Z'], 
                  Dw = d['Dw'], Cr = d['Cr'], C0r = d['C0r'],
                  material = Material.DictToObject(d['material']),  
                  contact_type = d['contact_type'],
                  name=d['name'], metadata=d['metadata'], mass=d['mass'])
        return obj
    
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
            return N(d = self.d, D = self.D, B = self.B, i = self.i, Z = self.Z, 
                     Dw = self.Dw, Cr = Cr, C0r = C0r,
                     oil = self.oil, material = self.material,  
                     contact_type = self.contact_type,
                     name = self.name, metadata = self.metadata)
        elif self.class_name == 'NU':
            return NU(d = self.d, D = self.D, B = self.B, i = self.i, Z = self.Z, 
                     Dw = self.Dw, Cr = Cr, C0r = C0r,
                     oil = self.oil, material = self.material,  
                     contact_type = self.contact_type,
                     name = self.name, metadata = self.metadata)
        elif self.class_name == 'NF':
            return NF(d = self.d, D = self.D, B = self.B, i = self.i, Z = self.Z, 
                     Dw = self.Dw, Cr = Cr, C0r = C0r,
                     oil = self.oil, material = self.material,  
                     contact_type = self.contact_type,
                     name = self.name, metadata = self.metadata)
        elif self.class_name == 'NUP':
            return NUP(d = self.d, D = self.D, B = self.B, i = self.i, Z = self.Z, 
                     Dw = self.Dw, Cr = Cr, C0r = C0r,
                     oil = self.oil, material = self.material,  
                     contact_type = self.contact_type,
                     name = self.name, metadata = self.metadata)
    
class NUP(RadialRollerBearing):
    symmetric = True
    taking_loads = 'both'
    generate_axial_load = False
    class_name = 'NUP'
    cost_coefficient = 0.5
    cost_constant = 2.7
    
    # TODO: remove alpha?
    def __init__(self, d, D, B, i=1, Z=None, Dw=None, Cr=None, C0r=None ,
                 material=material_iso, contact_type='linear_contact',
                 mass=None, name='', metadata={}):
        RadialRollerBearing.__init__(self, d, D, B, alpha=0, i = i, Z = Z, Dw = Dw, Cr=Cr,
                                     C0r=C0r,
                                     material=material, contact_type=contact_type,
                                     mass=mass, name=name, metadata=metadata)
        
    def InternalRingContour(self, direction=1, sign_V=1):


        d1 = self.d1
        pbi2 = vm.Point2D((-direction*self.B/2., sign_V*d1/2.))
        pbi1 = vm.Point2D((-direction*(self.B/2. - self.h), sign_V*d1/2.))
        pbi0 = vm.Point2D((-direction*(self.B/2. - self.h), sign_V*(self.F/2.)))
        pbi3 = vm.Point2D((-direction*self.B/2., sign_V*self.d/2.))
        pbi4 = vm.Point2D((direction*self.B/2., sign_V*self.d/2.))
        pbi5 = vm.Point2D((direction*self.B/2., sign_V*d1/2.))
        pbi6 = vm.Point2D((direction*(self.B/2. - self.h), sign_V*d1/2.))
        pbi7 = vm.Point2D((direction*(self.B/2. - self.h), sign_V*(self.F/2.)))
        irc = primitives2D.RoundedLineSegments2D([pbi0, pbi1, pbi2, pbi3, pbi4, pbi5, pbi6, pbi7],
                           {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius, 
                            5: self.radius, 6: self.radius}, True, adapt_radius = True)

        return vm.Contour2D([irc])


        
    def ExternalRingContour(self, direction=1, sign_V=1):
        
        D1 = self.D1
        pbe2 = vm.Point2D((-direction*self.B/2., sign_V*D1/2.))
        pbe1 = vm.Point2D((-direction*(self.B/2. - self.h), sign_V*D1/2.))
        pbe0 = vm.Point2D((-direction*(self.B/2. - self.h), sign_V*(self.E/2.)))
        pbe3 = vm.Point2D((-direction*self.B/2., sign_V*self.D/2.))
        pbe4 = vm.Point2D((direction*self.B/2., sign_V*self.D/2.))
        pbe5 = vm.Point2D((direction*self.B/2., sign_V*D1/2.))
        pbe6 = vm.Point2D((direction*(self.B/2. - self.h), sign_V*D1/2.))
        pbe7 = vm.Point2D((direction*(self.B/2. - self.h), sign_V*(self.E/2.)))
        be1 = primitives2D.RoundedLineSegments2D([pbe0, pbe1, pbe2, pbe3, pbe4, pbe5, pbe6, pbe7],
                           {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius, 
                            5: self.radius, 6: self.radius}, True)

        erc = vm.Contour2D([be1])
        return erc
    
    @classmethod
    def Graph(cls, list_node, direction=1):
        
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
    symmetric = True
    taking_loads = 'free'
    generate_axial_load = False
    class_name = 'N'
    cost_coefficient = 0.5
    cost_constant = 2.5
    
    def __init__(self, d, D, B, i=1, Z=None, Dw=None, Cr=None, C0r=None ,
                 material=material_iso, contact_type='linear_contact',
                 mass=None, name='', metadata={}):
        RadialRollerBearing.__init__(self, d, D, B, alpha=0, i = i, Z = Z, Dw = Dw, Cr=Cr,
                                     C0r=C0r,
                                     material=material, contact_type=contact_type,
                                     mass=mass, name=name, metadata=metadata)
        
    def InternalRingContour(self, direction=1, sign_V=1):

        d1 = self.d1
        pbi2 = vm.Point2D((-direction*self.B/2., sign_V*d1/2.))
        pbi1 = vm.Point2D((-direction*(self.B/2. - self.h), sign_V*d1/2.))
        pbi0 = vm.Point2D((-direction*(self.B/2. - self.h), sign_V*(self.F/2.)))
        pbi3 = vm.Point2D((-direction*self.B/2., sign_V*self.d/2.))
        pbi4 = vm.Point2D((direction*self.B/2., sign_V*self.d/2.))
        pbi5 = vm.Point2D((direction*self.B/2., sign_V*d1/2.))
        pbi6 = vm.Point2D((direction*(self.B/2. - self.h), sign_V*d1/2.))
        pbi7 = vm.Point2D((direction*(self.B/2. - self.h), sign_V*(self.F/2.)))
        irc = primitives2D.RoundedLineSegments2D([pbi0, pbi1, pbi2, pbi3, pbi4, pbi5, pbi6, pbi7],
                           {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius, 
                            5: self.radius, 6: self.radius}, True, adapt_radius = True)
        
        return vm.Contour2D([irc])


        
    def ExternalRingContour(self, direction=1, sign_V=1):
        
        D1 = self.E + 0.1*(self.D - self.E)
        pbe2 = vm.Point2D((-direction*self.B/2., sign_V*D1/2.))
        pbe1 = vm.Point2D((-direction*(self.B/2. - self.h), sign_V*(self.E/2.)))
        pbe3 = vm.Point2D((-direction*self.B/2., sign_V*self.D/2.))
        pbe4 = vm.Point2D((direction*self.B/2., sign_V*self.D/2.))
        pbe5 = vm.Point2D((direction*self.B/2., sign_V*D1/2.))
        pbe6 = vm.Point2D((direction*(self.B/2. - self.h), sign_V*(self.E/2.)))
        be1 = primitives2D.RoundedLineSegments2D([pbe1, pbe2, pbe3, pbe4, pbe5, pbe6], {1: self.radius, 
                           2: self.radius, 3: self.radius, 4: self.radius}, True)
        
        erc = vm.Contour2D([be1])
        return erc
    
    @classmethod
    def Graph(cls, list_node, direction=1):
        
        graph = nx.DiGraph()
        graph.add_edges_from([(list_node[4], list_node[0])])
        graph.add_edges_from([(list_node[1], list_node[5])])
        graph.add_edges_from([(list_node[6], list_node[2])])
        graph.add_edges_from([(list_node[3], list_node[7])])
            
        return graph
        
        
class NF(RadialRollerBearing):
    symmetric = True
    taking_loads = 'right'
    generate_axial_load = False
    class_name = 'NF'
    cost_coefficient = 0.5
    cost_constant = 2.5

    
    def __init__(self, d, D, B, i=1, Z=None, Dw=None, Cr=None, C0r=None,
                 material=material_iso, contact_type='linear_contact',
                 mass=None, name='', metadata={}):
        RadialRollerBearing.__init__(self, d, D, B, alpha=0, i = i, Z = Z, Dw = Dw, Cr=Cr,
                                     C0r=C0r,
                                     material=material, contact_type=contact_type,
                                     mass=mass, name=name, metadata=metadata)
        
    def InternalRingContour(self, direction=1, sign_V=1):

        d1 = self.d1
        pbi2 = vm.Point2D((-direction*self.B/2., sign_V*d1/2.))
        pbi1 = vm.Point2D((-direction*(self.B/2. - self.h), sign_V*d1/2.))
        pbi0 = vm.Point2D((-direction*(self.B/2. - self.h), sign_V*(self.F/2.)))
        pbi3 = vm.Point2D((-direction*self.B/2., sign_V*self.d/2.))
        pbi4 = vm.Point2D((direction*self.B/2., sign_V*self.d/2.))
        pbi5 = vm.Point2D((direction*self.B/2., sign_V*d1/2.))
        pbi6 = vm.Point2D((direction*(self.B/2. - self.h), sign_V*d1/2.))
        pbi7 = vm.Point2D((direction*(self.B/2. - self.h), sign_V*(self.F/2.)))
        irc = primitives2D.RoundedLineSegments2D([pbi0, pbi1, pbi2, pbi3, pbi4, pbi5, pbi6, pbi7],
                           {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius, 
                            5: self.radius, 6: self.radius}, True, adapt_radius = True)
        
        return vm.Contour2D([irc])


        
    def ExternalRingContour(self, direction=1, sign_V=1):
    
        D1 = self.D1
        D2 = self.E + 0.1*(self.D - self.E)
        pbe2 = vm.Point2D((-direction*self.B/2., sign_V*D1/2.))
        pbe1 = vm.Point2D((-direction*(self.B/2. - self.h), sign_V*D1/2.))
        pbe0 = vm.Point2D((-direction*(self.B/2. - self.h), sign_V*(self.E/2.)))
        pbe3 = vm.Point2D((-direction*self.B/2., sign_V*self.D/2.))
        pbe4 = vm.Point2D((direction*self.B/2., sign_V*self.D/2.))
        pbe5 = vm.Point2D((direction*self.B/2., sign_V*D2/2.))
        pbe6 = vm.Point2D((direction*(self.B/2. - self.h), sign_V*(self.E/2.)))
        be1 = primitives2D.RoundedLineSegments2D([pbe0, pbe1, pbe2, pbe3, pbe4, pbe5, pbe6],
                           {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius, 
                            5: self.radius}, True)

        erc = vm.Contour2D([be1])

        return erc
    
    @classmethod
    def Graph(cls, list_node, direction=1):
        
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
    symmetric = True
    taking_loads = 'free'
    generate_axial_load = False
    class_name = 'NU'
    cost_coefficient = 0.5
    cost_constant = 2.5

    def __init__(self, d, D, B, i=1, Z=None, Dw=None, Cr=None, C0r=None,
                 material=material_iso, contact_type='linear_contact',
                 mass=None, name='', metadata={}):
        RadialRollerBearing.__init__(self, d, D, B, alpha=0, i = i, Z = Z, Dw = Dw, Cr=Cr,
                                     C0r=C0r,
                                     material=material, contact_type=contact_type,
                                     mass=mass, name=name, metadata=metadata)

    def InternalRingContour(self, direction=1, sign_V=1):
        d1 = self.F - 0.1*(self.F - self.d)
        pbi2 = vm.Point2D((-direction*self.B/2., sign_V*d1/2.))
        pbi1 = vm.Point2D((-direction*(self.B/2. - self.h), sign_V*(self.F/2.)))
        pbi3 = vm.Point2D((-direction*self.B/2., sign_V*self.d/2.))
        pbi4 = vm.Point2D((direction*self.B/2., sign_V*self.d/2.))
        pbi5 = vm.Point2D((direction*self.B/2., sign_V*d1/2.))
        pbi6 = vm.Point2D((direction*(self.B/2. - self.h), sign_V*(self.F/2.)))
        irc = primitives2D.RoundedLineSegments2D([pbi1, pbi2, pbi3, pbi4, pbi5, pbi6], {1: self.radius, 
                           2: self.radius, 3: self.radius, 4: self.radius}, 
                           True, adapt_radius = True)
        
        return vm.Contour2D([irc])

        
    def ExternalRingContour(self, direction=1, sign_V=1):
        D1 = self.D1
        pbe2 = vm.Point2D((-direction*self.B/2., sign_V*D1/2.))
        pbe1 = vm.Point2D((-direction*(self.B/2. - self.h), sign_V*D1/2.))
        pbe0 = vm.Point2D((-direction*(self.B/2. - self.h), sign_V*(self.E/2.)))
        pbe3 = vm.Point2D((-direction*self.B/2., sign_V*self.D/2.))
        pbe4 = vm.Point2D((direction*self.B/2., sign_V*self.D/2.))
        pbe5 = vm.Point2D((direction*self.B/2., sign_V*D1/2.))
        pbe6 = vm.Point2D((direction*(self.B/2. - self.h), sign_V*D1/2.))
        pbe7 = vm.Point2D((direction*(self.B/2. - self.h), sign_V*(self.E/2.)))
        be1 = primitives2D.RoundedLineSegments2D([pbe0, pbe1, pbe2, pbe3, pbe4, pbe5, pbe6, pbe7],
                           {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius, 
                            5: self.radius, 6: self.radius}, True)
        erc = vm.Contour2D([be1])

        return erc
    
    @classmethod
    def Graph(cls, list_node, direction=1):
        
        graph = nx.DiGraph()
        graph.add_edges_from([(list_node[4], list_node[0])])
        graph.add_edges_from([(list_node[1], list_node[5])])
        graph.add_edges_from([(list_node[6], list_node[2])])
        graph.add_edges_from([(list_node[3], list_node[7])])
            
        return graph
    
class TaperedRollerBearing(RadialRollerBearing, AngularBallBearing):
    symmetric = False
    taking_loads = 'right'
    linkage = 'cylindric_joint'
    generate_axial_load = True
    coeff_baselife = 10/3.
    class_name = 'TaperedRollerBearing'
    cost_coefficient = 0.4
    cost_constant = 2
        
        
    _jsonschema = {
        "definitions": {
                "RadialBearing": RadialBearing._jsonschema},
        "allOf": [
            { "$ref": "#/definitions/RadialBearing" },
            {'type': 'object',
             "properties": {
                'i': {'const': 1},
                'contact_type': {'const': 'linear_contact'},
              },
            "required": ['contact_type']
            }
          ]
      }
    
    def __init__(self, d, D, B, alpha, i=1, Z=None, Dw=None, Cr=None, C0r=None,
                 material=material_iso, contact_type='linear_contact',
                 mass=None, name='', metadata={}):
        
        if Dw is None:
            self.Dw = (D - d)/7.*math.cos(alpha)

        RadialRollerBearing.__init__(self, d, D, B, alpha=alpha, i = i, Z = Z,
                                     Dw = Dw, Cr=Cr, C0r=C0r,
                                     material=material, contact_type=contact_type,
                                     mass=mass, name=name, metadata=metadata)
        
        # estimation for the graph 2D description
        self.Dpw = (self.d + self.D)/2.
        self.Lw = 0.7*self.B
        self.beta = math.atan(self.Dw/self.Dpw*math.sin(self.alpha))
        
    def InternalRingContour(self, direction=1, sign_V=1):
        
        shift_bi = 5e-4
#        shift_be = 1e-3

        p0 = vm.Point2D((0, sign_V*self.Dpw/2.))
        p1 = p0.Translation(vm.Vector2D((math.cos(self.alpha), -direction*sign_V*math.sin(self.alpha))), True)
        l1 = vm.Line2D(p0, p1)
        l1.Rotation(p0, direction*sign_V*self.beta)
        l1.Translation(vm.Vector2D((-0.8*direction*self.Dw/2.*math.sin(self.alpha), -sign_V*0.8*self.Dw/2.*math.cos(self.alpha))))
        l2 = l1.Translation(vm.Vector2D((-0.2*direction*self.Dw/2.*math.sin(self.alpha), -sign_V*0.2*self.Dw/2.*math.cos(self.alpha))), True)
        pbi3 = vm.Point2D((direction*(self.B/2. - shift_bi), sign_V*self.d/2.))
        pbi3T = pbi3.Translation(vm.Vector2D((0, 1)))
        pbi4 = vm.Point2D((-direction*(self.B/2.), sign_V*self.d/2.))
        pbi4T = pbi4.Translation(vm.Vector2D((0, 1)))
        l3 = vm.Line2D(pbi3, pbi3T)
        l4 = vm.Line2D(pbi4, pbi4T)
        pbi2 = vm.Point2D.LinesIntersection(l1, l3)
        pbi5 = vm.Point2D.LinesIntersection(l1, l4)
        l5 = vm.Line2D(vm.Point2D((direction*self.Lw/2.,0)), vm.Point2D((direction*self.Lw/2.,1)))
        l5.Rotation(vm.Point2D((0,sign_V*self.Dpw/2.)), -sign_V*direction*self.alpha)
        l6 = vm.Line2D(vm.Point2D((-direction*self.Lw/2.,0)), vm.Point2D((-direction*self.Lw/2.,1)))
        l6.Rotation(vm.Point2D((0,sign_V*self.Dpw/2.)),-sign_V*direction*self.alpha)
        pbi1 = vm.Point2D.LinesIntersection(l1, l5)
        pbi0 = vm.Point2D.LinesIntersection(l2, l5)
        pbi6 = vm.Point2D.LinesIntersection(l1, l6)
        pbi7 = vm.Point2D.LinesIntersection(l2, l6)
        
        bi1 = primitives2D.RoundedLineSegments2D([pbi0, pbi1, pbi2, pbi3, pbi4, pbi5, pbi6, pbi7], 
                                                 {1: self.radius, 2: self.radius, 3: self.radius,
                                                  4: self.radius, 5: self.radius,
                                                  6: self.radius}, True, adapt_radius=True)
        
        irc = vm.Contour2D([bi1])
        
        return irc
        
    def ExternalRingContour(self, direction=1, sign_V=1):
        
#        shift_bi = 5e-4
        shift_be = 1e-3

        p0 = vm.Point2D((0, sign_V*self.Dpw/2.))
        p1 = p0.Translation(vm.Vector2D((math.cos(self.alpha), -direction*sign_V*math.sin(self.alpha))), True)
        l0 = vm.Line2D(p0, p1)
        l0.Rotation(p0, -direction*sign_V*self.beta)
        l0.Translation(vm.Vector2D((direction*self.Dw/2.*math.sin(self.alpha), sign_V*self.Dw/2.*math.cos(self.alpha))))
        pbe3 = vm.Point2D((direction*self.B/2., sign_V*self.D/2.))
        pbe3T = pbe3.Translation(vm.Vector2D((0, 1)))
        pbe4 = vm.Point2D((-direction*(self.B/2. - shift_be), sign_V*self.D/2.))
        pbe4T = pbe4.Translation(vm.Vector2D((0, 1)))
        l3 = vm.Line2D(pbe3, pbe3T)
        l4 = vm.Line2D(pbe4, pbe4T)
        pbe2 = vm.Point2D.LinesIntersection(l0, l3)
        pbe5 = vm.Point2D.LinesIntersection(l0, l4)
        be1 = primitives2D.RoundedLineSegments2D([pbe2, pbe3, pbe4, pbe5], {0: self.radius, 
                                             1: self.radius, 2: self.radius, 3: self.radius}, True, adapt_radius=True)
        erc = vm.Contour2D([be1])
        
        return erc
    
    def CADVolumes(self, center = vm.o3D, axis = vm.x3D):
        axis.Normalize()
        
        y = axis.RandomUnitNormalVector()
#        y.vector = npy.round(y.vector,3)
#        y.vector = y.vector/y.Norm() 
        
#        z=vm.Vector3D(npy.cross(x.vector,y.vector))
#        z = axis.Cross(y)
        
        #Internal Ring
        IRC=self.InternalRingContour()        
        irc=primitives3D.RevolvedProfile(center, axis, y, [IRC], center,
                                         axis, angle=2*math.pi, name='Internal Ring')
        #External Ring
        ERC=self.ExternalRingContour()
        erc=primitives3D.RevolvedProfile(center, axis, y, [ERC], center,
                                         axis, angle=2*math.pi,name='External Ring')
        
        volumes = [irc, erc]
        return volumes


    
    def RollingContour(self, direction=1, sign_V=1):
        
        r1 = vm.Point2D((--direction*self.Lw/2., self.Dw/2. - self.Lw/2.*math.tan(self.beta)))
        r2 = vm.Point2D((-direction*self.Lw/2., self.Dw/2. + self.Lw/2.*math.tan(self.beta)))
        r3 = vm.Point2D((-direction*self.Lw/2., -self.Dw/2. - self.Lw/2.*math.tan(self.beta)))
        r4 = vm.Point2D((--direction*self.Lw/2., -self.Dw/2. + self.Lw/2.*math.tan(self.beta)))
        rol = primitives2D.RoundedLineSegments2D([r1, r2, r3, r4], {0: self.radius, 
                                             1: self.radius, 2: self.radius, 3: self.radius}, True)

        bg = vm.Contour2D([rol])
        return bg
    
    def PlotData(self, pos=0, direction=1, quote=True, constructor=True):
        
        be_sup = self.ExternalRingContour(direction = direction, sign_V = 1)
        be_sup1 = be_sup.Translation(vm.Vector2D((pos, 0)), True)
        export_D3 = [be_sup1.PlotData('be_sup', fill = 'url(#diagonal-stripe-1)')]
                
        be_inf = self.ExternalRingContour(direction = direction, sign_V = -1)
        be_inf1 = be_inf.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(be_inf1.PlotData('be_inf', fill = 'url(#diagonal-stripe-1)'))
                
        bi_sup = self.InternalRingContour(direction = direction, sign_V = 1)
        bi_sup1 = bi_sup.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(bi_sup1.PlotData('bi_sup', fill = 'url(#diagonal-stripe-1)'))
                
        bi_inf = self.InternalRingContour(direction = direction, sign_V = -1)
        bi_inf1 = bi_inf.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(bi_inf1.PlotData('bi_inf', fill = 'url(#diagonal-stripe-1)'))
                
        roller_sup = self.RollingContour(direction = direction, sign_V = 1)
        roller_sup = roller_sup.Rotation(vm.Point2D((0, 0)), -direction*self.alpha, True)
        roller_sup = roller_sup.Translation(vm.Vector2D((0, self.Dpw/2.)), True)
        roller_sup1 = roller_sup.Translation(vm.Vector2D((pos, 0)), True)
        export_D3.append(roller_sup1.PlotData('roller_sup', fill = 'none'))
        
        roller_inf = self.RollingContour(direction = direction, sign_V = -1)
        roller_inf = roller_inf.Rotation(vm.Point2D((0, 0)), direction*self.alpha, True)
        roller_inf = roller_inf.Translation(vm.Vector2D((0, -self.Dpw/2.)), True)
        roller_inf1 = roller_inf.Translation(vm.Vector2D((pos, 0)), True)
        
        if constructor:
            line1 = vm.LineSegment2D(vm.Point2D((-self.B/2., self.d/2.)), vm.Point2D((-self.B/2., -self.d/2.)))
            line1.Translation(vm.Vector2D((pos, 0)))
            li_data = [line1.PlotData(color = (0,0,0), stroke_width = 0.05, dash = False, marker = None)]
            line2 = vm.LineSegment2D(vm.Point2D((self.B/2., self.d/2.)), vm.Point2D((self.B/2., -self.d/2.)))
            line2.Translation(vm.Vector2D((pos, 0)))
            li_data.append(line2.PlotData(color = (0,0,0), stroke_width = 0.05, dash = False, marker = None))
            pt_data = {}
            pt_data['name'] = 'constructor line'
            pt_data['type'] = 'line'
            pt_data['plot_data'] = li_data
            export_D3.append(pt_data)
        
        export_D3.append(roller_inf1.PlotData('roller_inf', fill = 'none'))
        if quote:
            export_D3.extend(self.PlotDataQuote(pos))
            
        return export_D3
    
    def PlotContour(self, direction=1):
        
        be_sup = self.ExternalRingContour(direction = direction, sign_V = 1)
        be_inf = self.ExternalRingContour(direction = direction, sign_V = -1)
        bi_sup = self.InternalRingContour(direction = direction, sign_V = 1)
        bi_inf = self.InternalRingContour(direction = direction, sign_V = -1)
        roller_sup = self.RollingContour(direction = direction, sign_V = 1)
        roller_sup = roller_sup.Rotation(vm.Point2D((0, 0)), -direction*self.alpha, True)
        roller_sup = roller_sup.Translation(vm.Vector2D((0, self.Dpw/2.)), True)
        roller_inf = self.RollingContour(direction = direction, sign_V = -1)
        roller_inf = roller_inf.Rotation(vm.Point2D((0, 0)), direction*self.alpha, True)
        roller_inf = roller_inf.Translation(vm.Vector2D((0, -self.Dpw/2.)), True)
        
        bg = vm.Contour2D([be_sup, bi_sup, roller_sup, be_inf, bi_inf, roller_inf])
        return bg
    
    @classmethod
    def DictToObject(cls, d):
        if 'Cr' not in d.keys():
            d['Cr'] = None
        if 'C0r' not in d.keys():
            d['C0r'] = None
        obj = cls(d = d['d'], D = d['D'], B = d['B'], alpha = d['alpha'], i = d['i'], Z = d['Z'], 
                  Dw = d['Dw'], Cr = d['Cr'], C0r = d['C0r'],
                  material = Material.DictToObject(d['material']),  
                  contact_type = d['contact_type'],
                  name=d['name'], metadata=d['metadata'], mass=d['mass'])
        return obj
    
    def Copy(self):
        if not hasattr(self, 'Cr'):
            Cr = None
        else:
            Cr = self.Cr
        if not hasattr(self, 'C0r'):
            C0r = None
        else:
            C0r = self.C0r
        obj = TaperedRollerBearing(d = self.d, D = self.D, B = self.B, alpha = self.alpha, 
                            i = self.i, Z = self.Z, 
                            Dw = self.Dw, Cr = Cr, C0r = C0r,
                            oil = self.oil, material = self.material,  
                            contact_type = self.contact_type,
                            name = self.name, metadata = self.metadata)
        return obj
                 
class BearingCatalog:
    dessia_db_attributes = [{'name':'bearings',
                             'class':'mechanical_components.bearings.RadialBearing',
                             'type':'list'}]
    def __init__(self, bearings, name=''):
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
                
        self.name = name
        
    def __eq__(self, other_eb):
        
        equal = True
        for bg, other_bg in zip(self.bearings, other_eb.bearings):
            equal = equal and bg == other_bg
        return equal
    
    def __hash__(self):
        
        catalog_hash = 0.
        for bg in self.bearings:
            catalog_hash = catalog_hash + hash(bg)
        return int(catalog_hash%100000)
            
    @classmethod
    def LoadFromDataframe(cls, dataframe, catalog_name):
        bearings = []
        for index in dataframe.index:
            typ_rlt = dataframe.loc[index,'typ_bearing']
            d = round(dataframe.loc[index,'d'], 3)
            D = round(dataframe.loc[index,'D'], 3)
            B = round(dataframe.loc[index,'B'],3)
            i = dataframe.loc[index,'i']
            Z = dataframe.loc[index,'Z']
            Dw = dataframe.loc[index,'Dw']
            mass = None
            alp = dataframe.loc[index,'alpha']
            if str(alp) == 'nan':
                alpha = 0
            else:
                alpha = alp
            Cr = dataframe.loc[index,'Cr']*1e3
            C0r = dataframe.loc[index,'C0r']

            if typ_rlt == 'radial_roller_bearing':
                typ = dataframe.loc[index,'mounting']
                if typ == 'NUP':
                    bearings.append(NUP(d, D, B, i, 
                                                      Z, Dw, Cr, C0r, 
                                                      mass = mass))
                elif typ == 'N':
                    bearings.append(N(d, D, B, i, 
                                                      Z, Dw, Cr, C0r, 
                                                      mass = mass))
                elif typ == 'NF':
                    bearings.append(NF(d, D, B, i, 
                                                      Z, Dw, Cr, C0r, 
                                                      mass = mass,))
                elif typ == 'NU':
                    bearings.append(NU(d, D, B, i, 
                                                      Z, Dw, Cr, C0r, 
                                                      mass = mass,))
                    
            elif typ_rlt == 'radial_ball_bearing':
                bearings.append(RadialBallBearing(d, D, B, i, Z, 
                                                    Dw, Cr, C0r, mass = mass))
                
            elif typ_rlt == 'angular_ball_bearing':
                bearings.append(AngularBallBearing(d, D, B, i, Z, Dw, alpha,
                                                 Cr, C0r, mass = mass))
                
            elif typ_rlt == 'spherical_ball_bearing':
                bearings.append(SphericalBallBearing(d, D, B, i, Z, Dw, alpha, 
                                                   Cr, C0r, mass = mass))
                
            elif typ_rlt == 'tapered_roller_bearing':
                bearings.append(TaperedRollerBearing(d, D, B, i, Z, Dw, 
                                                   Cr, C0r, mass = mass))
        return cls(bearings, catalog_name)
                
    def Dict(self):
        d = {'name': self.name}
        bearings_dicts = []
        for bearing in self.bearings:
            bearings_dicts.append(bearing.Dict())
    
        d['bearings'] = bearings_dicts
        return d
    
    @classmethod
    def DictToObject(cls, dict_):
        bearings = [RadialBearing.DictToObject(b) for b in dict_['bearings']]
        return cls(bearings, dict_['name'])
        
    def SaveToFile(self, filepath, indent = 0):
        with open(filepath+'.json', 'w') as file:
            json.dump(self.Dict(), file, indent = indent)
        
    @classmethod
    def LoadFromFile(cls, filepath):
        if type(filepath) is str:            
            with open(filepath, 'r') as file:
                d = json.loads(file)
        else:
            d = json.loads(filepath.read().decode('utf-8'))
        return cls.DictToObject(d) 
    
    def SearchBearingCatalog(self, bearing_class, d, D):
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

    def NextBearingCatalog(self, bearing_class, d, D):
        try:
            next_bearings = self.bearings_by_dict[bearing_class][str(d)][str(D)]
            return next_bearings
        except KeyError:
            return False
    
    def Check(self):
        for bearing_class, bearings in self.bearings.items():
            for bearing in bearings:
                if not bearing.Check():
                        return False
        return True
    
    def InvalidBearings(self):
        invalid_bearings = []
        for bearing in self.bearings:
            if not bearing.Check():                    
                invalid_bearings.append(bearing)
        return invalid_bearings

    def Plot(self):
        d = [b.d for b in self.bearings]
        D = [b.D for b in self.bearings]
        B = [b.B for b in self.bearings]
        Cr = [b.Cr for b in self.bearings]
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
        ax1.plot(d ,D, 'o')
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
        
        
with pkg_resources.resource_stream(pkg_resources.Requirement('mechanical_components'),
                           'mechanical_components/catalogs/schaeffler.json') as schaeffler_json:
    schaeffler_catalog = BearingCatalog.LoadFromFile(schaeffler_json)
    

bearing_classes = [RadialBallBearing, AngularBallBearing,
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

class ConceptualBearingCombination:
    def __init__(self, bearing_classes, directions, mounting):
        
        self.bearing_classes = bearing_classes
        self.directions = directions
        self.mounting = mounting
        
    def __eq__(self, other_eb):
        
        equal = True
        for bg, other_bg in zip(self.bearing_classes, other_eb.bearing_classes):
            equal = equal and bg == other_bg
        equal = equal and self.directions == other_eb.directions
        equal = equal and self.mounting == other_eb.mounting
        return equal
    
    def __hash__(self):
        h = 3*len(self.mounting) + 127*sum(self.directions)
        for class_ in self.bearing_classes:
            h += len(class_.__name__)
        return h
        
    def BearingCombination(self, bearings):
        if self.mounting == 'both':
            connection_bi=['left', 'right']
            connection_be=['left', 'right']
        elif self.mounting == 'free':
            connection_bi=['left', 'right']
            connection_be=[]
        elif self.mounting == 'right':
            connection_bi=['left']
            connection_be=['right']
        elif self.mounting == 'left':
            connection_bi=['right']
            connection_be=['left']
        return BearingCombination(bearings, self.directions, radial_load_linkage = [True]*len(bearings),
                                  connection_bi = connection_bi, connection_be = connection_be,
                                  behavior_link = self.mounting)

    def CheckKinematic(self):
        
        list_node_output = list(range(8))
        list_node_bearings = []
        for i, bearing_classe in enumerate(self.bearing_classes):
            list_node_bearings.append(npy.arange(8*(i + 1), 8*(i + 2)))
    
        nx_graph = self.Graph(list_node_bearings, list_node_output)
        
        check = self.CheckViabilityAngularBearing(nx_graph, list_node_bearings,
                                                  list_node_output)
        return check
    
    def Graph(self, list_node_bearings, list_node_output):

        list_left_bearing = list_node_bearings[0]
        list_right_bearing = list_node_bearings[-1]
        nx_graph = self.bearing_classes[0].Graph(list_left_bearing, self.directions[0])
        
        bg0_load = self.bearing_classes[0].taking_loads
        if self.directions[0] == -1:
            if self.bearing_classes[0].taking_loads == 'left':
                bg0_load = 'right'
            elif self.bearing_classes[0].taking_loads == 'right':
                bg0_load = 'left'             
        if (self.mounting == 'left') or (self.mounting == 'both'):
            if bg0_load in ['left', 'both']:
                nx_graph.add_edges_from([(list_node_output[1], list_left_bearing[1])])
                nx_graph.add_edges_from([(list_left_bearing[0], list_node_output[0])])
        if (self.mounting == 'right') or (self.mounting == 'both'):
            if bg0_load in ['right', 'both']:
                nx_graph.add_edges_from([(list_node_output[3], list_left_bearing[3])])
                nx_graph.add_edges_from([(list_left_bearing[2], list_node_output[2])])
    
        for num_bg, (dir1, dir2, bg1, bg2) in enumerate(zip(self.directions[0: -1], self.directions[1:], \
                    self.bearing_classes[0: -1], self.bearing_classes[1:])):
            list_nd1 = list_node_bearings[num_bg]
            list_nd2 = list_node_bearings[num_bg + 1]
            nx_graph = nx.compose(bg2.Graph(list_nd2, dir2), nx_graph)
            
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
        if (self.mounting == 'left') or (self.mounting == 'both'):
            if bg_end_load in ['left', 'both']:
                nx_graph.add_edges_from([(list_node_output[6], list_right_bearing[6])])
                nx_graph.add_edges_from([(list_right_bearing[7], list_node_output[7])])
        if (self.mounting == 'right') or (self.mounting == 'both'):
            if bg_end_load in ['right', 'both']:
                nx_graph.add_edges_from([(list_node_output[4], list_right_bearing[4])])
                nx_graph.add_edges_from([(list_right_bearing[5], list_node_output[5])])
                
        return nx_graph
    
    def CheckViabilityAngularBearing(self, nx_graph, list_node_bearings, li_node_output):
        # axial load generate by angular_bearing
        node_axial_ring = []
        if self.mounting in ['left', 'both']:
            node_axial_ring.append(li_node_output[0])
            node_axial_ring.append(li_node_output[7])
        if self.mounting in ['right', 'both']:
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
        if self.mounting in ['right', 'both']:
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
        if self.mounting in ['left', 'both']:
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
        
        if self.mounting == 'right':
            if valid:
                valid = AnalyseConnection(list_node_bearings[0][4], li_node_output[2])
            if valid:
                valid = AnalyseConnection(list_node_bearings[-1][3], li_node_output[5])
                
        if self.mounting == 'left':
            if valid:
                valid = AnalyseConnection(list_node_bearings[0][6], li_node_output[0])
            if valid:
                valid = AnalyseConnection(list_node_bearings[-1][1], li_node_output[7])
        
        return valid

class BearingCombination:
    
    _jsonschema = {
        "$schema": "http://json-schema.org/draft-07/schema#",
        "definitions": {
                "RadialBearing": RadialBearing._jsonschema},
        'type': 'object',
        "properties": {
            'bearings': {
               "type": "array",
               "items": { "$ref": "#/definitions/RadialBearing"},
                         },
          },
        "required": ['bearings']
    }
    
    dessia_db_attributes = [{'name':'bearings',
                             'class':'mechanical_components.bearings.RadialBearing',
                             'type':'list'}]

    def __init__(self, bearings, directions, radial_load_linkage, internal_pre_load=0, 
                 connection_bi=['left', 'right'], connection_be=['left', 'right'], behavior_link='both'):
        self.bearings = bearings
        self.radial_load_linkage = radial_load_linkage
        self.internal_pre_load = internal_pre_load
        self.connection_be = connection_be
        self.connection_bi = connection_bi
        self.behavior_link = behavior_link
        self.mass = 0
        self.cost = 0
        self.directions = directions
        
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
            
    def __eq__(self, other_eb):
        
        equal = True
        for bg, other_bg in zip(self.bearings, other_eb.bearings):
            equal = equal and bg == other_bg
        equal = (equal and self.directions == other_eb.directions
                       and self.radial_load_linkage == other_eb.radial_load_linkage
                       and self.internal_pre_load == other_eb.internal_pre_load
                       and self.connection_be == other_eb.connection_be
                       and self.connection_bi == other_eb.connection_bi
                       and self.behavior_link == other_eb.behavior_link)
        if hasattr(self, 'axial_positions') and hasattr(other_eb, 'axial_positions'):
            equal = (equal and self.axial_positions == other_eb.axial_positions
                           and self.internal_diameters == other_eb.internal_diameters
                           and self.external_diameters == other_eb.external_diameters
                           and self.length == other_eb.length)
        elif hasattr(self, 'axial_positions') or hasattr(other_eb, 'axial_positions'):
            equal = False
        return equal
    
    def __hash__(self):
        h = 0
        for bg in self.bearings:
            h += hash(bg)
        h += sum(self.directions)
        return h
        
    def PlotGraph(self):
                
        for gp in self.graph:
            list_graph = []
            pos_m = -self.B/2.
            for bg in gp:
                graph = bg.PlotGraph()
                graph = graph.Translation((pos_m + bg.B/2., 0), True)
                pos_m += bg.B
                list_graph.append(graph)
            list_graph = vm.Contour2D(list_graph)
    #        list_graph = list_graph.Translation((pos, 0), True)
            f,a = list_graph.MPLPlot(style='--b',arrow= True)
            
    def Update(self, axial_positions, internal_diameters, external_diameters, length):
        # TODO Why axial position is not in init?
        self.axial_positions = axial_positions# TODO: move this in bearing assembly, for plots pass a center parameter
        self.internal_diameters = internal_diameters
        self.external_diameters = external_diameters
        self.length = length
    
    def ExternalBearing(self, sign=1):
        B = self.B
        d = self.d
        D = self.D
        ep = min(0.1*D, 0.1*d)
        De = D + ep
        Dg, Dd = D, D
        if 'left' in self.connection_be:
            Dg = Dg - ep
        if 'right' in self.connection_be:
            Dd = Dd -ep
        be = vm.Polygon2D([vm.Point2D((-B/2., sign*D/2.)), vm.Point2D((-B/2., sign*Dg/2.)),
                           vm.Point2D((-B/2. - ep, sign*Dg/2.)), vm.Point2D((-B/2. - ep, sign*De/2.)),
                           vm.Point2D((B/2. + ep, sign*De/2.)), vm.Point2D((B/2. + ep, sign*Dd/2.)),
                           vm.Point2D((B/2., sign*Dd/2.)), vm.Point2D((B/2., sign*D/2.))])
        return be
    
    def InternalBearing(self, sign=1):
        B = self.B
        d = self.d
        D = self.D
        ep = min(0.1*D, 0.1*d)
        di = d - ep
        dg, dd = d, d
        if 'left' in self.connection_bi:
            dg = dg + ep
        if 'right' in self.connection_bi:
            dd = dd + ep
        bi = vm.Polygon2D([vm.Point2D((-B/2., sign*d/2.)), vm.Point2D((-B/2., sign*dg/2.)),
                           vm.Point2D((-B/2. - ep, sign*dg/2.)), vm.Point2D((-B/2. - ep, sign*di/2.)),
                           vm.Point2D((B/2. + ep, sign*di/2.)), vm.Point2D((B/2. + ep, sign*dd/2.)),
                           vm.Point2D((B/2., sign*dd/2.)), vm.Point2D((B/2., sign*d/2.))])
        return bi
    
    def BearingBox(self, sign=1):
        box = vm.Polygon2D([vm.Point2D((self.axial_positions, sign*self.internal_diameters/2.)),
                      vm.Point2D((self.axial_positions, sign*self.external_diameters/2.)),
                      vm.Point2D((self.axial_positions + self.length, sign*self.external_diameters/2.)),
                      vm.Point2D((self.axial_positions + self.length, sign*self.internal_diameters/2.))])
        return box
    
    def PlotData(self, pos=0, box=False, typ=None, bearing_combination_result=None, quote=False, constructor=True):
        
        be_sup = vm.Contour2D([self.ExternalBearing(sign = 1)]).Translation(vm.Vector2D((pos, 0)), True)
        export_data = [be_sup.PlotData('be_sup', fill = 'url(#diagonal-stripe-1)')]
        be_inf = vm.Contour2D([self.ExternalBearing(sign = -1)]).Translation(vm.Vector2D((pos, 0)), True)
        export_data.append(be_inf.PlotData('be_inf', fill = 'url(#diagonal-stripe-1)'))
        bi_sup = vm.Contour2D([self.InternalBearing(sign = 1)]).Translation(vm.Vector2D((pos, 0)), True)
        export_data.append(bi_sup.PlotData('bi_sup', fill = 'url(#diagonal-stripe-1)'))
        bi_inf = vm.Contour2D([self.InternalBearing(sign = -1)]).Translation(vm.Vector2D((pos, 0)), True)
        export_data.append(bi_inf.PlotData('bi_inf', fill = 'url(#diagonal-stripe-1)'))
        
#        contour = []
        pos_m = -self.B/2.
        for bg, di in zip(self.bearings, self.directions):
            cont = bg.PlotData(pos = pos_m + bg.B/2. + pos, constructor = constructor, 
                               quote = False, direction = di)
#            cont1 = cont.Translation(vm.Vector2D((pos_m + bg.B/2. + pos, 0)), True)
#            cont_bg = vm.Contour2D([cont1])
            pos_m += bg.B
#            export = cont_bg.Plot3D()
            export_data.extend(cont)
            
        if typ == 'Load':
            pos_m = -self.B/2.
            for bg_ref, bg_simu in zip(self.bearings, bearing_combination_result):
                export_data.extend(bg_simu.PlotDataLoad(pos = pos + pos_m + bg_ref.B/2., d = bg_ref.d, D = bg_ref.D, 
                            B = bg_ref.B, d1 = bg_ref.d1, D1 = bg_ref.D1))
                pos_m += bg_ref.B
            
#        if typ == 'Load':
#            pos_m = -self.B/2.
#            for bg_ref, bg_simu in zip(self.bearings, bearing_combination_result):
#                bg_simu.PlotLoad(a, pos = pos + pos_m + bg_ref.B/2., d = bg_ref.d, D = bg_ref.D, 
#                            B = bg_ref.B, d1 = bg_ref.d1, D1 = bg_ref.D1)
#                pos_m += bg_ref.B
            
        if box:
            box_sup = vm.Contour2D([self.BearingBox(1)]).Translation(vm.Vector2D((pos, 0)), True)
            export_data.append(box_sup.PlotData('box_sup', fill = 'none', color='red', stroke_width = 0.3, opacity = 0.3))
            box_inf = vm.Contour2D([self.BearingBox(-1)]).Translation(vm.Vector2D((pos, 0)), True)
            export_data.append(box_inf.PlotData('box_inf', fill = 'none', color = 'red', stroke_width = 0.3, opacity = 0.3))
        
        return export_data
    
    def PlotContour2D(self, pos=0, a=None, box=True, typ='Graph'):
        be_sup = self.ExternalBearing(sign = 1)
        be_inf = self.ExternalBearing(sign = -1)
        bi_sup = self.InternalBearing(sign = 1)
        bi_inf = self.InternalBearing(sign = -1)
        contour = [be_sup, be_inf, bi_sup, bi_inf]
        linkage_area = vm.Contour2D(contour)
        linkage_area = linkage_area.Translation(vm.Vector2D((pos, 0)), True)
        
        contour = []
        pos_m = -self.B/2.
        for bg, di in zip(self.bearings, self.directions):
            cont = bg.PlotContour(direction = di)
            cont = cont.Translation(vm.Vector2D((pos_m + bg.B/2., 0)), True)
            pos_m += bg.B
            contour.append(cont)
        assembly_bg = vm.Contour2D(contour)
        assembly_bg = assembly_bg.Translation(vm.Vector2D((pos, 0)), True)
        
        return linkage_area, assembly_bg
    
    def SolveAxialLoad(self):
        shaft = unidimensional.Body(0, -1, name='Shaft')
        ground = unidimensional.Body(0, 0.5, name='Ground')
        
        p_shaft = unidimensional.Load(shaft, 500)
        id_ground = unidimensional.ImposedDisplacement(ground, 0.)
        imposed_displacements = [id_ground]
        loads = [p_shaft]
        
        
        component, nonlinear_linkages = self.ElementaryAxialLoad(ground, shaft, 0)
        bodies = [ground, shaft]
        for bir, bor in component:
            bodies.append(bir)
            bodies.append(bor)
        
        sm = unidimensional.UnidimensionalModel(bodies, [], nonlinear_linkages, loads,
                         imposed_displacements)
        result = sm.Solve(500)
        result.Plot(intensity_factor=1e-5)
    
    def ElementaryAxialLoad(self, ground, shaft, pos, radial_load, bearing_result, axial_load=None):
        nb_bg_radial = sum([1 if (p is True) else 0 for p in self.radial_load_linkage])
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
                    loads.append(unidimensional.Load(bor, -Fp))
                    loads.append(unidimensional.Load(bir, Fp))
                elif self.directions[num_bg] == 1:
                    global_axial_load += -Fp
                    link = unidimensional.CompressionSpring(bir, bor, k1, -j1, 'bearing {}'.format(num_bg))
                    nonlinear_linkages.append(link)
                    axial_bearings.append([link])
                    loads.append(unidimensional.Load(bor, Fp))
                    loads.append(unidimensional.Load(bir, -Fp))
            elif bg.taking_loads != 'free':
                if self.directions[num_bg] == -1:
                    link = unidimensional.CompressionSpring(bor, bir, k1, -j1, 'bearing {}'.format(num_bg))
                    nonlinear_linkages.append(link)
                    axial_bearings.append([link])
                elif self.directions[num_bg] == 1:
                    link = unidimensional.CompressionSpring(bir, bor, k1, -j1, 'bearing {}'.format(num_bg))
                    nonlinear_linkages.append(link)
                    axial_bearings.append([link])
            check_radial_linkage = self.radial_load_linkage[num_bg]
            if check_radial_linkage:
                bg_result.radial_load.append(radial_load/nb_bg_radial)
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
        if 'left' in self.connection_be:
            bor = component[0][1]
            nonlinear_linkages.append(unidimensional.UnilateralContact(ground, bor, bor.initial_position - ground.initial_position, name='Outer rings'))
        if 'right' in self.connection_be:
            bor = component[-1][1]
            nonlinear_linkages.append(unidimensional.UnilateralContact(bor, ground, ground.initial_position - bor.initial_position, name='Outer rings'))
        if 'left' in self.connection_bi:
            bir = component[0][0]
            nonlinear_linkages.append(unidimensional.UnilateralContact(shaft, bir, bir.initial_position - shaft.initial_position, name='Inner rings'))
        if 'right' in self.connection_bi:
            bir = component[-1][0]
            nonlinear_linkages.append(unidimensional.UnilateralContact(bir, shaft, shaft.initial_position - bir.initial_position, name='Inner rings'))
        return component, nonlinear_linkages, loads, axial_bearings, check_axial_load
    
    @classmethod
    def EstimateBaseLifeTime(cls, L10s):
        sum_L10_inv = 0
        for L10 in L10s:
            sum_L10_inv += (1/L10)**1.5
        return sum_L10_inv**(-1/1.5)
    
    def BaseLifeTime(self, bearing_combination_simulation_result):
        
        for bearing_result in bearing_combination_simulation_result.bearing_simulation_results:
            bearing_result.radial_load = []
            bearing_result.axial_load = []
                
        nb_bg_radial = sum([1 if (p is True) else 0 for p in self.radial_load_linkage])
        result_bgs = bearing_combination_simulation_result.bearing_simulation_results
        for radial_load, axial_load in zip(bearing_combination_simulation_result.radial_loads, 
                               bearing_combination_simulation_result.axial_loads):
            
            if (self.behavior_link != 'free') and (abs(axial_load) >= 1e-4):
                check_axial_load = self.AxialLoad(axial_load, radial_load, bearing_combination_simulation_result)
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
            L10 = bg.BaseLifeTime(Fr = bg_result.radial_load, Fa = bg_result.axial_load, 
                            N = speed, t = time, Cr = bg.Cr)
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
            
    def AxialLoad(self, axial_load, radial_load, bearing_combination_simulation_result):
        
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
            = self.ElementaryAxialLoad(ground, shaft, 0, radial_load, result_bgs, axial_load)
        loads = loads + loads_iter
        
        for bir, bor in component:
            bodies.append(bir)
            bodies.append(bor)
        nonlinear_linkages.extend(nonlinear_linkages_iter)
        
        sm = unidimensional.UnidimensionalModel(bodies, [], nonlinear_linkages, loads,
                         imposed_displacements)

        if check_axial_load: 
            result_sm = sm.Solve(500)
            bearing_combination_simulation_result.axial_load_model = result_sm
                    
            for num_bg, (axial_linkage, (bir, bor)) in enumerate(zip(axial_bearings, component)):
                for link in axial_linkage:
                    if link in result_sm.activated_nonlinear_linkages:
                        positions = (result_sm.positions[bir], result_sm.positions[bor])
                        result_bgs[num_bg].axial_load.append(link.Strains(positions))
            return True
        else:
            return False
    
    def Plot(self, pos=0, a=None, box=True, typ=None, ind_load_case=0):
        """
        Generate a Plot
        
        :param pos: axial position of bearing assembly center
        :param a: object define the append graph (default value is None)
        :param box: draw the box parameter of the bearing assembly
        :param typ: define the aditionnal draw (default is None), 'Graph' draw the graph connection between bearing, 'Load' define the load
        """
        linkage_area, assembly_bg = self.PlotContour2D(pos, a, box, typ)
        
        if a is None:
            f, a = linkage_area.MPLPlot(style = '-g')
        else:
            linkage_area.MPLPlot(a,'-g')
        
        
        assembly_bg.MPLPlot(a,'-k')

        if typ == 'Graph':
            list_graph = []
            pos_m = -self.B/2.
            for bg_ref, bg_simu in zip(self.bearings, bearing_combination_result):
                graph = bg_simu.PlotGraph(d = bg_ref.d, D = bg_ref.D, 
                                     B = bg_ref.B, d1 = bg_ref.d1, D1 = bg_ref.D1)
                graph = graph.Translation(vm.Vector2D((pos_m + bg_ref.B/2., 0)), True)
                pos_m += bg_ref.B
                list_graph.append(graph)
            list_graph = vm.Contour2D(list_graph)
            list_graph = list_graph.Translation(vm.Vector2D((pos, 0)), True)
            list_graph.MPLPlot(a, '--b', True)
            
        elif typ == 'Load':
            pos_m = -self.B/2.
            
            max_load = 0
            for bg in self.bearings:
                for nd in bg.load_bearing_results[ind_load_case].list_node:
                    if nd.load is not None:
                        max_load = max(nd.load, max_load)
                
            for bg in self.bearings:
                bg.PlotLoad(a, pos = pos + pos_m + bg.B/2., d = bg.d, D = bg.D, 
                            B = bg.B, d1 = bg.d1, D1 = bg.D1,
                            ind_load_case = ind_load_case,
                            max_load = max_load)
                pos_m += bg.B
        
        if box:
            box_sup = self.BearingBox(1)
            box_inf = self.BearingBox(-1)
            cont_box = [box_sup, box_inf]
            contour_box = vm.Contour2D(cont_box)
            contour_box = contour_box.Translation(vm.Vector2D((pos, 0)), True)
            contour_box.MPLPlot(a,'-r')
        
                
    def VolumeModel(self, center = vm.Point3D((0,0,0)), axis = vm.Vector3D((1,0,0))):
        groups = []
#        position = self.axial_positions
        center_bearing = center+0.5*(self.bearings[0].B -self.B)*axis
        for bearing in self.bearings:
            groups.append((bearing.name, bearing.CADVolumes(center=center_bearing)))
            center_bearing += bearing.B*axis
        model=vm.VolumeModel(groups)        
        return model   
    
    def Dict(self, subobjects_id={}, stringify_keys=True):
        """
        Export dictionary
        """
        d={}
        d['directions'] = self.directions
        
        d['radial_load_linkage'] = self.radial_load_linkage
        d['internal_pre_load'] = self.internal_pre_load 
        d['connection_bi'] = self.connection_bi
        d['connection_be'] = self.connection_be
        d['behavior_link'] = self.behavior_link

        bearings = []
        for bearing in self.bearings:
            if bearing in subobjects_id:
                bearings.append(subobjects_id[bearing])
            else:
                bearings.append(bearing.Dict())
        d['bearings'] = bearings

        if stringify_keys:
            return StringifyDictKeys(d)
        return d
    
    @classmethod
    def DictToObject(cls, d):
        bearings = []
        for bearing_s in d['bearings']:
            bearing = RadialBearing.DictToObject(bearing_s)
            bearings.append(bearing)
        obj = cls(bearings = bearings, directions = d['directions'], radial_load_linkage = d['radial_load_linkage'], 
                  internal_pre_load = 0, connection_bi = d['connection_bi'], 
                  connection_be = d['connection_be'], behavior_link = d['behavior_link'])
        return obj
        
class BearingAssembly:
    
    dessia_db_attributes = [{'name':'bearing_combinations',
                             'class':'mechanical_components.bearings.BearingCombination',
                             'type':'list'}]

    def __init__(self, bearing_combinations, pre_load=0, axial_positions=None):
        
        self.bearing_combinations = bearing_combinations
        self.mass = self.Mass()
        self.Cr_equ = self.Cr_equ()
        self.pre_load = pre_load
        self.load_bearing_assembly_results = None
        if axial_positions is not None:
            self.axial_positions = axial_positions
        self.B = 0
        for bc in bearing_combinations:
            self.B += bc.B
        self.D = 0
        for bc in bearing_combinations:
            self.D = max(bc.D, self.D)
        self.d = 0
        for bc in bearing_combinations:
            self.d = max(bc.d, self.d)
        self.cost = 0
        for bc in bearing_combinations:
            self.cost += bc.cost
            
    def __eq__(self, other_eb):
        equal = True
        for bc, other_bc in zip(self.bearing_combinations, other_eb.bearing_combinations):
            equal = equal and bc == other_bc
        equal = equal and self.pre_load == other_eb.pre_load
        if hasattr(self, 'axial_positions') and hasattr(other_eb, 'axial_positions'):
            equal = equal and list(self.axial_positions) == list(other_eb.axial_positions)
        elif hasattr(self, 'axial_positions') or hasattr(other_eb, 'axial_positions'):
            equal = False
        return equal
    
    def __hash__(self):
        h = 0
        for bc in self.bearing_combinations:
            h += hash(bc)
        return h
        
    def Update(self, axial_positions, internal_diameters, axial_pos, 
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
            assembly_bg.Update(pos, self.internal_diameters[num_linkage], self.external_diameters[num_linkage],
                               self.length[num_linkage])
            
    def Cr_equ(self):
        Cr_equ = 0
        for li_bg in self.bearing_combinations:
            for bg in li_bg.bearings:
                Cr_equ += (bg.Cr)
        return (Cr_equ)
        
    def Mass(self):
        mass = 0
        for li_bg in self.bearing_combinations:
            mass += li_bg.mass
        return mass
    
    def Shaft(self):
        d1 = self.bearing_combinations[0].d
        d2 = self.bearing_combinations[1].d
        B1 = self.bearing_combinations[0].B
        B2 = self.bearing_combinations[1].B
        pos_mid = (self.axial_positions[0] + self.axial_positions[1])/2.
        shaft = vm.Polygon2D([vm.Point2D((self.axial_positions[0] - B1/2. - 5e-3, -d1/2.)),
                      vm.Point2D((self.axial_positions[0] - B1/2. - 5e-3, d1/2.)),
                      vm.Point2D((pos_mid, d1/2.)),
                      vm.Point2D((pos_mid, d2/2.)),
                      vm.Point2D((self.axial_positions[1] + B2/2. + 5e-3, d2/2.)),
                      vm.Point2D((self.axial_positions[1] + B2/2. + 5e-3, -d2/2.)),
                      vm.Point2D((pos_mid, -d2/2.)),
                      vm.Point2D((pos_mid, -d1/2.)),])
        return shaft
    
    def PlotData(self, box=True, typ=None, constructor=False):
        
        shaft = self.Shaft()
        contour_shaft = vm.Contour2D([shaft])
        export_data = [contour_shaft.PlotData('contour_shaft', fill = 'none')]
        
        for assembly_bg, pos in zip(self.bearing_combinations, self.axial_positions):
            export_data.extend(assembly_bg.PlotData(pos, box, quote = False, constructor = constructor))
        return export_data
    
    def Plot(self, box=True, typ=None, ind_load_case=0):
        
        shaft = self.Shaft()
        contour_shaft = vm.Contour2D([shaft])
        f, a = contour_shaft.MPLPlot(style = '-k')
        
        for assembly_bg, pos in zip(self.bearing_combinations, self.axial_positions):
            assembly_bg.Plot(pos, a, box, typ, ind_load_case)
            
    def Graph(self):
        for li_bg in self.bearing_combinations:
            G, positions, li_axial_link, nd_axial_load = li_bg.Graph()
#            plt.figure()
#            nx.draw_networkx(G, pos = positions)
            
    def CheckLoad(self, bearing_assembly_simulation_result):
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
    def QuickShaftLoad(cls, positions, loads):              
    
        (pos1, pos2) = positions
        ground = genmechanics.Part('ground')
        shaft1 = genmechanics.Part('shaft1')
        p1 = npy.array([pos1, 0, 0])
        p2 = npy.array([pos2, 0, 0])
        bearing1 = linkages.FrictionlessBallLinkage(ground,shaft1,p1,[0,0,0],'bearing1')
        bearing2 = linkages.FrictionlessLinearAnnularLinkage(ground,shaft1,p2,[0,0,0],'bearing2')
        
        bearing_combination_results= []
        for ind_load_case, load_cases in enumerate(loads):
            load1 = []
            for pos, ld, tq in load_cases:
                load1.append(gm_loads.KnownLoad(shaft1, pos, [0,0,0], ld, tq, 'input'))
            load2 = gm_loads.SimpleUnknownLoad(shaft1, [(pos1 + pos2)/2,0,0], [0,0,0], [], [0], 'output torque')
            imposed_speeds = [(bearing1, 0, 100)]

            mech = genmechanics.Mechanism([bearing1, bearing2], ground, imposed_speeds, load1, [load2])
            
            axial_load = 0
            for pos, ld, tq in load_cases:
                axial_load += ld[0]
            bearing_combination_result = []
            for bg in [bearing1, bearing2]:
                tensor = mech.GlobalLinkageForces(bg,1)
                fr = (tensor[1]**2 + tensor[2]**2)**(0.5)
                bearing_combination_result.append(fr)
            bearing_combination_result.append(axial_load)
            bearing_combination_results.append(bearing_combination_result)
        return bearing_combination_results
                
    def ShaftLoad(self, positions, bearing_assembly_simulation_result):
        
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
        ground = genmechanics.Part('ground')
        shaft1 = genmechanics.Part('shaft1')
        p1 = npy.array([pos1, 0, 0])
        p2 = npy.array([pos2, 0, 0])
        bearing1 = linkages.FrictionlessBallLinkage(ground,shaft1,p1,[0,0,0],'bearing1')
        bearing2 = linkages.FrictionlessLinearAnnularLinkage(ground,shaft1,p2,[0,0,0],'bearing2')
        
        for ind_load_case, load_cases in enumerate(loads):
            load1 = []
            for pos, ld, tq in load_cases:
                load1.append(gm_loads.KnownLoad(shaft1, pos, [0,0,0], ld, tq, 'input'))
            load2 = gm_loads.SimpleUnknownLoad(shaft1, [(pos1 + pos2)/2,0,0], [0,0,0], [], [0], 'output torque')
            imposed_speeds = [(bearing1, 0, 100)]

            mech = genmechanics.Mechanism([bearing1, bearing2], ground, imposed_speeds, load1, [load2])
            
            axial_load = 0
            for pos, ld, tq in load_cases:
                axial_load += ld[0]
            for li_bg, bg, bearing_combination_result in zip(self.bearing_combinations, 
                                                                  [bearing1, bearing2], 
                                                                  result_bcs):
                tensor = mech.GlobalLinkageForces(bg,1)
                fr = (tensor[1]**2 + tensor[2]**2)**(0.5)
                bearing_combination_result.radial_loads.append(fr)
            self.AxialLoad(positions, axial_load, bearing_assembly_simulation_result)
        try:
            self.BaseLifeTime(bearing_assembly_simulation_result)
        except BearingL10Error:
            raise BearingL10Error()
        
    @classmethod
    def EstimateBaseLifeTime(cls, L10s):
        sum_L10_inv = 0
        for L10 in L10s:
            sum_L10_inv += (1/L10)**1.5
        return sum_L10_inv**(-1/1.5)
        
    def BaseLifeTime(self, bearing_assembly_simulation_result):
        
        result_bcs = bearing_assembly_simulation_result.bearing_combination_simulation_results
        result_bgs = [bg.bearing_simulation_results for bg in result_bcs]
        for bearing_combination, bearing_result in zip(self.bearing_combinations, result_bgs):
            for bg, bg_result in zip(bearing_combination.bearings, bearing_result):
                time = bearing_assembly_simulation_result.operating_times
                speed = bearing_assembly_simulation_result.speeds
                try:
                    L10 = bg.BaseLifeTime(Fr = bg_result.radial_load, Fa = bg_result.axial_load, 
                                    N = speed, t = time, Cr = bg.Cr)
#                if (str(L10) != 'nan') and (L10 != False):
                    bg_result.L10 = L10
                except BearingL10Error:
                    raise BearingL10Error()
                
        sum_L10_inv = 0
        valid_L10 = True
        for bearing_result in result_bgs:
            for bg_result in bearing_result:
#                if bg_result.L10 != False:
                sum_L10_inv += (1/bg_result.L10)**1.5
#                else:
#                    valid_L10 = False
                
        bearing_assembly_simulation_result.L10 = sum_L10_inv**(-1/1.5)
        
    def AxialLoad(self, positions, axial_load, bearing_assembly_simulation_result):
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
        for num_linkage, (bearing_combination, load_bearing_combination_result) in enumerate(zip(self.bearing_combinations, 
                                                      result_bcs)):
            
            pos = positions[num_linkage]
            radial_load = load_bearing_combination_result.radial_loads[-1]
            bearing_result = result_bgs[num_linkage]
            component, nonlinear_linkages_iter, loads_iter, axial_bearings, __\
                = bearing_combination.ElementaryAxialLoad(ground, shaft, pos, \
                                                          radial_load, bearing_result)
            if (bearing_combination.behavior_link != 'free') and (axial_load != 0):
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
            bearing_assembly_simulation_result.axial_load_model = result_sm
            for num_bc, (axial_linkages, component) in enumerate(zip(bc_axial_bearings, components)):
                if len(axial_linkages) == 0:
                    for result_bg in result_bgs[num_bc]:
                        result_bg.axial_load.append(0)
                for num_bg, (axial_linkage, (bir, bor)) in enumerate(zip(axial_linkages, component)):
                    for link in axial_linkage:
                        if link in result_sm.activated_nonlinear_linkages:
                            positions = (result_sm.positions[bir], result_sm.positions[bor])
                            result_bgs[num_bc][num_bg].axial_load.append(link.Strains(positions))

        except unidimensional.ModelConvergenceError:
            print('Convergence Error')
            pass
#            sm.PlotGraph()
#            sm.Plot()
#            raise unidimensional.ModelConvergenceError()
            
    
    def VolumeModel(self, center = (0,0,0), axis = (1,0,0)):
        groups = []
        
        for combination, combination_position in zip(self.bearing_combinations, self.axial_positions):
            position = combination_position
            
            for bearing in combination.bearings:
                groups.append((bearing.name, bearing.CADVolumes(center=vm.Point3D((position, 0, 0)))))
                position += bearing.B
        model=vm.VolumeModel(groups)        
        return model   
    
    def FreeCADExport(self, fcstd_filepath='An unamed bearing assembly', python_path='python', 
                      freecad_lib_path='/usr/lib/freecad/lib', export_types=['fcstd']):
        model = self.VolumeModel()
        model.FreeCADExport(fcstd_filepath, python_path=python_path,
                            freecad_lib_path=freecad_lib_path, export_types=export_types)
    
    def Dict(self, subobjects_id = {}, stringify_keys=True):
        """
        Export dictionary
        """
        d={}
        for k,v in self.__dict__.items():
            tv=type(v)
            if tv == npy.int64:
                d[k]=int(v)
            elif tv == npy.float64:
                d[k]=round(float(v), 5)
            else:
                d[k]=v
                
        li_bg = []
        for bearing_combination in self.bearing_combinations:
            if bearing_combination in subobjects_id:
                li_bg.append(subobjects_id[bearing_combination])
            else:
                li_bg.append(bearing_combination.Dict())
        d['bearing_combinations'] = li_bg
        d['axial_positions'] = list(self.axial_positions)
        
        if stringify_keys:
            return StringifyDictKeys(d)

        return d
        
    @classmethod
    def DictToObject(cls, d):
        li_bg = []
        for li in d['bearing_combinations']:
            BA = BearingCombination.DictToObject(li)
            li_bg.append(BA)
                 
        obj = cls(bearing_combinations = li_bg, pre_load = d['pre_load'], 
                  axial_positions = d['axial_positions'])
        if 'axial_pos' in d.keys():
            obj.Update(axial_positions = d['axial_positions'], 
                       internal_diameters = d['internal_diameters'],
                       axial_pos = d['axial_pos'], 
                       external_diameters = d['external_diameters'], 
                       length = d['length'])
        return obj

#class DetailedRadialRollerBearing(RadialRollerBearing):
#    #Roulement à rouleaux
#    def __init__(self, d, D, B, i, Z, Dw, d1=None, D1=None, Lw=None, radius=None, E=None, F=None, 
#                 alpha=0, bm=1.1, oil=oil_iso_vg_1500, material=material_iso, direction=1, typ='N'):
#        RadialRollerBearing.__init__(self, d, D, B, i, Z, Dw, alpha ,
#                                oil = oil, material = material, contact_type = 'linear_contact',
#                                direction = direction, typ = typ)
#
#        #diametre rouleau moyen
#        self.Dwe = Dw
#        self.bm = bm
#        
#        if d1 is not None:
#            self.d1 = d1
#        if D1 is not None:
#            self.D1 = D1
#        if Lw is not None:
#            self.Lw = Lw
#        if E is not None:
#            self.E = E
#        if F is not None:
#            self.F = F
#        if radius is not None:
#            self.radius = radius
#                
#        self.Dpw, self.Lwe, self.slack, self.ep = self.DefParam()
#        self.mass = self.Mass()
#        
#    def __str__(self):
#        s = '{}\n'.format(self.__class__.__name__)
#        s += 'Dimensions d:{} D:{} B:{} Dw:{} Z:{}'.format(self.d, self.D, self.B, self.Dw, self.Z)
##        s += 'L10: {} Lnm: {}
#        return s
#        
#    def Update(self, d1, D1, E, F, Z):
#        self.d1 = d1
#        self.D1 = D1
#        self.E = E
#        self.F = F
#        self.Z = Z
#        self.Dpw,self.Lwe,self.slack,self.ep = self.DefParam()
#        self.mass = self.Mass()
#        
#    def DefParam(self):
#        Dpw = (self.E+self.F)/2.
#        Lwe = self.Lw-2*self.radius
#        slack = (self.E-self.F-2*self.Dw)/4.
#        ep = (self.B-self.Lw-2*slack)/2.
#        return Dpw,Lwe,slack,ep
#        
#    def Mass(self):
#        rho = 7800
#        m = self.Z*math.pi*(self.Dw)**2/4.*self.Lw*rho
#        m += (math.pi*(self.D)**2/4.-math.pi*(self.E)**2/4.)*self.B*rho
#        m += (math.pi*(self.F)**2/4.-math.pi*(self.d)**2/4.)*self.B*rho
#        return m
#    
#    def BaseStaticLoad(self):
#        #Charge radiale statique de base
#        #besoin de convertir les dimensions en mm pour les formules ISO
#        C0r = 44*(1-(self.Dw*1e3*math.cos(self.alpha))/(self.Dpw*1e3))*self.i*self.Z \
#                *self.Lwe*1e3*self.Dw*1e3*math.cos(self.alpha)
#        return C0r
#    
#    def BaseDynamicLoad(self):
#        #Charge radiale dynamique de base
#        mu = float((self.Dwe*1e3)*math.cos(self.alpha)/(self.Dpw*1e3))
#        fc = 0.377*self.material.mu_delta*1/((2**((self.material.weibull_c\
#                +self.material.weibull_h-1)/(self.material.weibull_c-self.material.weibull_h+1)))\
#                *(0.5**(2*self.material.weibull_e/(self.material.weibull_c-self.material.weibull_h+1))))\
#                *self.material.B1*((1-mu)**((self.material.weibull_c+self.material.weibull_h-3)/\
#                (self.material.weibull_c-self.material.weibull_h+1))/((1+mu)**(2*self.material.weibull_e\
#                /(self.material.weibull_c-self.material.weibull_h+1))))\
#                *(mu**(2/(self.material.weibull_c-self.material.weibull_h+1)))\
#                *(1+(1.04*((1-mu)/(1+mu))**((self.material.weibull_c+self.material.weibull_h\
#                +2*self.material.weibull_e-3)/(self.material.weibull_c-self.material.weibull_h+1)))\
#                **((self.material.weibull_c-self.material.weibull_h+1)/2.))\
#                **(-2/(self.material.weibull_c-self.material.weibull_h+1))
#        Cr = fc*self.bm*self.i*((self.Lwe*1e3)*math.cos(self.alpha))\
#                **((self.material.weibull_c-self.material.weibull_h-1)/(self.material.weibull_c\
#                -self.material.weibull_h+1))*self.Z**((self.material.weibull_c\
#                -self.material.weibull_h-2*self.material.weibull_e+1)/(self.material.weibull_c\
#                -self.material.weibull_h+1))*(self.Dwe*1e3)**((self.material.weibull_c\
#                -self.material.weibull_h-3)/(self.material.weibull_c-self.material.weibull_h+1))
#        return Cr
#
#    def Dict(self):
#
#        d = {}
#        for k,v in self.__dict__.items():
#            tv = type(v)
#            if tv == npy.int64:
#                d[k] = int(v)
#            elif tv == npy.float64:
#                d[k] = round(float(v), 5)
#            else:
#                d[k] = v
#
#        d['oil'] = self.oil.Dict()
#        d['material'] = self.material.Dict()
#        return d

#class DetailedDrawnCupNeedleRollerBearing(DetailedRadialRollerBearing):
#    #Douille à aiguilles
#
#    def __init__(self, typ, B, d, D, d1, D1, Lw, Dw, radius, E, F, Z, i,
#                 alpha,bm=1, weibull_e=9/8., weibull_c=31/3., weibull_h=7/3.,
#                 B1=551.13373/0.483, mu_delta=0.83, c_gamma=0.05,
#                 oil_name='iso_vg_100'):
#        DetailedRadialRollerBearing.__init__(self, typ, B, d, D, d1, D1, Lw, Dw, radius, E,
#                                     F, Z, i, alpha, bm, weibull_e, weibull_c,
#                                     weibull_h, B1, mu_delta, c_gamma, oil_name)
#        
#class DetailedNeedleRollerBearing(DetailedRadialRollerBearing):
#    #Cage à aiguilles
#    def __init__(self, typ, B, d, D, d1, D1, Lw, Dw, radius, E, F, Z, i,
#                 alpha, bm=1, weibull_e=9/8., weibull_c=31/3., weibull_h=7/3.,
#                 B1=551.13373/0.483, mu_delta=0.83, c_gamma=0.05,
#                 oil_name='iso_vg_100'):
#        DetailedRadialRollerBearing.__init__(self, typ, B, d, D, d1, D1, Lw, Dw, radius, E,
#                                     F, Z, i, alpha, bm, weibull_e, weibull_c,
#                                     weibull_h, B1, mu_delta, c_gamma, oil_name)
#
#class DetailedSphericalRollerBearing(DetailedRadialRollerBearing):
#    #Roulement à rotule à rouleaux
#    def __init__(self, typ, B, d, D, d1, D1, Lw, Dw, radius, E, F, Z, i,
#                 alpha,bm=1.15, weibull_e=9/8., weibull_c=31/3., weibull_h=7/3.,
#                 B1=551.13373/0.483, mu_delta=0.83, c_gamma=0.05,
#                 oil_name='iso_vg_100'):
#        DetailedRadialRollerBearing.__init__(self, typ, B, d, D, d1, D1, Lw, Dw, radius, E,
#                                     F, Z, i, alpha, bm, weibull_e, weibull_c,
#                                     weibull_h, B1, mu_delta, c_gamma,
#                                     oil_name)
        
        
#class BearingAssemblyOptimizationResults:
#    
#    dessia_db_attributes = [{'name':'bearing_assemblies',
#                             'class':'mechanical_components.bearings.BearingAssembly',
#                             'type':'list'}]
#    
#    def __init__(self, loads, speeds, times,
#                 internal_diameters, axial_positions, external_diameters, length,
#                 linkage_types, mounting_types, bearing_numbers,
#                 sort, sort_arg, number_solutions, bearing_assemblies, results):
#        self.loads = loads
#        self.speeds = speeds
#        self.times = times
#        self.internal_diameters = internal_diameters
#        self.axial_positions = axial_positions
#        self.external_diameters = external_diameters
#        self.length = length
#        self.linkage_types = linkage_types
#        self.mounting_types = mounting_types
#        self.bearing_numbers = bearing_numbers
#        self.sort = sort
#        self.sort_arg = sort_arg
#        self.number_solutions = number_solutions
#        self.bearing_assemblies = bearing_assemblies
#        self.results = results
#
##    def PlotData(self):
##        plot_data = []
##        for ba in self.bearing_assemblies:
##            data = {}
##            data['bearing_assemblies'] = {'plot_data': [ba.PlotData()]}
##            data['bearing_assemblies']['bearing_combinations'] = {'plot_data': []}
##            data['bearing_assemblies']['bearing_combinations']['bearings'] = {'plot_data': []}
##            for num_bc, bc in enumerate(ba.bearing_combinations):
##                bc_result = self.results[ba][0]['bearing_combinations'][num_bc]
##                plot_data.append(bc.PlotData(typ='Load', bearing_combination_result = bc_result))
##                li_data = []
##                for bg in bc.bearings:
##                    li_data.append(bg.PlotData())
##                data['bearing_assemblies']['bearing_combinations']['bearings']['plot_data'].append(li_data)
##            plot_data.append(data)
#            
#            
##        return plot_data
#
#    def Dict(self, subobjects_id={}, stringify_keys=True):
#
#        """
#        Export dictionary
#        """
#        d={}
#        for k,v in self.__dict__.items():
#            tv=type(v)
#            if tv==npy.int64:
#                d[k]=int(v)
#            elif tv==npy.float64:
#                d[k]=round(float(v), 5)
#            else:
#                d[k]=v
#                
#        d['bearing_assemblies'] = []
#        for bearing_assembly in self.bearing_assemblies:
#            if bearing_assembly in subobjects_id:
#                d['bearing_assemblies'].append(subobjects_id[bearing_assembly])
#            else:                
#                d['bearing_assemblies'].append(bearing_assembly.Dict())
#                
#        dict_result = {}
#        for ba, results in self.results.items():
#            key_ba = self.bearing_assemblies.index(ba)
#            dict_result[int(key_ba)] = []
#            for result in results:
#                di_result_1 = {}
#                di_result_1['loads'] = result['loads']
#                di_result_1['axial_loads'] = result['axial_loads']
#                di_result_1['radial_loads'] = result['radial_loads']
#                di_result_1['positions'] = result['positions']
#                li_bcs = []
#                for num_bc, bc in enumerate(result['bearing_combinations']):
#                    li_bgs = []
#                    for bg in bc:
#                        li_bgs.append(bg.Dict())
#                    li_bcs.append(li_bgs)
#                di_result_1['bearing_combinations'] = li_bcs
#            dict_result[key_ba].append(di_result_1)
#        d['results'] = dict_result
#            
#        if 'max' in d['sort']:
#            del d['sort']['max']    
#        
#        if stringify_keys:
#            return StringifyDictKeys(d)
#        return d
#        
#    @classmethod
#    def DictToObject(cls, d):
#        
#        li_ba = []
#        for bearing_assembly in d['bearing_assemblies']:
#            li_ba.append(BearingAssembly.DictToObject(bearing_assembly))
#        res = {}
#        for key_ba, results in d['results'].items():
#            res[li_ba[int(key_ba)]] = []
#            for result in results:
#                dict_temp = {}
#                dict_temp['loads'] = result['loads']
#                dict_temp['axial_loads'] = result['axial_loads']
#                dict_temp['radial_loads'] = result['radial_loads']
#                dict_temp['positions'] = result['positions']
#                dict_temp['bearing_combinations'] = []
#                for num_bc, bc in enumerate(result['bearing_combinations']):
#                    bgs = []
#                    for bg in bc:
#                        bgs.append(LoadBearing.DictToObject(bg))
#                    dict_temp['bearing_combinations'].append(bgs)
#                res[li_ba[int(key_ba)]].append(dict_temp)
#        obj = cls(loads = d['loads'], 
#                  speeds = d['speeds'],
#                  times = d['times'],
#                  internal_diameters = d['internal_diameters'],
#                  axial_positions = d['axial_positions'],
#                  external_diameters = d['external_diameters'],
#                  length = d['length'],
#                  linkage_types = d['linkage_types'],
#                  mounting_types = d['mounting_types'],
#                  bearing_numbers = d['bearing_numbers'],
#                  sort = d['sort'],
#                  sort_arg = d['sort_arg'],
##                  path = d['path'],
#                 number_solutions = d['number_solutions'],
#                  bearing_assemblies = li_ba, results = res)
#                
#        return obj
#    
##    def DefOptimizer(self):
##        obj = BearingAssemblyOptimizer.DefOptimizer(self.list_pos_unknown, 
##                 self.loads, 
##                 self.torques, self.speeds, self.list_time,
##                 self.internal_diameters, self.axial_positions, self.external_diameters, self.length,
##                 self.linkage_type, self.mounting_type, self.bearing_numbers,
##                 self.sort, self.sort_arg, self.path, self.nb_sol)
##        return obj

class BearingSimulationResult:
    def __init__(self, axial_load=None, radial_load=None, L10=None):
        if axial_load is None:
            self.axial_load = []
        else:
            self.axial_load = axial_load
        if radial_load is None:
            self.radial_load = []
        else:
            self.radial_load = radial_load
        self.L10 = L10
        
    def __eq__(self, other_eb):
        equal = (self.axial_load == other_eb.axial_load
                 and self.radial_load == other_eb.radial_load
                 and self.L10 == other_eb.L10)
        return equal
    
    def __hash__(self):
        h = int(self.L10 % 230)
        return h
    
    def Dict(self, subobjects_id = {}, stringify_keys=True):
        """Export dictionary
        """
        d = {}
        for k,v in self.__dict__.items():
            tv = type(v)
            if tv == npy.int64:
                d[k] = int(v)
            elif tv==npy.float64:
                d[k]=round(float(v), 5)
            else:
                d[k] = v
                
        if stringify_keys:
            return StringifyDictKeys(d)

        return d
    
    @classmethod
    def DictToObject(cls, d):                 
        obj = cls(axial_load = d['axial_load'], 
                  radial_load = d['radial_load'], 
                  L10 = d['L10'])
        return obj

class BearingCombinationSimulationResult:
    dessia_db_attributes = [{'name':'bearing_simulation_results',
                             'class':'mechanical_components.bearings.BearingSimulationResult',
                             'type':'list'}]

    def __init__(self, bearing_simulation_results, axial_loads=None, radial_loads=None,
                 speeds=None, operating_times=None, axial_load_model=None):
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
        if axial_load_model is not None:
            self.axial_load_model = axial_load_model
            
    def __eq__(self, other_eb):
        equal = (self.axial_loads == other_eb.axial_loads
                 and self.radial_loads == other_eb.radial_loads)
        if hasattr(self, 'speeds') and hasattr(other_eb, 'speeds'):
            equal = equal and self.speeds == other_eb.speeds
        elif hasattr(self, 'speeds') or hasattr(other_eb, 'speeds'):
            equal = False
        if hasattr(self, 'operating_times') and hasattr(other_eb, 'operating_times'):
            equal = equal and self.operating_times == other_eb.operating_times
        elif hasattr(self, 'operating_times') or hasattr(other_eb, 'operating_times'):
            equal = False
        for bearing_simulation_result, other_bearing_simulation_result in zip(self.bearing_simulation_results, other_eb.bearing_simulation_results):
            equal = equal and bearing_simulation_result == other_bearing_simulation_result
        return equal
    
    def __hash__(self):
        h = 0
        for bsr in self.bearing_simulation_results:
            h += hash(bsr)
        return int(h % 2091)
            
    def Dict(self, subobjects_id = {}, stringify_keys=True):
        """
        Export dictionary
        """
        d={}
        for k,v in self.__dict__.items():
            tv=type(v)
            if tv==npy.int64:
                d[k]=int(v)
            elif tv==npy.float64:
                d[k]=round(float(v), 5)
            else:
                d[k]=v
        bgs = []
        for bg in self.bearing_simulation_results:
            if bg in subobjects_id:
                bgs.append(subobjects_id[bg])
            else:
                bgs.append(bg.Dict())
        d['bearing_simulation_results'] = bgs
        if stringify_keys:
            return StringifyDictKeys(d)

        return d
    
    @classmethod
    def DictToObject(cls, d):     
        bc_result = []
        for bg in d['bearing_simulation_results']:
            bc_result.append(BearingSimulationResult.DictToObject(bg))
        obj = cls(bearing_simulation_results = bc_result,
                  axial_loads = d['axial_loads'], 
                  radial_loads = d['radial_loads'])    
        return obj
    
    
class BearingAssemblySimulationResult:
    dessia_db_attributes = [{'name':'bearing_combination_simulation_results',
                             'class':'mechanical_components.bearings.BearingCombinationSimulationResult',
                             'type':'list'}]

    def __init__(self, bearing_combination_simulation_results, 
                 loads, speeds, operating_times, 
                 axial_load_model=None, L10=None):
        self.loads = loads
        self.speeds = speeds
        self.operating_times = operating_times
        self.axial_load_model = axial_load_model
        self.L10 = L10
        self.bearing_combination_simulation_results = bearing_combination_simulation_results
        
    def __eq__(self, other_eb):
        equal = (self.loads == other_eb.loads
                 and self.speeds == other_eb.speeds
                 and self.operating_times == other_eb.operating_times)
        if hasattr(self, 'L10') and hasattr(other_eb, 'L10'):
            equal = equal and self.L10 == other_eb.L10
        elif hasattr(self, 'L10') or hasattr(other_eb, 'L10'):
            equal = False
        for bearing_combination_simulation_result, other_bearing_combination_simulation_result \
                in zip(self.bearing_combination_simulation_results,\
                       other_eb.bearing_combination_simulation_results):
            equal = equal and bearing_combination_simulation_result == other_bearing_combination_simulation_result
        return equal
    
    def __hash__(self):
        h = int(self.L10 % 9397)
        return h
        
    def PlotAxialModel(self):
        self.axial_load_model.Plot(intensity_factor=1e-5)
        
    def Dict(self, subobjects_id = {}, stringify_keys=True):
        """
        Export dictionary
        """
        d = {}
        for k,v in self.__dict__.items():
            tv = type(v)
            if tv == npy.int64:
                d[k] = int(v)
            elif tv==npy.float64:
                d[k]=round(float(v), 5)
            else:
                d[k] = v
        del d['axial_load_model']
        bcs = []
        for bc in self.bearing_combination_simulation_results:
            if bc in subobjects_id:
                bcs.append(subobjects_id[bc])
            else:
                bcs.append(bc.Dict())
        d['bearing_combination_simulation_results'] = bcs
                
        if stringify_keys:
            return StringifyDictKeys(d)

        return d
    
    @classmethod
    def DictToObject(cls, d):
        li_bc = []
        for bc in d['bearing_combination_simulation_results']:
            li_bc.append(BearingCombinationSimulationResult.DictToObject(bc))
        obj = cls(bearing_combination_simulation_results = li_bc,
                  loads = d['loads'],
                  speeds = d['speeds'],
                  operating_times = d['operating_times'],
                  L10 = d['L10'])
        return obj

    
class BearingAssemblySimulation:
    dessia_db_attributes = [{'name':'bearing_assembly',
                             'class':'mechanical_components.bearings.BearingAssembly',
                             'type':'object'},
                            {'name':'bearing_assembly_simulation_result',
                             'class':'mechanical_components.bearings.BearingAssemblySimulationResult',
                             'type':'object'}]
    def __init__(self, bearing_assembly,
                 bearing_assembly_simulation_result):
        self.bearing_assembly = bearing_assembly
        self.bearing_assembly_simulation_result = bearing_assembly_simulation_result
        
    def __eq__(self, other_eb):
        equal = (self.bearing_assembly == other_eb.bearing_assembly
                 and self.bearing_assembly_simulation_result == other_eb.bearing_assembly_simulation_result)
        return equal
    
    def __hash__(self):
        br_hash = hash(self.bearing_assembly) \
                + hash(self.bearing_assembly_simulation_result)
        return br_hash
        
    def Dict(self, subobjects_id = {}, stringify_keys=True):
        """
        Export dictionary
        """
        d = {}
        for k,v in self.__dict__.items():
            tv = type(v)
            if tv == npy.int64:
                d[k] = int(v)
            elif tv==npy.float64:
                d[k]=round(float(v), 5)
            else:
                d[k] = v
                
        if self.bearing_assembly in subobjects_id:
            d['bearing_assembly'] = subobjects_id[self.bearing_assembly]
        else:
            d['bearing_assembly'] = self.bearing_assembly.Dict(subobjects_id)
        
        if self.bearing_assembly_simulation_result in subobjects_id:
            d['bearing_assembly_simulation_result'] = subobjects_id[self.bearing_assembly_simulation_result]
        else:
            d['bearing_assembly_simulation_result'] = self.bearing_assembly_simulation_result.Dict(subobjects_id)
                
        if stringify_keys:
            return StringifyDictKeys(d)

        return d
    
    @classmethod
    def DictToObject(cls, d):
        BA = BearingAssembly.DictToObject(d['bearing_assembly'])
        BAS = BearingAssemblySimulationResult.DictToObject(d['bearing_assembly_simulation_result'])
            
        obj = cls(bearing_assembly = BA, 
                  bearing_assembly_simulation_result = BAS)
        
        return obj

class BearingCombinationSimulation:
    dessia_db_attributes = [{'name':'bearing_combination',
                             'class':'mechanical_components.bearings.BearingCombination',
                             'type':'object'},
                            {'name':'bearing_combination_simulation_result',
                             'class':'mechanical_components.bearings.BearingCombinationSimulationResult',
                             'type':'object'}]
    def __init__(self, bearing_combination,
                 bearing_combination_simulation_result):
        self.bearing_combination = bearing_combination
        self.bearing_combination_simulation_result = bearing_combination_simulation_result
    
    def __eq__(self, other_eb):
        equal = (self.bearing_combination == other_eb.bearing_combination
                 and self.bearing_combination_simulation_result == other_eb.bearing_combination_simulation_result)
        return equal
    
    def __hash__(self):
        br_hash = hash(self.bearing_combination) \
                + hash(self.bearing_combination_simulation_result)
        return br_hash
        
    def Dict(self, subobjects_id = {}, stringify_keys=True):
        """
        Export dictionary
        """
        d = {}
        for k,v in self.__dict__.items():
            tv = type(v)
            if tv == npy.int64:
                d[k] = int(v)
            elif tv==npy.float64:
                d[k]=round(float(v), 5)
            else:
                d[k] = v
                
        if self.bearing_combination in subobjects_id:
            d['bearing_combination'] = subobjects_id[self.bearing_combination]
        else:
            d['bearing_combination'] = self.bearing_combination.Dict(subobjects_id)
        
        if self.bearing_combination_simulation_result in subobjects_id:
            d['bearing_combination_simulation_result'] = subobjects_id[self.bearing_combination_simulation_result]
        else:
            d['bearing_combination_simulation_result'] = self.bearing_combination_simulation_result.Dict(subobjects_id)
                
        if stringify_keys:
            return StringifyDictKeys(d)

        return d
    
    @classmethod
    def DictToObject(cls, d):
        BC = BearingCombination.DictToObject(d['bearing_combination'])
        BCS = BearingCombinationSimulationResult.DictToObject(d['bearing_combination_simulation_result'])
            
        obj = cls(bearing_combination = BC, 
                  bearing_combination_simulation_result = BCS)
        
        return obj
    

class BearingL10Error(Exception):
    def __init__(self):
        super().__init__('L10 simulation failed due to dynamic equivalent load equal 0')
        
class CatalogSearchError(Exception):
    def __init__(self):
        super().__init__('No available bearings with d and D defined')