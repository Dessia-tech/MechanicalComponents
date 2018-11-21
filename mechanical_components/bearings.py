import numpy as npy
#import math as mt
from scipy import interpolate
#import os
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math
from copy import deepcopy
import copy

from scipy.optimize import minimize,fsolve
import networkx as nx
import matplotlib.pyplot as plt
from itertools import product

import volmdlr as vm
import volmdlr.primitives2D as primitives2D
import volmdlr.primitives3D as primitives3D

from mechanical_components.bearings_snr import RadialRollerBearingSNR
from mechanical_components.catalogs.ISO_bearings \
    import dico_rlts_iso,dico_roller_iso,dico_radial_clearance_iso,dico_rules

import genmechanics
import genmechanics.linkages as linkages
import genmechanics.loads as loads

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
                        0.1:{npy.inf:{1:1,2:0.85,3:0.7,4:0.5,5:0.3,6:0.05,7:0}}}

class Oil:
    def __init__(self,oil_data,dict_oil_contamination):
        self.oil_data = oil_data
        self.oil_kinematic_viscosity_curve = self.KinematicViscosity(oil_data)
        self.dict_oil_contamination = dict_oil_contamination
    
    def Dict(self):
        return self.__dict__
    
    def FunCoeff(self,x,data,type_x='Linear',type_y='Linear'):
        if type_x == 'Log': 
            x = npy.log10(x)
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
                A = (npy.log10(npy.log10(0.6+val_np[0,1]))-npy.log10(npy.log10(0.6+val_np[-1,1])))/(npy.log10(val_np[0,0])-npy.log10(val_np[-1,0]))
                B = npy.log10(npy.log10(0.6+val_np[0,1]))-A*npy.log10(val_np[0,0])
                oil_kinematic_viscosity_curve[key]['A'] = A
                oil_kinematic_viscosity_curve[key]['B'] = B
        return oil_kinematic_viscosity_curve['data']
    
    def OilParameterContamination(self,Dpw,grade):
        for k,v in self.dict_oil_contamination.items():
            if (Dpw>=k) and (Dpw<list(v.keys())[0]):
                return list(v.values())[0][grade]
            
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
    def Dict2Obj(cls, d):
        obj = cls(oil_data = d['oil_data'],dict_oil_contamination = d['dict_oil_contamination'])
        return obj
    
oil_iso_vg_1500=Oil(iso_vg_1500,dict_oil_contamination)
oil_iso_vg_1000=Oil(iso_vg_1000,dict_oil_contamination)
oil_iso_vg_680=Oil(iso_vg_680,dict_oil_contamination)
oil_iso_vg_460=Oil(iso_vg_460,dict_oil_contamination)
oil_iso_vg_320=Oil(iso_vg_320,dict_oil_contamination)
oil_iso_vg_220=Oil(iso_vg_220,dict_oil_contamination)
oil_iso_vg_150=Oil(iso_vg_150,dict_oil_contamination)
oil_iso_vg_100=Oil(iso_vg_100,dict_oil_contamination)
oil_iso_vg_68=Oil(iso_vg_68,dict_oil_contamination)
oil_iso_vg_46=Oil(iso_vg_46,dict_oil_contamination)
oil_iso_vg_32=Oil(iso_vg_32,dict_oil_contamination)
oil_iso_vg_22=Oil(iso_vg_22,dict_oil_contamination)
oil_iso_vg_15=Oil(iso_vg_15,dict_oil_contamination)
oil_iso_vg_10=Oil(iso_vg_10,dict_oil_contamination)
        
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
    def Dict2Obj(cls, d):
        obj = cls(weibull_e = d['weibull_e'], weibull_c = d['weibull_c'], 
                  weibull_h = d['weibull_h'],
                  B1 = d['B1'], mu_delta = d['mu_delta'], c_gamma = d['c_gamma'])
        return obj

material_iso=Material()

class LoadNode:
    def __init__(self, num, ext_load=None, load=None):
        self.load = load
        self.ext_load = ext_load
        self.num = num
        if self.ext_load is None:
            self.load = None
        if self.ext_load is not None:
            self.load = 0
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

class LoadBearing:
    
    def CheckNone(self):
        input_node = [self.list_node[i] for i in [1, 3, 4, 6]]
        output_node = [self.list_node[i] for i in [0, 2, 5, 7]]
        for ind_input, ind_output in zip([1, 3, 4, 6], [0, 2, 5, 7]):
            if self.list_node[ind_input].load is None:
                self.list_node[ind_output].load = None
            elif self.list_node[ind_output].load is None:
                self.list_node[ind_input].load = None
    
    def PropagateNone(self):
        if (self.direction == 1 or self.sub_direction == 1):
            if self.list_node[3].load is None:
                self.list_node[6].load = None
            if self.list_node[4].load is None:
                self.list_node[1].load = None
            if (self.list_node[4].load is None) and (self.list_node[6].load is None):
                self.list_node[3].load = None
            if (self.list_node[1].load is None) and (self.list_node[3].load is None):
                self.list_node[4].load = None
        elif (self.direction == -1 or self.sub_direction == -1):
            if self.list_node[1].load is None:
                self.list_node[4].load = None
            if self.list_node[6].load is None:
                self.list_node[3].load = None
            if (self.list_node[1].load is None) and (self.list_node[3].load is None):
                self.list_node[6].load = None
            if (self.list_node[4].load is None) and (self.list_node[6].load is None):
                self.list_node[1].load = None
        elif self.direction == 0:
            if self.list_node[1].load is None:
                self.list_node[4].load = None
            if self.list_node[6].load is None:
                self.list_node[3].load = None
            if self.list_node[3].load is None:
                self.list_node[6].load = None
            if self.list_node[4].load is None:
                self.list_node[1].load = None
            
    def Load(self):
        
        valid = False
        if (self.direction == 1) or (self.sub_direction == 1):
            bi_input = [3, 6]
            bi_output = [2, 7]
            be_input = [4, 1]
            be_output = [5, 0]
            valid = True
        elif (self.direction == -1) or (self.sub_direction == -1):
            bi_input = [6, 3]
            bi_output = [7, 2]
            be_input = [1, 4]
            be_output = [0, 5]
            valid = True
        else:
            bi_axial_load = None
            be_axial_load = None
            bi_input = [3, 6]
            bi_output = [2, 7]
            be_input = [4, 1]
            be_output = [5, 0]
        if valid:
            if self.list_node[bi_input[0]].load is None:
                bi_axial_load = None
            elif self.list_node[bi_input[1]].load is not None:
                bi_axial_load = (self.list_node[bi_input[0]].ext_load + self.list_node[bi_input[0]].load)/2.
                bi_axial_load -= self.list_node[bi_input[1]].ext_load + self.list_node[bi_input[1]].load
            else:
                bi_axial_load = (self.list_node[bi_input[0]].ext_load + self.list_node[bi_input[0]].load)
            if self.list_node[be_input[0]].load is None:
                be_axial_load = None
            elif self.list_node[be_input[1]].load is not None:
                be_axial_load = (self.list_node[be_input[0]].ext_load + self.list_node[be_input[0]].load)/2.
                be_axial_load -= self.list_node[be_input[1]].ext_load + self.list_node[be_input[1]].load
            else:
                be_axial_load = (self.list_node[be_input[0]].ext_load + self.list_node[be_input[0]].load)
            if (bi_axial_load is not None) and (be_axial_load is not None):
                if (bi_axial_load < 0) and (be_axial_load < 0):
                    for n_input, n_output in zip(be_input, be_output[::-1]):
                        self.list_node[n_output].load = self.list_node[n_input].ext_load + self.list_node[n_input].load
                    for n_input, n_output in zip(bi_input, bi_output[::-1]):
                        self.list_node[n_output].load = self.list_node[n_input].ext_load + self.list_node[n_input].load
                        
                elif (bi_axial_load < 0) and (be_axial_load >= 0):
                    self.list_node[bi_output[0]].load = be_axial_load
                    if self.list_node[bi_input[1]].load is not None:
                        self.list_node[bi_output[0]].load += self.list_node[bi_input[1]].load + self.list_node[bi_input[1]].ext_load
                    if self.list_node[be_output[1]].load is not None:
                        self.list_node[be_output[1]].load = self.list_node[be_input[0]].load + self.list_node[be_input[0]].ext_load - be_axial_load
                    if self.list_node[bi_output[1]].load is not None:
                        self.list_node[bi_output[1]].load = self.list_node[bi_input[0]].load + self.list_node[bi_input[0]].ext_load
                elif (bi_axial_load >= 0) and (be_axial_load < 0):
                    self.list_node[be_output[0]].load = bi_axial_load
                    if self.list_node[be_input[1]].load is not None:
                        self.list_node[be_output[0]].load += self.list_node[be_input[1]].load + self.list_node[be_input[1]].ext_load
                    if self.list_node[bi_output[1]].load is not None:
                        self.list_node[bi_output[1]].load = self.list_node[bi_input[0]].load + self.list_node[bi_input[0]].ext_load - bi_axial_load
                    if self.list_node[be_output[1]].load is not None:
                        self.list_node[be_output[1]].load = self.list_node[be_input[0]].load + self.list_node[be_input[0]].ext_load
                elif (bi_axial_load >= 0) and (be_axial_load >= 0):
                    self.list_node[bi_output[0]].load = be_axial_load
                    if self.list_node[bi_input[1]].load is not None:
                        self.list_node[bi_output[0]].load += self.list_node[bi_input[1]].load + self.list_node[bi_input[1]].ext_load
                    self.list_node[be_output[0]].load = bi_axial_load
                    if self.list_node[be_input[1]].load is not None:
                        self.list_node[be_output[0]].load += self.list_node[be_input[1]].load + self.list_node[be_input[1]].ext_load
                    if self.list_node[be_output[1]].load is not None:
                        self.list_node[be_output[1]].load = self.list_node[be_input[0]].load + self.list_node[be_input[0]].ext_load - be_axial_load
                    if self.list_node[bi_output[1]].load is not None:
                        self.list_node[bi_output[1]].load = self.list_node[bi_input[0]].load + self.list_node[bi_input[0]].ext_load - bi_axial_load
        
        if (bi_axial_load is None) and (be_axial_load is not None):
            if self.list_node[be_input[0]].load is not None:
                self.list_node[be_output[1]].load = self.list_node[be_input[0]].load + self.list_node[be_input[0]].ext_load
            if self.list_node[be_input[1]].load is not None:
                self.list_node[be_output[0]].load = self.list_node[be_input[1]].load + self.list_node[be_input[1]].ext_load
        if (be_axial_load is None) and (bi_axial_load is not None):
            if self.list_node[bi_input[0]].load is not None:
                self.list_node[bi_output[1]].load = self.list_node[bi_input[0]].load + self.list_node[bi_input[0]].ext_load
            if self.list_node[bi_input[1]].load is not None:
                self.list_node[bi_output[0]].load = self.list_node[bi_input[1]].load + self.list_node[bi_input[1]].ext_load

class RadialBearing(LoadBearing):
    def __init__(self, d, D, B, i, Z, Dw, alpha, Cr=None, C0r=None ,oil=oil_iso_vg_1500, 
                 material=material_iso, typ_contact=None, mass=None):
        self.d = d
        self.D = D
        self.B = B
        self.Dpw = (d + D)/2.
        self.i = i
        self.Z = Z
        self.Dw = Dw
        self.alpha = alpha
        self.oil = oil
        self.material = material
        self.typ_contact = typ_contact
        self.mass = mass
        self.direction = 0
        self.sub_direction = 0
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
        self.jeu = (self.E-self.F-2*self.Dw)/4.
    
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

        if total_cycles == 0.:
            return math.inf
        else:
            Pr = (Pr / total_cycles) ** (1/self.coeff_baselife)
            L10 = (Cr/Pr)**(self.coeff_baselife)
            return L10
        
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
#        print(Fr, Fa)
        for fr, fa, n, ti, Ti  in zip(Fr, Fa, N, t, T):
            if (((fr != 0.) or (fa != 0.)) and (n > 0.)):
                cycles = n * ti * 2 * math.pi
                total_cycles += cycles
                
                a1 = ((1-self.material.c_gamma)
                     * (npy.log(1/S)/npy.log(100/90.))**(1/self.material.weibull_e)
                     + self.material.c_gamma)
                L10 = self.BaseLifeTime([fr], [fa], [n], [ti], self.Cr)
                Pr = self.EquivalentDynamicLoad(fr, fa)
                # viscosité cinématique de référence
                if n < (1000*2*npy.pi/60.):
                    nu1 = 45000*(n*60/(2*npy.pi))**(-0.83)*(self.Dpw*1e3)**(-0.5)
                else:
                    nu1 = 4500*(n*60/(2*npy.pi))**(-0.5)*(self.Dpw*1e3)**(-0.5)
                
                coeff_oil = self.oil.oil_kinematic_viscosity_curve
                nu = 10**(10**(coeff_oil['A']*npy.log10(Ti)+coeff_oil['B']))-0.6
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
    
    def CheckFNRRules(self, Fr, Fa, N):
        check_rules = True
        val_rules = npy.inf
        for fr, fa, n  in zip(Fr, Fa, N):
            if self.typ_bearing == 'radial_roller_bearing':
                rules_snr = RadialRollerBearingSNR(self.d, self.D, self.B, self.Z, self.alpha, self.Dpw)
                check_rules_iter, val_rules_iter = rules_snr.RuleAxialLoad(fr, fa, n, level_axial_load='constant_load')
                val_rules = min(val_rules_iter, val_rules)
                if check_rules_iter == False:
                    check_rules = False
        return check_rules, val_rules

    def VolumeModel(self, center = (0,0,0), axis = (1,0,0)):
        center = vm.Point3D(npy.round(center,6))
        x = vm.Vector3D(axis)
        x.vector = x.vector/x.Norm()
        
        y = x.RandomUnitNormalVector()
        y.vector = npy.round(y.vector,3)
        y.vector = y.vector/y.Norm() 
        
        z=vm.Vector3D(npy.cross(x.vector,y.vector))
        
        #Internal Ring
        IRC=self.InternalRingContour()        
        irc=primitives3D.RevolvedProfile(center, x, z, [IRC], center,
                                         x,angle=2*math.pi, name='Internal Ring')
        #External Ring
        ERC=self.ExternalRingContour()
        erc=primitives3D.RevolvedProfile(center,x, z, [ERC], center,
                                         x, angle=2*math.pi,name='External Ring')
        #roller
        ROL=self.RollingContourCAD()
        
        radius=self.F/2.+self.jeu+self.Dw/2.
        rol=[]
        theta=2*npy.pi/self.Z
        
        for zi in range(int(self.Z)):
            center_roller = center + radius*math.cos(zi*theta) * y + radius*math.sin(zi*theta) * z
            rol.append(primitives3D.RevolvedProfile(center_roller, x, z, [ROL],
                                                    center_roller, x,
                                                    angle=2*math.pi,name='Roller {}'.format(zi+1)))
        
        tot=[irc,erc]+rol
        model=vm.VolumeModel([('', tot)])
        return model   
    
    def FreeCADExport(self, fcstd_filepath, python_path='python', 
                      path_lib_freecad='/usr/lib/freecad/lib',export_types=['fcstd']):
        model = self.VolumeModel()
        model.FreeCADExport(fcstd_filepath,python_path,path_lib_freecad,export_types)
        
    def Graph(self, list_node=None, num=None, pos_x=0):
        
        self.num = num
        if num is None:
            num_init = 0
        else:
            num_init = num
            
        if list_node is None:
            li = []
            for n in npy.arange(num_init, num_init + 8, 1):
                li.append(LoadNode(n, 0)) 
            list_node = li
        
        self.next = {}
        for nd in list_node:
            self.next[nd] = []
        self.next[list_node[4]].append(list_node[0])
        self.next[list_node[1]].append(list_node[5])
        self.next[list_node[6]].append(list_node[2])
        self.next[list_node[3]].append(list_node[7])
        self.list_node = list_node
        
    def ResumeLoad(self):
        
        if (self.direction == 1) or (self.sub_direction == 1):
            if self.list_node[3].load is not None:
                if self.list_node[7].load is not None:
                    transversale_load = self.list_node[3].load + self.list_node[3].ext_load - self.list_node[7].load
                    internal_ring_load = self.list_node[7].load
                else:
                    transversale_load = self.list_node[3].load + self.list_node[3].ext_load
                    internal_ring_load = 0
            else:
                transversale_load = 0
                internal_ring_load = 0
            if self.list_node[4].load is not None:
                if self.list_node[0].load is not None:
                    external_ring_load = max(self.list_node[0].load, self.list_node[1].load)
                else:
                    external_ring_load = 0
            else:
                external_ring_load = 0
            
        elif (self.direction == -1) or (self.sub_direction == -1):
            if self.list_node[6].load is not None:
                if self.list_node[2].load is not None:
                    transversale_load = self.list_node[6].load + self.list_node[6].ext_load - self.list_node[2].load
                    internal_ring_load = self.list_node[2].load
                else:
                    transversale_load = self.list_node[6].load + self.list_node[6].ext_load
                    internal_ring_load = 0
            else:
                transversale_load = 0
                internal_ring_load = 0
            if self.list_node[1].load is not None:
                if self.list_node[5].load is not None:
                    external_ring_load = max(self.list_node[4].load, self.list_node[5].load)
                else:
                    external_ring_load = 0
            else:
                external_ring_load = 0
        else:
            transversale_load = 0
            if self.list_node[5].load is not None:
                external_ring_load = self.list_node[5].load
            else:
                external_ring_load = 0
            if self.list_node[2].load is not None:
                external_ring_load = self.list_node[2].load
            else:
                external_ring_load = 0
                
        self.transversale_load = transversale_load
        self.external_ring_load = external_ring_load
        self.internal_ring_load = internal_ring_load

    def PlotLoad(self, a, pos=0):
        
        delta_1 = (self.D - self.D1)/10.
        delta_2 = (self.d1 - self.d)/10.
        positions = {}
        positions[self.list_node[0]] = vm.Point2D((-self.B/2., self.D/2. - delta_1))
        positions[self.list_node[1]] = vm.Point2D((-self.B/2., self.D/2. - 3*delta_1))
        positions[self.list_node[2]] = vm.Point2D((-self.B/2., self.d/2. + 3*delta_2))
        positions[self.list_node[3]] = vm.Point2D((-self.B/2., self.d/2. + delta_2))
        positions[self.list_node[4]] = vm.Point2D((self.B/2., self.D/2. - delta_1))
        positions[self.list_node[5]] = vm.Point2D((self.B/2., self.D/2. - 3*delta_1))
        positions[self.list_node[6]] = vm.Point2D((self.B/2., self.d/2. + 3*delta_2))
        positions[self.list_node[7]] = vm.Point2D((self.B/2., self.d/2. + delta_2))
        
        max_load = 0
        for nd in self.list_node:
            if nd.load is not None:
                max_load = max(nd.load, max_load)

        if (self.direction == 1) or (self.sub_direction == 1):
            if self.transversale_load != 0:
                list_line = [vm.LineSegment2D(positions[self.list_node[3]], positions[self.list_node[5]])]
                list_line.append(vm.LineSegment2D(positions[self.list_node[4]], positions[self.list_node[2]]))
                load_arrow = vm.Contour2D(list_line)
                load_arrow = load_arrow.Translation((pos, 0), True)
                load_arrow.MPLPlot(a,'-b', True, width = self.B/10./max_load*self.transversale_load)
        elif (self.direction == -1) or (self.sub_direction == -1):
            if self.transversale_load != 0:
                list_line = [vm.LineSegment2D(positions[self.list_node[6]], positions[self.list_node[0]])]
                list_line.append(vm.LineSegment2D(positions[self.list_node[1]], positions[self.list_node[7]]))
                load_arrow = vm.Contour2D(list_line)
                load_arrow = load_arrow.Translation((pos, 0), True)
                load_arrow.MPLPlot(a,'-b', True, width = self.B/10./max_load*self.transversale_load)
        if self.internal_ring_load != 0:
            list_line = [vm.LineSegment2D(positions[self.list_node[3]], positions[self.list_node[7]])]
            list_line.append(vm.LineSegment2D(positions[self.list_node[6]], positions[self.list_node[2]]))
            load_arrow = vm.Contour2D(list_line)
            load_arrow = load_arrow.Translation((pos, 0), True)
            load_arrow.MPLPlot(a,'-b', True, width = self.B/10./max_load*self.internal_ring_load)
        if self.external_ring_load != 0:
            list_line = [vm.LineSegment2D(positions[self.list_node[1]], positions[self.list_node[5]])]
            list_line.append(vm.LineSegment2D(positions[self.list_node[4]], positions[self.list_node[0]]))
            load_arrow = vm.Contour2D(list_line)
            load_arrow = load_arrow.Translation((pos, 0), True)
            load_arrow.MPLPlot(a,'-b', True, width = self.B/10./max_load*self.external_ring_load)

    
    def PlotGraph(self):
        
        delta_1 = (self.D - self.D1)/10.
        delta_2 = (self.d1 - self.d)/10.
        positions = {}
        positions[self.list_node[0]] = vm.Point2D((-self.B/2., self.D/2. - delta_1))
        positions[self.list_node[1]] = vm.Point2D((-self.B/2., self.D/2. - 2*delta_1))
        positions[self.list_node[2]] = vm.Point2D((-self.B/2., self.d/2. + 2*delta_2))
        positions[self.list_node[3]] = vm.Point2D((-self.B/2., self.d/2. + delta_2))
        positions[self.list_node[4]] = vm.Point2D((self.B/2., self.D/2. - delta_1))
        positions[self.list_node[5]] = vm.Point2D((self.B/2., self.D/2. - 2*delta_1))
        positions[self.list_node[6]] = vm.Point2D((self.B/2., self.d/2. + 2*delta_2))
        positions[self.list_node[7]] = vm.Point2D((self.B/2., self.d/2. + delta_2))
        
        list_line = []
        check_line = []
        for nd in self.list_node:
            np = self.next[nd]
            for n in np:
                if (n in self.list_node) and [nd, n] not in check_line:
                    n1 = positions[nd]
                    n2 = positions[n]
                    list_line.append(vm.LineSegment2D(n1, n2))
                    check_line.append([nd, n])
                    
        load_arrow = vm.Contour2D(list_line)
        return load_arrow
    
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

        li_nd = []
        for nd in self.list_node:
            li_nd.append(nd.Dict())
        d['list_node'] = li_nd
        
        d['material'] = self.material.Dict()
        d['oil'] = self.oil.Dict()
        
        di_next = {}
        for key, val in self.next.items():
            pos_obj = self.list_node.index(key)
            di_next[pos_obj] = []
            for v in val:
                pos_v = self.list_node.index(v)
                di_next[pos_obj].append(pos_v)
        d['next'] = di_next
        
        di_position = {}
        for key, val in self.positions.items():
            pos_obj = self.list_node.index(key)
            di_position[pos_obj] = val
        d['positions'] = di_position
        
        return d
        
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
class ConceptRadialBallBearing(RadialBearing):
    def __init__(self, d, D, B, i, Z, Dw, alpha, Cr=None, C0r=None ,oil=oil_iso_vg_1500, 
                 material=material_iso, typ_contact=None, mass=None):
        RadialBearing.__init__(self, d, D, B, i, Z, Dw, alpha, Cr, C0r, oil, material, typ_contact, mass)
        self.coeff_baselife = 3.
        self.name = 'RadialBallBearing'
        
        # estimation for the graph 2D description
        h1 = self.Dw/2. - (self.E - self.D1)/2.
        self.h = self.B/2. - self.Dw/2.*npy.sin(npy.arccos(h1/(self.Dw/2.))) - 1e-4
        
        
    def EquivalentStaticLoad(self, fr, fa=None):
        #Charge radiale statique équivalente
        X0 = 0.6
        Y0 = 0.5
        P0r = max(fr,X0*fr+Y0*fa)
        return P0r
    
    def EquivalentDynamicLoad(self, fr, fa = 0):
        alphap = fsolve((lambda alphap:npy.cos(5/180*npy.pi)/npy.cos(alphap) \
                    -(1+0.012534*(fa/(self.i*self.Z*((self.Dw*1e3)**2) \
                    *npy.sin(alphap)))**(2/3.))),self.alpha+1)[0]
#        alphap = 0.001
        ksi = 1.05
        nu = 1-npy.sin(5/180.*npy.pi)/2.5
        e = ksi*npy.tan(alphap)
        X1 = 1-0.4*ksi/nu
        X3 = 1
        Y1 = 0.4/nu*1/npy.tan(alphap)
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
        pbi1 = pbi2.Translation((self.h, 0))
        pbi3 = vm.Point2D((-self.B/2., self.d/2.))
        pbi4 = vm.Point2D((self.B/2., self.d/2.))
        pbi5 = vm.Point2D((self.B/2., self.d1/2.))
        pbi6 = pbi5.Translation((-self.h, 0))
        bi1 = primitives2D.RoundedLineSegments2D([pbi1, pbi2, pbi3, pbi4, pbi5, pbi6], {1: self.radius, 
                                                 2: self.radius, 3: self.radius, 4: self.radius}, False)
        cbi1 = vm.Arc2D(pbi1, vm.Point2D((0, self.F/2)), pbi6)
        return vm.Contour2D([bi1, cbi1])
        
    def ExternalRingContour(self):
        
        pbe2 = vm.Point2D((-self.B/2., self.D1/2.))
        pbe1 = pbe2.Translation((self.h, 0))
        pbe3 = vm.Point2D((-self.B/2., self.D/2.))
        pbe4 = vm.Point2D((self.B/2., self.D/2.))
        pbe5 = vm.Point2D((self.B/2., self.D1/2.))
        pbe6 = pbe5.Translation((-self.h, 0))
        be1 = primitives2D.RoundedLineSegments2D([pbe1, pbe2, pbe3, pbe4, pbe5, pbe6], {1: self.radius, 
                                                 2: self.radius, 3: self.radius, 4: self.radius}, False)
        cbe1 = vm.Arc2D(pbe1, vm.Point2D((0, self.E/2)), pbe6)
        return vm.Contour2D([be1, cbe1])
    
    def RollingContour(self):
        
        p0 = vm.Point2D((0, 0))
        c1 = vm.Circle2D(p0, self.Dw/2.) 
        return vm.Contour2D([c1])
    
    def Plot(self):
        
        be_sup = self.ExternalRingContour()
        bi_sup = self.InternalRingContour()
        ball_sup = self.RollingContour()
        ball_sup.Translation((0, self.Dpw/2.))
        
        bearing_sup = vm.Contour2D([be_sup, bi_sup, ball_sup])
        bearing_inf = bearing_sup.Rotation(vm.Point2D((0, 0)), npy.pi, True)
        
        bg = vm.Contour2D([bearing_sup, bearing_inf])
        return bg
    
    def Graph(self, list_node=None, num=None, pos_x=0, sub_direction=1):
                
        RadialBearing.Graph(self, list_node, num, pos_x)
        decal_x = 0
        decal_y = 0.1
    
        li = self.list_node
        self.sub_direction = sub_direction
        if self.sub_direction == 1:
            self.next[list_node[4]].append(list_node[2])
            self.next[list_node[3]].append(list_node[5])
        elif self.sub_direction == -1:
            self.next[list_node[1]].append(list_node[7])
            self.next[list_node[6]].append(list_node[0])
        else:
            self.next[list_node[4]].append(list_node[2])
            self.next[list_node[3]].append(list_node[5])
            self.next[list_node[1]].append(list_node[7])
            self.next[list_node[6]].append(list_node[0])
            
        li_pos = {li[0]: (pos_x, 1), li[1]: (pos_x, 1 - decal_y), 
          li[2]: (pos_x + decal_x, decal_y), li[3]: (pos_x + decal_x, 0),
          li[4]: (pos_x + 1 - decal_x, 1), li[5]: (pos_x + 1 - decal_x, 1 - decal_y), 
          li[6]: (pos_x + 1, decal_y), li[7]: (pos_x + 1, 0)}
        self.positions = li_pos
        
    @classmethod
    def Dict2Obj(cls, d):
        if 'Cr' not in d.keys():
            d['Cr'] = None
        if 'C0r' not in d.keys():
            d['C0r'] = None
        obj = cls(d = d['d'], D = d['D'], B = d['B'], i = d['i'], Z = d['Z'], 
                  Dw = d['Dw'], alpha = d['alpha'], Cr = d['Cr'], C0r = d['C0r'],
                  oil = Oil.Dict2Obj(d['oil']), 
                  material = Material.Dict2Obj(d['material']),  
                  typ_contact = d['typ_contact'], mass = d['mass'])
        return obj
        
class ConceptAngularBallBearing(RadialBearing):
    def __init__(self, d, D, B, i, Z, Dw, alpha, Cr=None, C0r=None ,oil=oil_iso_vg_1500, 
                 material=material_iso, typ_contact=None, direction=1, mass=None):
        RadialBearing.__init__(self, d, D, B, i, Z, Dw, alpha, Cr, C0r, oil, material, 
                               typ_contact, mass)
        self.coeff_baselife = 3.
        self.direction = direction
        self.name = 'AngularBallBearing'
        
        # estimation for the graph 2D description
        h1 = self.Dw/2. - (self.E - self.D1)/2.
        self.h1 = self.B/2. - self.Dw/2.*npy.sin(npy.arccos(h1/(self.Dw/2.))) - 1e-4
        h2 = 0.95*self.Dw/2.
        self.h2 = self.B/2. - self.Dw/2.*npy.sin(npy.arccos(h2/(self.Dw/2.))) - 1e-4
        self.D2 = 0.6*(self.D - self.E) + self.E
        self.d2 = 0.6*(self.F - self.d) + self.d
        
    def EquivalentStaticLoad(self, fr, fa=None):
        #Charge radiale statique équivalente
        if self.i == 1:
            X0 = 0.5
            evol_Y0 = [0.52, 0.5, 0.46, 0.42, 0.38, 0.33, 0.29, 0.26, 0.22]
            evol_alpha = npy.array([5, 10, 15, 20, 25, 30, 35, 40, 45])/180*npy.pi
            f = interpolate.interp1d(list(evol_alpha),evol_Y0, fill_value='extrapolate')
            Y0 = float(f(self.alpha))
        elif self.i == 2:
            X0 = 1
            evol_Y0 = [1.04, 1, 0.92, 0.84, 0.76, 0.66, 0.58, 0.52, 0.44]
            evol_alpha = npy.array([5, 10, 15, 20, 25, 30, 35, 40, 45])/180*npy.pi
            f = interpolate.interp1d(list(evol_alpha),evol_Y0, fill_value='extrapolate')
            Y0 = float(f(self.alpha))
        P0r = max(fr,X0*fr+Y0*fa)
        return P0r
    
    def EquivalentDynamicLoad(self, fr, fa = 0):
        alphap = fsolve((lambda alphap:npy.cos(self.alpha)/npy.cos(alphap) \
                    -(1+0.012534*(fa/(self.i*self.Z*((self.Dw*1e3)**2)*npy.sin(alphap)))**(2/3.))),self.alpha)[0]
        if self.alpha <= 5/180.*npy.pi:
            if self.i == 1:
                ksi = 1.05
            elif self.i == 2:
                ksi = 1.25
        else:
            ksi = 1.25
        if self.alpha <= 5/180.*npy.pi:
            nu = 1-npy.sin(5/180.*npy.pi)/2.5
        elif self.alpha <= 15/180.*npy.pi:
            nu = 1-npy.sin(self.alpha)/2.5
        else:
            nu = 1-npy.sin(self.alpha)/2.75
        if self.alpha <= 15/180.*npy.pi:
            e = ksi*npy.tan(alphap)
        X1 = 1-0.4*ksi/nu
        X3 = 1
        if self.alpha <= 15/180.*npy.pi:
            Y1 = 0.4/nu*1/npy.tan(alphap)
        else:
            alpha1 = npy.arccos(npy.cos(self.alpha)*0.9724)
            Y1 = fsolve((lambda Y1:Y1-(0.4/npy.tan(alpha1))/(1-(1/3.)*npy.sin(alpha1))),1)[0]
        if self.alpha > 15/180.*npy.pi:
            e = (1-X1)/Y1
            Y3 = 0.625/e
        elif self.alpha <= 15/180.*npy.pi:
            Y3 = 0.625/ksi*(1/npy.tan(alphap))
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
    
    def InternalRingContour(self, sign_V=1):
        
        def graph(sign_V, sign_H):
            pbi2 = vm.Point2D((-sign_H*self.B/2., sign_V*self.d2/2.))
            pbi1 = vm.Point2D((-sign_H*(self.B/2. - self.h2), sign_V*(self.Dpw/2. - self.Dw/2.*0.95)))
            pbi3 = vm.Point2D((-sign_H*self.B/2., sign_V*self.d/2.))
            pbi4 = vm.Point2D((sign_H*self.B/2., sign_V*self.d/2.))
            pbi5 = vm.Point2D((sign_H*self.B/2., sign_V*self.d1/2.))
            pbi6 = pbi5.Translation((-sign_H*self.h1, 0))
            bi1 = primitives2D.RoundedLineSegments2D([pbi1, pbi2, pbi3, pbi4, pbi5, pbi6], {1: self.radius, 
                                                 2: self.radius, 3: self.radius, 4: self.radius}, False, adapt_radius = True)
            cbi1 = vm.Arc2D(pbi1, vm.Point2D((0, sign_V*self.F/2)), pbi6)
            bearing = vm.Contour2D([bi1, cbi1])
            return bearing
        sup_contour = graph(sign_V, sign_H = -self.direction)
        bearing1 = vm.Contour2D([sup_contour])
        return bearing1
        
    def ExternalRingContour(self, sign_V=1):
        
        def graph(sign_V, sign_H):
            pbe2 = vm.Point2D((-sign_H*self.B/2., sign_V*self.D1/2.))
            pbe1 = pbe2.Translation((sign_H*self.h1, 0))
            pbe3 = vm.Point2D((-sign_H*self.B/2., sign_V*self.D/2.))
            pbe4 = vm.Point2D((sign_H*self.B/2., sign_V*self.D/2.))
            pbe5 = vm.Point2D((sign_H*self.B/2., sign_V*self.D2/2.))
            pbe6 = vm.Point2D((sign_H*(self.B/2. - self.h2), sign_V*(self.Dpw/2. + self.Dw/2.*0.95)))
            be1 = primitives2D.RoundedLineSegments2D([pbe1, pbe2, pbe3, pbe4, pbe5, pbe6], {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius}, False, adapt_radius = True)
            cbe1 = vm.Arc2D(pbe1, vm.Point2D((0, sign_V*self.E/2)), pbe6)
            bearing = vm.Contour2D([be1, cbe1])
            return bearing
        sup_contour = graph(sign_V, sign_H = -self.direction)
        bearing1 = vm.Contour2D([sup_contour])
        return bearing1
    
    def RollingContour(self):
        
        p0 = vm.Point2D((0, 0))
        c1 = vm.Circle2D(p0, self.Dw/2.) 
        return vm.Contour2D([c1])
    
    def Plot(self):
        
        be_sup = self.ExternalRingContour(1)
        be_inf = self.ExternalRingContour(-1)
        bi_sup = self.InternalRingContour(1)
        bi_inf = self.InternalRingContour(-1)
        ball = self.RollingContour()
        ball_sup = ball.Translation((0, self.Dpw/2.), True)
        ball_inf = ball.Translation((0, -self.Dpw/2.), True)
        bg = vm.Contour2D([be_sup, bi_sup, ball_sup, be_inf, bi_inf, ball_inf])
        return bg
    
    def Graph(self, list_node=None, num=None, pos_x=0):
                
        RadialBearing.Graph(self, list_node, num, pos_x)
        decal_x = 0.3
        decal_y = 0.1
    
        li = self.list_node
        if self.direction == 1:
            self.next[list_node[4]].append(list_node[2])
            self.next[list_node[3]].append(list_node[5])
            li_pos = {li[0]: (pos_x + decal_x, 1), li[1]: (pos_x + decal_x, 1 - decal_y), 
              li[2]: (pos_x, decal_y), li[3]: (pos_x, 0),
              li[4]: (pos_x + 1, 1), li[5]: (pos_x + 1, 1 - decal_y), 
              li[6]: (pos_x + 1 - decal_x, decal_y), li[7]: (pos_x + 1 - decal_x, 0)}
        elif self.direction == -1:
            self.next[list_node[1]].append(list_node[7])
            self.next[list_node[6]].append(list_node[0])
            li_pos = {li[0]: (pos_x, 1), li[1]: (pos_x, 1 - decal_y), 
              li[2]: (pos_x + decal_x, decal_y), li[3]: (pos_x + decal_x, 0),
              li[4]: (pos_x + 1 - decal_x, 1), li[5]: (pos_x + 1 - decal_x, 1 - decal_y), 
              li[6]: (pos_x + 1, decal_y), li[7]: (pos_x + 1, 0)}
        self.positions = li_pos
        
    @classmethod
    def Dict2Obj(cls, d):
        if 'Cr' not in d.keys():
            d['Cr'] = None
        if 'C0r' not in d.keys():
            d['C0r'] = None
        obj = cls(d = d['d'], D = d['D'], B = d['B'], i = d['i'], Z = d['Z'], 
                  Dw = d['Dw'], alpha = d['alpha'], Cr = d['Cr'], C0r = d['C0r'],
                  oil = Oil.Dict2Obj(d['oil']), 
                  material = Material.Dict2Obj(d['material']),  
                  typ_contact = d['typ_contact'], direction = d['direction'],
                  mass = d['mass'])
        return obj
            
class ConceptSphericalBallBearing(RadialBearing):
    def __init__(self, d, D, B, i, Z, Dw, alpha, Cr=None, C0r=None ,oil=oil_iso_vg_1500, 
                 material=material_iso, typ_contact=None, mass=None):
        RadialBearing.__init__(self, d, D, B, i, Z, Dw, alpha, Cr, C0r, oil, material, typ_contact, mass)
        self.coeff_baselife = 3.
        self.name = 'SphericalBallBearing'
        
    def EquivalentStaticLoad(self, fr, fa=None):
        #Charge radiale statique équivalente
        if self.i == 1:
            X0 = 0.5
            Y0 = 0.22/npy.tan(self.alpha)
        elif self.i == 2:
            X0 = 1
            Y0 = 0.44/npy.tan(self.alpha)
        P0r = max(fr,X0*fr+Y0*fa)
        return P0r
    
    def EquivalentDynamicLoad(self, fr, fa = 0):
        alphap = self.alpha
        ksi = 1.5
        nu = 1
        e = ksi*npy.tan(alphap)
        X1 = 1-0.4*ksi/nu
        X3 = 1
        Y1 = 0.4/nu*(1/npy.tan(alphap))
        Y3 = 0.625/ksi*(1/npy.tan(alphap))
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
    def Dict2Obj(cls, d):
        if 'Cr' not in d.keys():
            d['Cr'] = None
        if 'C0r' not in d.keys():
            d['C0r'] = None
        obj = cls(d = d['d'], D = d['D'], B = d['B'], i = d['i'], Z = d['Z'], 
                  Dw = d['Dw'], alpha = d['alpha'], Cr = d['Cr'], C0r = d['C0r'],
                  oil = Oil.Dict2Obj(d['oil']), 
                  material = Material.Dict2Obj(d['material']),  
                  typ_contact = d['typ_contact'], 
                  mass = d['mass'])
        return obj
            
class ConceptRadialRollerBearing(RadialBearing):
    def __init__(self, d, D, B, i, Z, Dw, alpha=0, Cr=None, C0r=None ,oil=oil_iso_vg_1500, 
                 material=material_iso, typ_contact='linear_contact', direction=1, typ='N', mass=None):
        RadialBearing.__init__(self, d, D, B, i, Z, Dw, alpha, Cr, C0r, oil, material, typ_contact, mass)
        self.coeff_baselife = 10/3.
        self.direction = direction
        self.typ = typ
        self.name = 'RadialRollerBearing'
        
        # estimation for the graph 2D description
        self.Dpw = (self.d + self.D)/2.
        self.Lw = 0.7*self.B
        self.h = self.B/2. - self.Lw/2. - 1e-4
        
    def EquivalentStaticLoad(self, fr, fa=None):
        #Charge radiale statique équivalente
        if self.alpha != 0:
            x0 = 0.5*self.i
            y0 = 0.22*1/npy.tan(self.alpha)*self.i
        else:
            x0,y0 = 1,0
        P0r = max(fr,X0*fr+Y0*fa)
        return P0r
    
    def EquivalentDynamicLoad(self, fr, fa = 0):
        
        ksi = 1.5 #param of the ISO 1281
        e = ksi*npy.tan(self.alpha)
        nu = 1-0.15*npy.sin(self.alpha)
        w = self.material.weibull_e*self.coeff_baselife
  
        if self.typ_contact == 'point_contact':
            if self.i == 1:
                Jr0p5 = 0.2288
                Ja0p5 = 0.2782
                J10p5 = 0.5625
                J20p5 = 0.5875
            else: #analyse if the parameters are true for i>2
                Jr0p5 = 0.4577
                Ja0p5 = 0.
                J10p5 = 0.6925
                J20p5 = 0.7233
        elif self.typ_contact == 'linear_contact':
            if self.i == 1:
                Jr0p5 = 0.2453
                Ja0p5 = 0.3090
                J10p5 = 0.6495
                J20p5 = 0.6744
            else: #analyse if the parameters are true for i>2
                Jr0p5 = 0.4906
                Ja0p5 = 0
                J10p5 = 0.7577
                J20p5 = 0.7867
        elif self.typ_contact == 'mixed_contact':
            if self.i == 1:
                Jr0p5 = 0.2369
                Ja0p5 = 0.2932
                J10p5 = 0.6044
                J20p5 = 0.6295
            else: #analyse if the parameters are true for i>2
                Jr0p5 = 0.4739
                Ja0p5 = 0
                J10p5 = 0.7244
                J20p5 = 0.7543
                
        X1 = 1-Jr0p5/(J10p5*J20p5)**0.5*ksi/nu        
        X2 = 2**(1-(1/w))*X1
        X3 = 1
        if self.alpha > 0:
            Y1 = Jr0p5*(1/npy.tan(self.alpha))/(J10p5*J20p5)**0.5*1/nu
            Y2 = 2**(1-(1/w))*Y1
            Y3 = 1/ksi*(2**(1-(1/w))-1)*(1/npy.tan(self.alpha))
            
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
    
    def InternalRingContour(self, sign_V=1):
        
        def graph(sign_V, sign_H):
            if self.typ in ['NUP', 'N', 'NF']:
                d1 = self.d1
                pbi2 = vm.Point2D((-sign_H*self.B/2., sign_V*d1/2.))
                pbi1 = vm.Point2D((-sign_H*(self.B/2. - self.h), sign_V*d1/2.))
                pbi0 = vm.Point2D((-sign_H*(self.B/2. - self.h), sign_V*(self.F/2.)))
                pbi3 = vm.Point2D((-sign_H*self.B/2., sign_V*self.d/2.))
                pbi4 = vm.Point2D((sign_H*self.B/2., sign_V*self.d/2.))
                pbi5 = vm.Point2D((sign_H*self.B/2., sign_V*d1/2.))
                pbi6 = vm.Point2D((sign_H*(self.B/2. - self.h), sign_V*d1/2.))
                pbi7 = vm.Point2D((sign_H*(self.B/2. - self.h), sign_V*(self.F/2.)))
                bi1 = primitives2D.RoundedLineSegments2D([pbi0, pbi1, pbi2, pbi3, pbi4, pbi5, pbi6, pbi7],
                                   {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius, 
                                    5: self.radius, 6: self.radius}, True, adapt_radius = True)
            elif self.typ == 'NJ':
                d1 = self.d1
                d2 = self.F - 0.1*(self.F - self.d)
                pbi2 = vm.Point2D((-sign_H*self.B/2., sign_V*d1/2.))
                pbi1 = vm.Point2D((-sign_H*(self.B/2. - self.h), sign_V*d1/2.))
                pbi0 = vm.Point2D((-sign_H*(self.B/2. - self.h), sign_V*(self.F/2.)))
                pbi3 = vm.Point2D((-sign_H*self.B/2., sign_V*self.d/2.))
                pbi4 = vm.Point2D((sign_H*self.B/2., sign_V*self.d/2.))
                pbi5 = vm.Point2D((sign_H*self.B/2., sign_V*d2/2.))
                pbi6 = vm.Point2D((sign_H*(self.B/2. - self.h), sign_V*(self.F/2.)))
                bi1 = primitives2D.RoundedLineSegments2D([pbi0, pbi1, pbi2, pbi3, pbi4, pbi5, pbi6],
                                   {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius, 
                                    5: self.radius}, True, adapt_radius = True)
            elif self.typ == 'NU':
                d1 = self.F - 0.1*(self.F - self.d)
                pbi2 = vm.Point2D((-sign_H*self.B/2., sign_V*d1/2.))
                pbi1 = vm.Point2D((-sign_H*(self.B/2. - self.h), sign_V*(self.F/2.)))
                pbi3 = vm.Point2D((-sign_H*self.B/2., sign_V*self.d/2.))
                pbi4 = vm.Point2D((sign_H*self.B/2., sign_V*self.d/2.))
                pbi5 = vm.Point2D((sign_H*self.B/2., sign_V*d1/2.))
                pbi6 = vm.Point2D((sign_H*(self.B/2. - self.h), sign_V*(self.F/2.)))
                bi1 = primitives2D.RoundedLineSegments2D([pbi1, pbi2, pbi3, pbi4, pbi5, pbi6], {1: self.radius, 
                                   2: self.radius, 3: self.radius, 4: self.radius}, 
                                   True, adapt_radius = True)
            
            bearing = vm.Contour2D([bi1])
            return bearing
        contour = graph(sign_V, sign_H = self.direction)
        bg = vm.Contour2D([contour])
        return bg
        
    def ExternalRingContour(self, sign_V=1):
        
        def graph(sign_V, sign_H):
            if self.typ in ['NU', 'NJ', 'NUP']:
                D1 = self.D1
                pbe2 = vm.Point2D((-sign_H*self.B/2., sign_V*D1/2.))
                pbe1 = vm.Point2D((-sign_H*(self.B/2. - self.h), sign_V*D1/2.))
                pbe0 = vm.Point2D((-sign_H*(self.B/2. - self.h), sign_V*(self.E/2.)))
                pbe3 = vm.Point2D((-sign_H*self.B/2., sign_V*self.D/2.))
                pbe4 = vm.Point2D((sign_H*self.B/2., sign_V*self.D/2.))
                pbe5 = vm.Point2D((sign_H*self.B/2., sign_V*D1/2.))
                pbe6 = vm.Point2D((sign_H*(self.B/2. - self.h), sign_V*D1/2.))
                pbe7 = vm.Point2D((sign_H*(self.B/2. - self.h), sign_V*(self.E/2.)))
                be1 = primitives2D.RoundedLineSegments2D([pbe0, pbe1, pbe2, pbe3, pbe4, pbe5, pbe6, pbe7],
                                   {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius, 
                                    5: self.radius, 6: self.radius}, True)
            elif self.typ == 'NF':
                D1 = self.D1
                D2 = self.E + 0.1*(self.D - self.E)
                pbe2 = vm.Point2D((-sign_H*self.B/2., sign_V*D1/2.))
                pbe1 = vm.Point2D((-sign_H*(self.B/2. - self.h), sign_V*D1/2.))
                pbe0 = vm.Point2D((-sign_H*(self.B/2. - self.h), sign_V*(self.E/2.)))
                pbe3 = vm.Point2D((-sign_H*self.B/2., sign_V*self.D/2.))
                pbe4 = vm.Point2D((sign_H*self.B/2., sign_V*self.D/2.))
                pbe5 = vm.Point2D((sign_H*self.B/2., sign_V*D2/2.))
                pbe6 = vm.Point2D((sign_H*(self.B/2. - self.h), sign_V*(self.E/2.)))
                be1 = primitives2D.RoundedLineSegments2D([pbe0, pbe1, pbe2, pbe3, pbe4, pbe5, pbe6],
                                   {1: self.radius, 2: self.radius, 3: self.radius, 4: self.radius, 
                                    5: self.radius}, True)
            elif self.typ == 'N':
                D1 = self.E + 0.1*(self.D - self.E)
                pbe2 = vm.Point2D((-sign_H*self.B/2., sign_V*D1/2.))
                pbe1 = vm.Point2D((-sign_H*(self.B/2. - self.h), sign_V*(self.E/2.)))
                pbe3 = vm.Point2D((-sign_H*self.B/2., sign_V*self.D/2.))
                pbe4 = vm.Point2D((sign_H*self.B/2., sign_V*self.D/2.))
                pbe5 = vm.Point2D((sign_H*self.B/2., sign_V*D1/2.))
                pbe6 = vm.Point2D((sign_H*(self.B/2. - self.h), sign_V*(self.E/2.)))
                be1 = primitives2D.RoundedLineSegments2D([pbe1, pbe2, pbe3, pbe4, pbe5, pbe6], {1: self.radius, 
                                   2: self.radius, 3: self.radius, 4: self.radius}, True)
            bearing = vm.Contour2D([be1])
            return bearing
        contour = graph(sign_V, sign_H = self.direction)
        bg = vm.Contour2D([contour])
        return bg
    
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
    
    def Plot(self):
        
        be_sup = self.ExternalRingContour(1)
        be_inf = self.ExternalRingContour(-1)
        bi_sup = self.InternalRingContour(1)
        bi_inf = self.InternalRingContour(-1)
        roller = self.RollingContour()
        roller_sup = roller.Translation((0, self.Dpw/2.), True)
        roller_inf = roller.Translation((0, -self.Dpw/2.), True)
        
        bg = vm.Contour2D([be_sup, bi_sup, roller_sup, be_inf, bi_inf, roller_inf])
        return bg
    
    def Graph(self, list_node=None, num=None, pos_x=0, sub_direction=1):
                
        RadialBearing.Graph(self, list_node, num, pos_x)
        decal_x = 0
        decal_y = 0.1
        self.sub_direction = sub_direction
        
        li = self.list_node
        if self.typ in ['NUP']:
            if self.sub_direction == 1:
                self.next[list_node[4]].append(list_node[2])
                self.next[list_node[3]].append(list_node[5])
            elif self.sub_direction == -1:
                self.next[list_node[1]].append(list_node[7])
                self.next[list_node[6]].append(list_node[0])
            else:
                self.next[list_node[4]].append(list_node[2])
                self.next[list_node[3]].append(list_node[5])
                self.next[list_node[1]].append(list_node[7])
                self.next[list_node[6]].append(list_node[0])
        elif (self.typ == 'NJ' and self.direction == 1) or (self.typ == 'NF' and self.direction == -1):
            self.next[list_node[4]].append(list_node[2])
            self.next[list_node[3]].append(list_node[5])
        elif (self.typ == 'NJ' and self.direction == -1) or (self.typ == 'NF' and self.direction == 1):
            self.next[list_node[1]].append(list_node[7])
            self.next[list_node[6]].append(list_node[0])
            
        li_pos = {li[0]: (pos_x, 1), li[1]: (pos_x, 1 - decal_y), 
          li[2]: (pos_x + decal_x, decal_y), li[3]: (pos_x + decal_x, 0),
          li[4]: (pos_x + 1 - decal_x, 1), li[5]: (pos_x + 1 - decal_x, 1 - decal_y), 
          li[6]: (pos_x + 1, decal_y), li[7]: (pos_x + 1, 0)}
        self.positions = li_pos
    
    @classmethod
    def Dict2Obj(cls, d):
        if 'Cr' not in d.keys():
            d['Cr'] = None
        if 'C0r' not in d.keys():
            d['C0r'] = None
        obj = cls(d = d['d'], D = d['D'], B = d['B'], i = d['i'], Z = d['Z'], 
                  Dw = d['Dw'], alpha = d['alpha'], Cr = d['Cr'], C0r = d['C0r'],
                  oil = Oil.Dict2Obj(d['oil']), 
                  material = Material.Dict2Obj(d['material']),  
                  typ_contact = d['typ_contact'], direction = d['direction'],
                  typ = d['typ'], mass = d['mass'])
        return obj
    
class ConceptTaperedRollerBearing(ConceptRadialRollerBearing, ConceptAngularBallBearing):
    def __init__(self, d, D, B, i, Z, Dw, alpha=0, Cr=None, C0r=None ,oil=oil_iso_vg_1500, 
                 material=material_iso, typ_contact='linear_contact', direction=1, mass=None):
        ConceptRadialRollerBearing.__init__(self, d, D, B, i, Z, Dw, alpha, Cr, C0r, oil, 
                                            material, typ_contact, direction, mass=mass)
        self.coeff_baselife = 10/3.
        self.name = 'TaperedRollerBearing'
        
        # estimation for the graph 2D description
        self.Dpw = (self.d + self.D)/2.
        self.Lw = 0.7*self.B
        self.beta = npy.arctan(self.Dw/self.Dpw*npy.sin(self.alpha))
    
    def InternalRingContour(self, sign_V=1):
        
        shift_bi = 5e-4
        shift_be = 1e-3
        def graph(sign_V, sign_H):
            p0 = vm.Point2D((0, sign_V*self.Dpw/2.))
            p1 = p0.Translation((npy.cos(self.alpha), sign_H*sign_V*npy.sin(self.alpha)), True)
            l1 = vm.Line2D(p0, p1)
            l1.Rotation(p0, -sign_H*sign_V*self.beta)
            l1.Translation((0.8*sign_H*self.Dw/2.*npy.sin(self.alpha), -sign_V*0.8*self.Dw/2.*npy.cos(self.alpha)))
            l2 = l1.Translation((0.2*sign_H*self.Dw/2.*npy.sin(self.alpha), -sign_V*0.2*self.Dw/2.*npy.cos(self.alpha)), True)
            pbi3 = vm.Point2D((-sign_H*(self.B/2. - shift_bi), sign_V*self.d/2.))
            pbi3T = pbi3.Translation((0, 1))
            pbi4 = vm.Point2D((sign_H*(self.B/2.), sign_V*self.d/2.))
            pbi4T = pbi4.Translation((0, 1))
            l3 = vm.Line2D(pbi3, pbi3T)
            l4 = vm.Line2D(pbi4, pbi4T)
            pbi2 = vm.Point2D.LinesIntersection(l1, l3)
            pbi5 = vm.Point2D.LinesIntersection(l1, l4)
            l5 = vm.Line2D(vm.Point2D((-sign_H*self.Lw/2.,0)), vm.Point2D((-sign_H*self.Lw/2.,1)))
            l5.Rotation(vm.Point2D((0,sign_V*self.Dpw/2.)),sign_V*sign_H*self.alpha)
            l6 = vm.Line2D(vm.Point2D((sign_H*self.Lw/2.,0)), vm.Point2D((sign_H*self.Lw/2.,1)))
            l6.Rotation(vm.Point2D((0,sign_V*self.Dpw/2.)),sign_V*sign_H*self.alpha)
            pbi1 = vm.Point2D.LinesIntersection(l1, l5)
            pbi0 = vm.Point2D.LinesIntersection(l2, l5)
            pbi6 = vm.Point2D.LinesIntersection(l1, l6)
            pbi7 = vm.Point2D.LinesIntersection(l2, l6)
            
            bi1 = primitives2D.RoundedLineSegments2D([pbi0, pbi1, pbi2, pbi3, pbi4, pbi5, pbi6, pbi7], 
                                                     {1: self.radius, 2: self.radius, 3: self.radius,
                                                      4: self.radius, 5: self.radius,
                                                      6: self.radius}, True)
            
            bearing = vm.Contour2D([bi1])
            return bearing
        contour = graph(sign_V, sign_H = -self.direction)
        bg = vm.Contour2D([contour])
        return bg
        
    def ExternalRingContour(self, sign_V=1):
        
        shift_bi = 5e-4
        shift_be = 1e-3
        def graph(sign_V, sign_H):
            p0 = vm.Point2D((0, sign_V*self.Dpw/2.))
            p1 = p0.Translation((npy.cos(self.alpha), sign_H*sign_V*npy.sin(self.alpha)), True)
            l0 = vm.Line2D(p0, p1)
            l0.Rotation(p0, sign_H*sign_V*self.beta)
            l0.Translation((-sign_H*self.Dw/2.*npy.sin(self.alpha), sign_V*self.Dw/2.*npy.cos(self.alpha)))
            pbe3 = vm.Point2D((-sign_H*self.B/2., sign_V*self.D/2.))
            pbe3T = pbe3.Translation((0, 1))
            pbe4 = vm.Point2D((sign_H*(self.B/2. - shift_be), sign_V*self.D/2.))
            pbe4T = pbe4.Translation((0, 1))
            l3 = vm.Line2D(pbe3, pbe3T)
            l4 = vm.Line2D(pbe4, pbe4T)
            pbe2 = vm.Point2D.LinesIntersection(l0, l3)
            pbe5 = vm.Point2D.LinesIntersection(l0, l4)
            be1 = primitives2D.RoundedLineSegments2D([pbe2, pbe3, pbe4, pbe5], {0: self.radius, 
                                                 1: self.radius, 2: self.radius, 3: self.radius}, True)
            bearing = vm.Contour2D([be1])
            return bearing
        contour = graph(sign_V, sign_H = -self.direction)
        bg = vm.Contour2D([contour])
        return bg
    
    def RollingContour(self, sign_V=1):
        
        def graph(sign_V, sign_H):
            r1 = vm.Point2D((-sign_H*self.Lw/2., self.Dw/2. - self.Lw/2.*npy.tan(self.beta)))
            r2 = vm.Point2D((sign_H*self.Lw/2., self.Dw/2. + self.Lw/2.*npy.tan(self.beta)))
            r3 = vm.Point2D((sign_H*self.Lw/2., -self.Dw/2. - self.Lw/2.*npy.tan(self.beta)))
            r4 = vm.Point2D((-sign_H*self.Lw/2., -self.Dw/2. + self.Lw/2.*npy.tan(self.beta)))
            rol = primitives2D.RoundedLineSegments2D([r1, r2, r3, r4], {0: self.radius, 
                                                 1: self.radius, 2: self.radius, 3: self.radius}, True)
            return rol
        contour = graph(sign_V, sign_H = -self.direction)
        bg = vm.Contour2D([contour])
        return bg
    
    def Plot(self):
        
        be_sup = self.ExternalRingContour(1)
        be_inf = self.ExternalRingContour(-1)
        bi_sup = self.InternalRingContour(1)
        bi_inf = self.InternalRingContour(-1)
        roller_sup = self.RollingContour(1)
        roller_sup = roller_sup.Rotation(vm.Point2D((0, 0)), self.direction*self.alpha, True)
        roller_sup = roller_sup.Translation((0, self.Dpw/2.), True)
        roller_inf = self.RollingContour(-1)
        roller_inf = roller_inf.Rotation(vm.Point2D((0, 0)), -self.direction*self.alpha, True)
        roller_inf = roller_inf.Translation((0, -self.Dpw/2.), True)
        
        bg = vm.Contour2D([be_sup, bi_sup, roller_sup, be_inf, bi_inf, roller_inf])
        return bg
    
    def Graph(self, list_node=None, num=None, pos_x=0):
        
        ConceptAngularBallBearing.Graph(self, list_node, num, pos_x)
        
    @classmethod
    def Dict2Obj(cls, d):
        if 'Cr' not in d.keys():
            d['Cr'] = None
        if 'C0r' not in d.keys():
            d['C0r'] = None
        obj = cls(d = d['d'], D = d['D'], B = d['B'], i = d['i'], Z = d['Z'], 
                  Dw = d['Dw'], alpha = d['alpha'], Cr = d['Cr'], C0r = d['C0r'],
                  oil = Oil.Dict2Obj(d['oil']), 
                  material = Material.Dict2Obj(d['material']),  
                  typ_contact = d['typ_contact'], direction = d['direction'],
                  mass = d['mass'])
        return obj
   
class ObjBearingAssembly:
    def __init__(self, list_bearing, num):
        self.list_bearing = list_bearing
        self.num = num
    def Load(self):
        for bg in self.list_bearing + self.list_bearing[::-1]:
            bg.CheckNone()
            bg.PropagateNone()
            bg.CheckNone()
        valid = True
        list_valid = [0]
        while valid:
            list_valid_m =list_valid
            list_valid = []
            for bg in self.list_bearing:
                bg.Load()
            for bg in self.list_bearing:
                for nd in bg.list_node:
                    list_valid.append(nd.load)
            if list_valid == list_valid_m:
                valid = False
 
class BearingAssembly:
    def __init__(self, list_bearing, radial_load_linkage, internal_pre_load=0, 
                 connection_bi=['n', 'p'], connection_be=['n', 'p'], behavior_link='pn'):
        self.list_bearing = list_bearing
        self.radial_load_linkage = radial_load_linkage
        self.connection_be = connection_be
        self.connection_bi = connection_bi
        self.behavior_link = behavior_link
        self.mass = 0
        for bg in list_bearing:
            if bg.mass is not None:
                self.mass += bg.mass
        self.B = 0
        for bg in list_bearing:
            self.B += bg.B
        self.D = 0
        for bg in list_bearing:
            self.D = max(self.D, bg.D)
        self.d = npy.inf
        for bg in list_bearing:
            self.d = min(self.d, bg.d)
        
        self.graph, self.li_case = self.Graph(list_bearing)
        if self.graph == []:
            self.check = False
        else:
            self.check = True
        
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
            
    def Update(self, pos, d_shaft_min, d_ext, length):
        self.pos_x = pos
        self.d_shaft_min = d_shaft_min
        self.d_ext = d_ext
        self.length = length
    
    def ExternalBearing(self, sign=1):
        B = self.B
        d = self.d
        D = self.D
        ep = min(0.1*D, 0.1*d)
        De = D + ep
        Dg, Dd = D, D
        if 'n' in self.connection_be:
            Dg = Dg - ep
        if 'p' in self.connection_be:
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
        if 'n' in self.connection_bi:
            dg = dg + ep
        if 'p' in self.connection_bi:
            dd = dd + ep
        bi = vm.Polygon2D([vm.Point2D((-B/2., sign*d/2.)), vm.Point2D((-B/2., sign*dg/2.)),
                           vm.Point2D((-B/2. - ep, sign*dg/2.)), vm.Point2D((-B/2. - ep, sign*di/2.)),
                           vm.Point2D((B/2. + ep, sign*di/2.)), vm.Point2D((B/2. + ep, sign*dd/2.)),
                           vm.Point2D((B/2., sign*dd/2.)), vm.Point2D((B/2., sign*d/2.))])
        return bi
    
    def BearingBox(self, sign=1):
        box = vm.Polygon2D([vm.Point2D((self.pos_x, sign*self.d_shaft_min/2.)),
                      vm.Point2D((self.pos_x, sign*self.d_ext/2.)),
                      vm.Point2D((self.pos_x + self.length, sign*self.d_ext/2.)),
                      vm.Point2D((self.pos_x + self.length, sign*self.d_shaft_min/2.))])
        return box
    
    def Plot(self, pos=0, a=None, box=True, typ='Graph'):
        """
        Generate a Plot
        
        :param pos: axial position of bearing assembly center
        :param a: object define the append graph (default value is None)
        :param box: draw the box parameter of the bearing assembly
        :param typ: define the aditionnal draw (default is None), 'Graph' draw the graph connection between bearing, 'Load' define the load
        """
        
        be_sup = self.ExternalBearing(sign = 1)
        be_inf = self.ExternalBearing(sign = -1)
        bi_sup = self.InternalBearing(sign = 1)
        bi_inf = self.InternalBearing(sign = -1)
        contour = [be_sup, be_inf, bi_sup, bi_inf]
        linkage_area = vm.Contour2D(contour)
        linkage_area = linkage_area.Translation((pos, 0), True)
        if a is None:
            f, a = linkage_area.MPLPlot(style = '-g')
        else:
            linkage_area.MPLPlot(a,'-g')
        
        contour = []
        pos_m = -self.B/2.
        for bg in self.list_bearing:
            cont = bg.Plot()
            cont = cont.Translation((pos_m + bg.B/2., 0), True)
            pos_m += bg.B
            contour.append(cont)
        assembly_bg = vm.Contour2D(contour)
        assembly_bg = assembly_bg.Translation((pos, 0), True)
        assembly_bg.MPLPlot(a,'-k')

        if typ == 'Graph':
            list_graph = []
            pos_m = -self.B/2.
            for bg in self.best_graph:
                graph = bg.PlotGraph()
                graph = graph.Translation((pos_m + bg.B/2., 0), True)
                pos_m += bg.B
                list_graph.append(graph)
            list_graph = vm.Contour2D(list_graph)
            list_graph = list_graph.Translation((pos, 0), True)
            list_graph.MPLPlot(a, '--b', True)
            
        elif typ == 'Load':
            pos_m = -self.B/2.
            for bg in self.best_graph:
                bg.PlotLoad(a, pos = pos + pos_m + bg.B/2.)
                pos_m += bg.B
        
        if box:
            box_sup = self.BearingBox(1)
            box_inf = self.BearingBox(-1)
            cont_box = [box_sup, box_inf]
            contour_box = vm.Contour2D(cont_box)
            contour_box = contour_box.Translation((pos, 0), True)
            contour_box.MPLPlot(a,'-r')
            
    def Graph(self, list_bearing):
        indice_bearing_x = []
        list_direction = []
        indice = 1
        for bg in list_bearing:
            if (bg.name in ['RadialBallBearing']) or (bg.name in ['RadialRollerBearing'] and bg.typ == 'NUP'):
                if max(indice_bearing_x + [0]) > 0:
                    if indice == indice_bearing_x[-1]:
                        indice_bearing_x.append(indice)
                    else:
                        indice += 1
                        indice_bearing_x.append(indice)
                else:
                    indice_bearing_x.append(indice)
                list_direction.append(0)
            else:
                indice_bearing_x.append(0)
                list_direction.append(bg.direction)
        list_name = [bg.name for bg in list_bearing]
                
        nb_bearing_x = max(indice_bearing_x)

        nb_bearing = len(list_bearing)
        case_graph = []
        if nb_bearing_x > 0:
            for list_sub_direction in product([-1, 1], repeat = nb_bearing_x):
                for ind_dir, sub_dir in enumerate(list_sub_direction):
                    pos_sub_dir = [pos_i for pos_i, i in enumerate(indice_bearing_x) if i == ind_dir + 1]
                    for pos_s in pos_sub_dir:
                        list_direction[pos_s] = sub_dir
                case_graph.append([*list_direction])
        else:
            case_graph.append(list_direction)
                
        graph = []
        li_case = []
        for case in case_graph:
            li = []
            for bg in list_bearing:
                bg_copy = deepcopy(bg)
                li.append(bg_copy)
            list_bg_export = self.ConnectionBearing(case, li)
            check = self.CheckViability(list_bg_export)
            if check:
                graph.append(deepcopy(list_bg_export))
                li_case.append(case)
        return graph, li_case
            
    def ConnectionBearing(self, case, list_bearing):
        compt = 0
        list_bg_export = []
        nb_bearing = len(list_bearing)
        list_name = [bg.name for bg in list_bearing]
        
        for i_direction, (direction, bg) in enumerate(zip(case, list_bearing)):
            list_nd = []
            if i_direction == 0:
                if 'n' in self.connection_be:
                    if bg.direction == 1:
                        list_nd.extend([LoadNode(compt), LoadNode(compt + 1)])
                    else:
                        list_nd.extend([LoadNode(compt, 0), LoadNode(compt + 1, 0)])
                    compt += 2
                else:
                    list_nd.extend([LoadNode(compt), LoadNode(compt + 1)])
                    compt += 2
                if 'n' in self.connection_bi:
                    if bg.direction == -1:
                        list_nd.extend([LoadNode(compt), LoadNode(compt + 1)])
                    else:
                        list_nd.extend([LoadNode(compt, 0), LoadNode(compt + 1, 0)])
                    compt += 2
                else:
                    list_nd.extend([LoadNode(compt), LoadNode(compt + 1)])
                    compt += 2
            elif nb_bearing > 1:
                list_nd = list_nd_m[4:8]
            if (nb_bearing > 1) and (i_direction < nb_bearing - 1):
                if list_name[i_direction] in ['RadialBallBearing', 'RadialRollerBearing']:
                    name = 'BallRoller'
                else:
                    name = 'Angular'
                if list_name[i_direction + 1] in ['RadialBallBearing', 'RadialRollerBearing']:
                    name_p = 'BallRoller'
                else:
                    name_p = 'Angular'
                if list_name[i_direction] == list_name[i_direction + 1] and direction == case[i_direction + 1]:
                    case_right = 1
                elif case[i_direction + 1] == -1 and (name == 'BallRoller' and name_p == 'Angular'):
                    case_right = 2
                elif case[i_direction + 1] == 1 and (name == 'BallRoller' and name_p == 'Angular'):
                    case_right = 3
                    
                elif direction == -1 and (name == 'Angular' and name_p == 'BallRoller'):
                    case_right = 3
                elif direction == 1 and (name == 'Angular' and name_p == 'BallRoller'):
                    case_right = 2
                
                elif (direction == -1 and case[i_direction + 1] == 1) and (name == 'Angular' and name_p == 'Angular'):
                    case_right = 3
                elif (direction == 1 and case[i_direction + 1] == -1) and (name == 'Angular' and name_p == 'Angular'):
                    case_right = 2
                    
                elif name == 'BallRoller' and name_p == 'BallRoller':
                    case_right = 1
                if case_right == 1:
                    list_nd.extend([LoadNode(i, 0) for i in range(compt, compt + 4)])
                    compt += 4
                if case_right == 2:
                    list_nd.extend([LoadNode(i, 0) for i in range(compt, compt + 2)])
                    compt += 2
                    list_nd.extend([LoadNode(i) for i in range(compt, compt + 2)])
                    compt += 2
                if case_right == 3:
                    list_nd.extend([LoadNode(i) for i in range(compt, compt + 2)])
                    compt += 2
                    list_nd.extend([LoadNode(i, 0) for i in range(compt, compt + 2)])
                    compt += 2
                    
            elif i_direction == nb_bearing - 1:
                if 'p' in self.connection_be:
                    if bg.direction == -1:
                        list_nd.extend([LoadNode(compt), LoadNode(compt + 1)])
                    else:
                        list_nd.extend([LoadNode(compt, 0), LoadNode(compt + 1, 0)])
                    compt += 2
                else:
                    list_nd.extend([LoadNode(compt), LoadNode(compt + 1)])
                    compt += 2
                if 'p' in self.connection_bi:
                    if bg.direction == 1:
                        list_nd.extend([LoadNode(compt), LoadNode(compt + 1)])
                    else:
                        list_nd.extend([LoadNode(compt, 0), LoadNode(compt + 1, 0)])
                    compt += 2
                else:
                    list_nd.extend([LoadNode(compt), LoadNode(compt + 1)])
                    compt += 2
            list_nd_m = [*list_nd]
            if (bg.name in ['RadialBallBearing']) or (bg.name in ['RadialRollerBearing'] and bg.typ == 'NUP'):
                bg.Graph(list_node = list_nd, num = i_direction, pos_x = i_direction, sub_direction = direction)
            else:
                bg.Graph(list_node = list_nd, num = i_direction, pos_x = i_direction)
#                b1 = copy.deepcopy(bg)
            list_bg_export.append(bg)
        return list_bg_export
    
    def CheckViability(self, graph):
        node_axial_ring = []
        if 'n' in self.connection_be:
            node_axial_ring.append(graph[0].list_node[0])
        if 'p' in self.connection_be:
            node_axial_ring.append(graph[-1].list_node[5])
        if 'n' in self.connection_bi:
            node_axial_ring.append(graph[0].list_node[2])
        if 'p' in self.connection_bi:
            node_axial_ring.append(graph[-1].list_node[7])
            
        node_axial_input = []
        for pos_bg, bg in enumerate(graph):
            if bg.name not in ['RadialRollerBearing', 'RadialBallBearing']:
                if bg.direction == 1:
                    node_axial_input.extend([bg.list_node[i] for i in [3, 4]])
                elif bg.direction == -1:
                    node_axial_input.extend([bg.list_node[i] for i in [1, 6]])
        
        valid = True
        for nd_axial_input in node_axial_input:
            valid_input = False
            for nd_axial_ring in node_axial_ring:
                valid_search = True
                n0 = [nd_axial_input]
                while valid_search:
                    np = []
                    for bg in graph:
                        for n1 in n0:
                            if n1 in list(bg.next.keys()):
                                if bg.next[n1] != []:
                                    for n1_add in bg.next[n1]:
                                        if n1_add.load is not None:
                                            np.append(n1_add)
                    n0 = np
                    if np == []:
                        valid_search = False
                    elif nd_axial_ring in np:
                        valid_search = False
                        valid_input = True
                        
            if valid_input == False:
                valid = False
        return valid
                
                            
    
    def BearingAssemblyLoad(self, fa, fr, results=None):
        if results is not None:
            results.Inputs(fa, fr)
            
        nb_radial_bearing = len([i for i in self.radial_load_linkage if i == True])
        fr_per_bearing = fr/(nb_radial_bearing*1.)
        
        for list_bearing in self.graph:
            # Load of angular and tapered bearing
            for ind_bg, bg in enumerate(list_bearing):
                if bg.name not in ['RadialRollerBearing', 'RadialBallBearing']:
                    if bg.direction == 1:
                        bg.list_node[3].ext_load = fr_per_bearing
                        bg.list_node[4].ext_load = fr_per_bearing
                    elif bg.direction == -1:
                        bg.list_node[1].ext_load = fr_per_bearing
                        bg.list_node[6].ext_load = fr_per_bearing
                    else:
                        continue
            # external axial load
            if fa > 0:
                list_bearing[0].list_node[3].ext_load = fa
            elif fa < 0:
                list_bearing[-1].list_node[6].ext_load = -fa
            # axial load solver
            
            valid = True
            list_valid = [0]
            while valid:
                list_valid_m =list_valid
                list_valid = []
                for bg in list_bearing:
                    bg.CheckNone()
                    bg.PropagateNone()
                    bg.CheckNone()
                    bg.PropagateNone()
                    bg.CheckNone()
                for bg in list_bearing:
                    bg.Load()
                for bg in list_bearing:
                    bg.ResumeLoad()
                    for nd in bg.list_node:
                        list_valid.append(nd.load)
                if list_valid == list_valid_m:
                    valid = False
        # search the best sub graph
        indice_m = 0
        best_graph = 0
        for best, list_bearing in enumerate(self.graph):
            indice = 0
            n0 = list_bearing[0].list_node[0].load
            n2 = list_bearing[0].list_node[2].load
            n5 = list_bearing[-1].list_node[5].load
            n7 = list_bearing[-1].list_node[7].load
            for n in [n0, n2, n5, n7]:
                if (n is not None) and (n != 0):
                    indice += 1
            if indice > indice_m:
                best_graph = best
                indice_m = indice
        self.best_graph = self.graph[best_graph]
        self.case = self.li_case[best_graph]
        #resume of the axial load
        n0 = self.best_graph[0].list_node[0].load
        n2 = self.best_graph[0].list_node[2].load
        n5 = self.best_graph[-1].list_node[5].load
        n7 = self.best_graph[-1].list_node[7].load
        fa = 0
        for n in [n0, n2, n5, n7]:
            if n is not None:
                fa += n
        
        if results is not None:
            results.Outputs(fa)
            return fa, results
        else:
            return fa
    
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
            
        del d['list_bearing']
        li_bg = []
        for bg in self.best_graph:
            li_bg.append(bg.Dict())
        d['list_bearing'] = li_bg
        del d['best_graph']
            
        del d['graph']
        return d
    
    @classmethod
    def Dict2Obj(cls, d):
        li_bg = []
        for li in d['list_bearing']:
            if li['name'] == 'RadialRollerBearing':
                BA = ConceptRadialRollerBearing.Dict2Obj(li)
            elif li['name'] == 'TaperedRollerBearing':
                BA = ConceptTaperedRollerBearing.Dict2Obj(li)
            elif li['name'] == 'SphericalBallBearing':
                BA = ConceptSphericalBallBearing.Dict2Obj(li)
            elif li['name'] == 'AngularBallBearing':
                BA = ConceptAngularBallBearing.Dict2Obj(li)
            elif li['name'] == 'RadialBallBearing':
                BA = ConceptRadialBallBearing.Dict2Obj(li)
            li_bg.append(BA)
        obj = cls(list_bearing = li_bg, radial_load_linkage = d['radial_load_linkage'], 
                  internal_pre_load = 0, connection_bi = d['connection_bi'], 
                  connection_be = d['connection_be'], behavior_link = d['behavior_link'])
        obj.best_graph = obj.ConnectionBearing(case = d['case'], list_bearing = li_bg)
        return obj
            
class CompositeBearingAssembly:
    def __init__(self, list_bearing_assembly, list_position=None, pre_load=0, pos_x=None):
        
        self.list_bearing_assembly = list_bearing_assembly
        self.mass = self.MassCalculate()
        if pos_x is not None:
            pos_x = list(pos_x)
        self.pos_x = pos_x
        
    def Update(self, pos_x, list_pos_unknown, list_load, d_shaft_min, axial_pos, 
               d_ext, length):
        if pos_x is not None:
            pos_x = list(pos_x)
        self.pos_x = pos_x
        self.list_pos_unknown = list_pos_unknown
        self.list_load = list_load
        self.d_shaft_min = d_shaft_min
        self.axial_pos = axial_pos
        self.d_ext = d_ext
        self.length = length
        for num_linkage, assembly_bg in enumerate(self.list_bearing_assembly):
            pos = self.axial_pos[num_linkage] - self.pos_x[num_linkage]
            assembly_bg.Update(pos, self.d_shaft_min[num_linkage], self.d_ext[num_linkage],
                               self.length[num_linkage])
        
    def MassCalculate(self):
        mass = 0
        for li_bg in self.list_bearing_assembly:
            mass += li_bg.mass
        return mass
    
    def Shaft(self):
        
        d1 = self.list_bearing_assembly[0].d
        d2 = self.list_bearing_assembly[1].d
        B1 = self.list_bearing_assembly[0].B
        B2 = self.list_bearing_assembly[1].B
        pos_mid = (self.pos_x[0] + self.pos_x[1])/2.
        shaft = vm.Polygon2D([vm.Point2D((self.pos_x[0] - B1/2. - 5e-3, -d1/2.)),
                      vm.Point2D((self.pos_x[0] - B1/2. - 5e-3, d1/2.)),
                      vm.Point2D((pos_mid, d1/2.)),
                      vm.Point2D((pos_mid, d2/2.)),
                      vm.Point2D((self.pos_x[1] + B2/2. + 5e-3, d2/2.)),
                      vm.Point2D((self.pos_x[1] + B2/2. + 5e-3, -d2/2.)),
                      vm.Point2D((pos_mid, -d2/2.)),
                      vm.Point2D((pos_mid, -d1/2.)),])
        return shaft
    
    def Plot(self, box=True, typ=None):
        
        shaft = self.Shaft()
        contour_shaft = vm.Contour2D([shaft])
        f, a = contour_shaft.MPLPlot(style = '-k')
        
        for assembly_bg, pos in zip(self.list_bearing_assembly, self.pos_x):
            assembly_bg.Plot(pos, a, box, typ)
    def Graph(self):
        for li_bg in self.list_bearing_assembly:
            G, positions, li_axial_link, nd_axial_load = li_bg.Graph()
            plt.figure()
            nx.draw_networkx(G, pos = positions)
    
    def ShaftLoad(self, pos1, pos2, list_pos_unknown, list_load, list_torque,
                  results=None):
        
        if results is not None:
            results.Inputs(pos1, pos2, list_pos_unknown, list_load, list_torque)
            li_results_bearing_assembly = []
        
        ground = genmechanics.Part('ground')
        shaft1 = genmechanics.Part('shaft1')
        p1 = npy.array([pos1, 0, 0])
        p2 = npy.array([pos2, 0, 0])
        bearing1 = linkages.FrictionlessBallLinkage(ground,shaft1,p1,[0,0,0],'bearing1')
        bearing2 = linkages.FrictionlessLinearAnnularLinkage(ground,shaft1,p2,[0,0,0],'bearing2')
        
        load1 = []
        for pos, ld, tq in zip(list_pos_unknown, list_load, list_torque):
            load1.append(loads.KnownLoad(shaft1, pos, [0,0,0], ld, tq, 'input'))
        load2 = loads.SimpleUnknownLoad(shaft1, [(pos1 + pos2)/2,0,0], [0,0,0], [], [0], 'output torque')
        imposed_speeds = [(bearing1, 0, 100)]
        
        mech = genmechanics.Mechanism([bearing1,bearing2],ground,imposed_speeds,load1,[load2])
        
        axial_load = 0
        for load in list_load:
            axial_load += load[0]
        
        li_fa = []
        li_fr = []
        for li_bg, bg in zip(self.list_bearing_assembly, [bearing1, bearing2]):
            tensor = mech.GlobalLinkageForces(bg,1)
            fr = (tensor[1]**2 + tensor[2]**2)**(0.5)
            if results is not None:
                re_bearing_assembly = ResultsBearingAssembly()
                fa, re_bearing_assembly = li_bg.BearingAssemblyLoad(fa = axial_load, fr = fr, 
                                               results = re_bearing_assembly)
                li_results_bearing_assembly.append(re_bearing_assembly)
            else:
                fa = li_bg.BearingAssemblyLoad(fa = axial_load, fr = fr)
            li_fa.append(fa)
            li_fr.append(fr)
        
        if results is not None:
            results.Outputs(li_fa, li_fr, li_results_bearing_assembly)
            return results
        else:
            return li_fa, li_fr
        
    
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
                
        li_bg = []
        for item in self.list_bearing_assembly:
            li_bg.append(item.Dict())
        d['list_bearing_assembly'] = li_bg
            
        return d
        
    @classmethod
    def Dict2Obj(cls, d):
        li_bg = []
        for li in d['list_bearing_assembly']:
            BA = BearingAssembly.Dict2Obj(li)
            li_bg.append(BA)
                 
        obj = cls(list_bearing_assembly = li_bg, list_position = None, pre_load = 0, 
                  pos_x = d['pos_x'])
        return obj
    
#class ResultsBearing:
#    
class ResultsBearingAssembly:
    def __init__(self):
        self.inputs = []
        self.outputs = []
        
    def Inputs(self, fa, fr):
        inputs = {}
        inputs['fa'] = fa
        inputs['fr'] = fr
        self.inputs.append(inputs)
        
    def Outputs(self, fa):
        self.outputs.append({'fa': fa})
    
class ResultsCompositeBearingAssembly:
    def __init__(self):
        self.inputs = []
        self.outputs = []
        
    def Inputs(self, pos1, pos2, list_pos_unknown, list_load, list_torque):
        inputs = {}
        inputs['pos1'] = pos1
        inputs['pos2'] = pos2
        inputs['list_pos_unknown'] = list_pos_unknown
        inputs['list_load'] = list_load
        inputs['list_torque'] = list_torque
        self.inputs.append(inputs)
        
    def Outputs(self,li_fa, li_fr, results_bearing_assembly):
        self.outputs.append({'li_fa': li_fa, 'li_fr': li_fr, 
                             'results_bearing_assembly': results_bearing_assembly})

class RadialRollerBearing(ConceptRadialRollerBearing):
    #Roulement à rouleaux
    def __init__(self, d, D, B, i, Z, Dw, d1=None, D1=None, Lw=None, radius=None, E=None, F=None, 
                 alpha=0, bm=1.1, oil=oil_iso_vg_1500, material=material_iso, direction=1, typ='N'):
        ConceptRadialRollerBearing.__init__(self, d, D, B, i, Z, Dw, alpha ,
                                oil = oil, material = material, typ_contact = 'linear_contact',
                                direction = direction, typ = typ)

        #diametre rouleau moyen
        self.Dwe = Dw
        self.bm = bm
        
        if d1 is not None:
            self.d1 = d1
        if D1 is not None:
            self.D1 = D1
        if Lw is not None:
            self.Lw = Lw
        if E is not None:
            self.E = E
        if F is not None:
            self.F = F
        if radius is not None:
            self.radius = radius
                
        self.Dpw, self.Lwe, self.jeu, self.ep = self.DefParam()
        self.mass = self.Mass()
        
    def __str__(self):
        s = '{}\n'.format(self.__class__.__name__)
        s += 'Dimensions d:{} D:{} B:{} Dw:{} Z:{}'.format(self.d, self.D, self.B, self.Dw, self.Z)
#        s += 'L10: {} Lnm: {}
        return s
        
    def Update(self, d1, D1, E, F, Z):
        self.d1 = d1
        self.D1 = D1
        self.E = E
        self.F = F
        self.Z = Z
        self.Dpw,self.Lwe,self.jeu,self.ep = self.DefParam()
        self.mass = self.Mass()
        
    def DefParam(self):
        Dpw = (self.E+self.F)/2.
        Lwe = self.Lw-2*self.radius
        jeu = (self.E-self.F-2*self.Dw)/4.
        ep = (self.B-self.Lw-2*jeu)/2.
        return Dpw,Lwe,jeu,ep
        
    def Mass(self):
        rho = 7800
        m = self.Z*npy.pi*(self.Dw)**2/4.*self.Lw*rho
        m += (npy.pi*(self.D)**2/4.-npy.pi*(self.E)**2/4.)*self.B*rho
        m += (npy.pi*(self.F)**2/4.-npy.pi*(self.d)**2/4.)*self.B*rho
        return m
    
    def BaseStaticLoad(self):
        #Charge radiale statique de base
        #besoin de convertir les dimensions en mm pour les formules ISO
        C0r = 44*(1-(self.Dw*1e3*npy.cos(self.alpha))/(self.Dpw*1e3))*self.i*self.Z \
                *self.Lwe*1e3*self.Dw*1e3*npy.cos(self.alpha)
        return C0r
    
    def BaseDynamicLoad(self):
        #Charge radiale dynamique de base
        mu = float((self.Dwe*1e3)*npy.cos(self.alpha)/(self.Dpw*1e3))
        fc = 0.377*self.material.mu_delta*1/((2**((self.material.weibull_c\
                +self.material.weibull_h-1)/(self.material.weibull_c-self.material.weibull_h+1)))\
                *(0.5**(2*self.material.weibull_e/(self.material.weibull_c-self.material.weibull_h+1))))\
                *self.material.B1*((1-mu)**((self.material.weibull_c+self.material.weibull_h-3)/\
                (self.material.weibull_c-self.material.weibull_h+1))/((1+mu)**(2*self.material.weibull_e\
                /(self.material.weibull_c-self.material.weibull_h+1))))\
                *(mu**(2/(self.material.weibull_c-self.material.weibull_h+1)))\
                *(1+(1.04*((1-mu)/(1+mu))**((self.material.weibull_c+self.material.weibull_h\
                +2*self.material.weibull_e-3)/(self.material.weibull_c-self.material.weibull_h+1)))\
                **((self.material.weibull_c-self.material.weibull_h+1)/2.))\
                **(-2/(self.material.weibull_c-self.material.weibull_h+1))
        Cr = fc*self.bm*self.i*((self.Lwe*1e3)*npy.cos(self.alpha))\
                **((self.material.weibull_c-self.material.weibull_h-1)/(self.material.weibull_c\
                -self.material.weibull_h+1))*self.Z**((self.material.weibull_c\
                -self.material.weibull_h-2*self.material.weibull_e+1)/(self.material.weibull_c\
                -self.material.weibull_h+1))*(self.Dwe*1e3)**((self.material.weibull_c\
                -self.material.weibull_h-3)/(self.material.weibull_c-self.material.weibull_h+1))
        return Cr

    def Dict(self):

        d = {}
        for k,v in self.__dict__.items():
            tv = type(v)
            if tv == npy.int64:
                d[k] = int(v)
            elif tv == npy.float64:
                d[k] = float(v)
            else:
                d[k] = v

        d['oil'] = self.oil.Dict()
        d['material'] = self.material.Dict()
        return d

class DrawnCupNeedleRollerBearing(RadialRollerBearing):
    #Douille à aiguilles

    def __init__(self, typ, B, d, D, d1, D1, Lw, Dw, radius, E, F, Z, i,
                 alpha,bm=1, weibull_e=9/8., weibull_c=31/3., weibull_h=7/3.,
                 B1=551.13373/0.483, mu_delta=0.83, c_gamma=0.05,
                 oil_name='iso_vg_100'):
        RadialRollerBearing.__init__(self, typ, B, d, D, d1, D1, Lw, Dw, radius, E,
                                     F, Z, i, alpha, bm, weibull_e, weibull_c,
                                     weibull_h, B1, mu_delta, c_gamma, oil_name)
        
class NeedleRollerBearing(RadialRollerBearing):
    #Cage à aiguilles
    def __init__(self, typ, B, d, D, d1, D1, Lw, Dw, radius, E, F, Z, i,
                 alpha, bm=1, weibull_e=9/8., weibull_c=31/3., weibull_h=7/3.,
                 B1=551.13373/0.483, mu_delta=0.83, c_gamma=0.05,
                 oil_name='iso_vg_100'):
        RadialRollerBearing.__init__(self, typ, B, d, D, d1, D1, Lw, Dw, radius, E,
                                     F, Z, i, alpha, bm, weibull_e, weibull_c,
                                     weibull_h, B1, mu_delta, c_gamma, oil_name)

class SphericalRollerBearing(RadialRollerBearing):
    #Roulement à rotule à rouleaux
    def __init__(self, typ, B, d, D, d1, D1, Lw, Dw, radius, E, F, Z, i,
                 alpha,bm=1.15, weibull_e=9/8., weibull_c=31/3., weibull_h=7/3.,
                 B1=551.13373/0.483, mu_delta=0.83, c_gamma=0.05,
                 oil_name='iso_vg_100'):
        RadialRollerBearing.__init__(self, typ, B, d, D, d1, D1, Lw, Dw, radius, E,
                                     F, Z, i, alpha, bm, weibull_e, weibull_c,
                                     weibull_h, B1, mu_delta, c_gamma,
                                     oil_name)