import numpy as npy
#import math as mt
from scipy import interpolate
#import os
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math
from copy import deepcopy

from scipy.optimize import minimize,fsolve
import networkx as nx
import matplotlib.pyplot as plt

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
        return self.__dict__

material_iso=Material()

class RadialBearing:
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
        
    def Graph(self, G=None, li_nd=[None, None, None, None, None, None, None, None]):

        nb_bearing = 0
        if G is not None:
            for tup in list(G.edges(data = True)):
                if 'bearing' in tup[2].keys():
                    nb_bearing = max(nb_bearing, tup[2]['bearing'])
        
        if G is None:
            G = nx.DiGraph()
            nd_m = 0
        else:
            nd_m = max([a for a in li_nd if a is not None] + [max(G.nodes)])
        li = li_nd
        for pos, nd in enumerate(li_nd):
            if nd is None:
                li[pos] = nd_m + 1
                nd_m = li[pos]
        G.add_edges_from([(li[4], li[0])], typ = 'ring', bearing = nb_bearing + 1)
        G.add_edges_from([(li[1], li[5])], typ = 'ring', bearing = nb_bearing + 1)
        G.add_edges_from([(li[6], li[2])], typ = 'ring', bearing = nb_bearing + 1)
        G.add_edges_from([(li[3], li[7])], typ = 'ring', bearing = nb_bearing + 1)

        return G, li

    
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
    
    def Graph(self, G=None, li_nd=[None, None, None, None, None, None, None, None], pos_x=0, sub_graph=None, copy=False):
        
        if copy:
            G = deepcopy(G)
        G, li = RadialBearing.Graph(self, G, li_nd)
        
        nb_bearing = 0
        if G is not None:
            for tup in list(G.edges(data = True)):
                if 'bearing' in tup[2].keys():
                    nb_bearing = max(nb_bearing, tup[2]['bearing'])
                    
        if sub_graph == 1:
            G.add_edges_from([(li[4], li[2])], typ = 'internal', bearing = nb_bearing)
            G.add_edges_from([(li[3], li[5])], typ = 'internal', bearing = nb_bearing)
        elif sub_graph == -1:
            G.add_edges_from([(li[6], li[0])], typ = 'internal', bearing = nb_bearing)
            G.add_edges_from([(li[1], li[7])], typ = 'internal', bearing = nb_bearing)
        else:
            G.add_edges_from([(li[4], li[2])], typ = 'internal', bearing = nb_bearing)
            G.add_edges_from([(li[3], li[5])], typ = 'internal', bearing = nb_bearing)
            G.add_edges_from([(li[6], li[0])], typ = 'internal', bearing = nb_bearing)
            G.add_edges_from([(li[1], li[7])], typ = 'internal', bearing = nb_bearing)
        decal_y = 0.1
        li_pos = {li[0]: (pos_x, 1), li[1]: (pos_x, 1 - decal_y), 
                  li[2]: (pos_x, decal_y), li[3]: (pos_x, 0),
                  li[4]: (pos_x + 1, 1), li[5]: (pos_x + 1, 1 - decal_y), 
                  li[6]: (pos_x + 1, decal_y), li[7]: (pos_x + 1, 0)}
        return G, li, li_pos
        
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
    
    def Graph(self, G=None, li_nd=[None, None, None, None, None, None, None, None], pos_x=0):
        
        nb_bearing = 0
        if G is not None:
            for tup in list(G.edges(data = True)):
                if 'bearing' in tup[2].keys():
                    nb_bearing = max(nb_bearing, tup[2]['bearing'])
                    
        G, li = RadialBearing.Graph(self, G, li_nd)
        if len(li) == 8:
            li.extend([max(li) + 1, max(li) + 2, max(li) + 3])
        decal_x = 0.3
        decal_y = 0.1
    
        if self.direction == 1:
            G.add_edges_from([(li[4], li[8])], typ = 'internal', bearing = nb_bearing)
            G.add_edges_from([(li[8], li[2])], typ = 'internal', bearing = nb_bearing)
            G.add_edges_from([(li[3], li[9])], typ = 'internal', bearing = nb_bearing)
            G.add_edges_from([(li[9], li[5])], typ = 'internal', bearing = nb_bearing)
            G.add_edges_from([(li[10], li[8])], typ = 'internal', bearing = nb_bearing)
            G.add_edges_from([(li[10], li[9])], typ = 'internal', bearing = nb_bearing)
            li_pos = {li[0]: (pos_x + decal_x, 1), li[1]: (pos_x + decal_x, 1 - decal_y), 
              li[2]: (pos_x, decal_y), li[3]: (pos_x, 0),
              li[4]: (pos_x + 1, 1), li[5]: (pos_x + 1, 1 - decal_y), 
              li[6]: (pos_x + 1 - decal_x, decal_y), li[7]: (pos_x + 1 - decal_x, 0)}
            li_pos[li[8]] = ((li_pos[li[4]][0] + li_pos[li[2]][0])/2., (li_pos[li[4]][1] + li_pos[li[2]][1])/2.)
            li_pos[li[9]] = ((li_pos[li[5]][0] + li_pos[li[3]][0])/2., (li_pos[li[5]][1] + li_pos[li[3]][1])/2.)
            li_pos[li[10]] = ((li_pos[li[3]][0] + li_pos[li[7]][0])/2., -0.5)
        elif self.direction == -1:
            G.add_edges_from([(li[1], li[8])], typ = 'internal', bearing = nb_bearing)
            G.add_edges_from([(li[8], li[7])], typ = 'internal', bearing = nb_bearing)
            G.add_edges_from([(li[6], li[9])], typ = 'internal', bearing = nb_bearing)
            G.add_edges_from([(li[9], li[0])], typ = 'internal', bearing = nb_bearing)
            G.add_edges_from([(li[10], li[8])], typ = 'internal', bearing = nb_bearing)
            G.add_edges_from([(li[10], li[9])], typ = 'internal', bearing = nb_bearing)
            li_pos = {li[0]: (pos_x, 1), li[1]: (pos_x, 1 - decal_y), 
              li[2]: (pos_x + decal_x, decal_y), li[3]: (pos_x + decal_x, 0),
              li[4]: (pos_x + 1 - decal_x, 1), li[5]: (pos_x + 1 - decal_x, 1 - decal_y), 
              li[6]: (pos_x + 1, decal_y), li[7]: (pos_x + 1, 0)}
            li_pos[li[8]] = ((li_pos[li[1]][0] + li_pos[li[7]][0])/2., (li_pos[li[1]][1] + li_pos[li[7]][1])/2.)
            li_pos[li[9]] = ((li_pos[li[6]][0] + li_pos[li[0]][0])/2., (li_pos[li[6]][1] + li_pos[li[0]][1])/2.)
            li_pos[li[10]] = ((li_pos[li[3]][0] + li_pos[li[7]][0])/2., -0.5)
        return G, li, li_pos
            
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
    
    def Graph(self, G=None, li_nd=[None, None, None, None, None, None, None, None], 
              pos_x=0, sub_graph=None, copy=False):
        
        if copy:
            G = deepcopy(G)
        G, li = RadialBearing.Graph(self, G, li_nd)
        
        nb_bearing = 0
        if G is not None:
            for tup in list(G.edges(data = True)):
                if 'bearing' in tup[2].keys():
                    nb_bearing = max(nb_bearing, tup[2]['bearing'])
                    
        if self.typ in ['NUP']:
            if sub_graph == 1:
                G.add_edges_from([(li[4], li[2])], typ = 'internal', bearing = nb_bearing)
                G.add_edges_from([(li[3], li[5])], typ = 'internal', bearing = nb_bearing)
            elif sub_graph == -1:
                G.add_edges_from([(li[6], li[0])], typ = 'internal', bearing = nb_bearing)
                G.add_edges_from([(li[1], li[7])], typ = 'internal', bearing = nb_bearing)
            else:
                G.add_edges_from([(li[4], li[2])], typ = 'internal', bearing = nb_bearing)
                G.add_edges_from([(li[3], li[5])], typ = 'internal', bearing = nb_bearing)
                G.add_edges_from([(li[6], li[0])], typ = 'internal', bearing = nb_bearing)
                G.add_edges_from([(li[1], li[7])], typ = 'internal', bearing = nb_bearing)
        elif (self.typ == 'NJ' and self.direction == 1) or (self.typ == 'NF' and self.direction == -1):
            G.add_edges_from([(li[4], li[2])], typ = 'internal', bearing = nb_bearing)
            G.add_edges_from([(li[3], li[5])], typ = 'internal', bearing = nb_bearing)
        elif (self.typ == 'NJ' and self.direction == -1) or (self.typ == 'NF' and self.direction == 1):
            G.add_edges_from([(li[6], li[0])], typ = 'internal', bearing = nb_bearing)
            G.add_edges_from([(li[1], li[7])], typ = 'internal', bearing = nb_bearing)
            
        decal_y = 0.1
        li_pos = {li[0]: (pos_x, 1), li[1]: (pos_x, 1 - decal_y), 
                  li[2]: (pos_x, decal_y), li[3]: (pos_x, 0),
                  li[4]: (pos_x + 1, 1), li[5]: (pos_x + 1, 1 - decal_y), 
                  li[6]: (pos_x + 1, decal_y), li[7]: (pos_x + 1, 0)}
        return G, li, li_pos
    
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
    
    def Graph(self, G=None, li_nd=[None, None, None, None, None, None, None, None], pos_x=0):
        
        G, li, li_pos = ConceptAngularBallBearing.Graph(self, G, li_nd, pos_x)
        return G, li, li_pos
    
class BearingAssembly:
    def __init__(self, list_bearing, radial_load_linkage, internal_pre_load=0, 
                 connection_bi=['n', 'p'], connection_be=['n', 'p'], behavior_link='pn'):
        self.list_bearing = list_bearing
        self.radial_load_linkage = radial_load_linkage
        self.connection_be = connection_be
        self.connection_bi = connection_bi
        self.behavior_link = behavior_link
        self.mass = 0
        for bg in self.list_bearing:
            if bg.mass is not None:
                self.mass += bg.mass
        self.B = 0
        for bg in self.list_bearing:
            self.B += bg.B
        self.D = 0
        for bg in self.list_bearing:
            self.D = max(self.D, bg.D)
        self.d = npy.inf
        for bg in self.list_bearing:
            self.d = min(self.d, bg.d)
        
        self.graph, self.positions, self.li_axial_link, self.nd_axial_load = self.Graph()
        
        li_link = []
        if 'n' in self.connection_be:
            li_link.append(self.li_axial_link[0])
        if 'p' in self.connection_be:
            li_link.append(self.li_axial_link[5])
        if 'n' in self.connection_bi:
            li_link.append(self.li_axial_link[2])
        if 'p' in self.connection_bi:
            li_link.append(self.li_axial_link[7])
        self.li_link = li_link
            
        try:
            check = self.CheckViability(self.li_link, self.nd_axial_load)
        except:
            self.PlotGraph()
            raise KeyError('No solution found for bearing mounting')
        
    def PlotGraph(self):
        for key_graph, val_graph in self.graph.items():
            plt.figure()
            nx.draw_networkx(val_graph, pos = self.positions)
            
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
    
    def Plot(self, pos=0, a=None, box=True):
        
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
        
        if box:
            box_sup = self.BearingBox(1)
            box_inf = self.BearingBox(-1)
            cont_box = [box_sup, box_inf]
            contour_box = vm.Contour2D(cont_box)
            contour_box = contour_box.Translation((pos, 0), True)
            contour_box.MPLPlot(a,'-r')
    
    def Graph(self):
        G = {}
        bg_0 = self.list_bearing[0]
        if (bg_0.name == 'RadialBallBearing') or (bg_0.name == 'RadialRollerBearing' and bg_0.typ == 'NUP'):
            G[(1,)], li, li_pos = bg_0.Graph(li_nd=[None, None, None, None, None, None, None, None], sub_graph = 1)
            G[(-1,)], li, li_pos = bg_0.Graph(li_nd=[None, None, None, None, None, None, None, None], sub_graph = -1)
        else:
            G[(0,)], li, li_pos = bg_0.Graph(li_nd=[None, None, None, None, None, None, None, None])
        positions = li_pos
        li_axial_link = li[0:8]
        nd_axial_load = {}
        if len(li) > 8:
            nd_axial_load[0] = li[-3:-1]
            
        for pos, bg in enumerate(self.list_bearing[1:]):
            if (bg.name == 'RadialBallBearing') or (bg.name == 'RadialRollerBearing' and bg.typ == 'NUP'):
                double_graph = True
            else:
                double_graph = False
                
            if bg_0.name in ['RadialRollerBearing', 'RadialBallBearing']:
                if bg.name in ['RadialRollerBearing', 'RadialBallBearing']:
                    li_m = li[4: 8] + [None]*4
                elif bg.direction == 1:
                    li_m = [None]*2 + li[6: 8] + [None]*4
                elif bg.direction == -1:
                    li_m = li[4: 6] + [None]*6
            elif bg_0.direction == 1:
                if bg.name in ['RadialRollerBearing', 'RadialBallBearing']:
                    li_m = li[4: 6] + [None]*6
                elif bg.direction == 1:
                    li_m = li[4: 8] + [None]*4
                elif bg.direction == -1:
                    li_m = li[4: 6] + [None]*6
            elif bg_0.direction == -1:
                if bg.name in ['RadialRollerBearing', 'RadialBallBearing']:
                    li_m = [None]*2 + li[6: 8] + [None]*4
                elif bg.direction == 1:
                    li_m = [None]*2 + li[6: 8] + [None]*4
                elif bg.direction == -1:
                    li_m = li[4: 8] + [None]*4
            Gp = {}
            for key_G, value_G in G.items():
                if double_graph:
                    if key_G[-1] == 1:
                        add = *key_G,1
                        Gp[add], li, li_pos = bg.Graph(value_G, li_m, pos + 1, sub_graph = 1)
                    elif key_G[-1] == -1:
                        add = *key_G,-1
                        Gp[add], li, li_pos = bg.Graph(value_G, li_m, pos + 1, sub_graph = -1)
                    else:
                        add = *key_G,1
                        Gp[add], li, li_pos = bg.Graph(value_G, li_m, pos + 1, sub_graph = 1, copy = True)
                        add = *key_G,-1
                        Gp[add], li, li_pos = bg.Graph(value_G, li_m, pos + 1, sub_graph = -1, copy = True)
                else:
                    add = *key_G,0
                    Gp[add], li, li_pos = bg.Graph(value_G, li_m, pos + 1)
            G = dict(list(Gp.items()))
            li_axial_link[4:8] = li[4:8]
            if len(li) > 8:
                nd_axial_load[pos + 1] = li[-3:-1]
            positions = dict(list(li_pos.items()) + list(positions.items()))
            bg_0 = bg
        return G, positions, li_axial_link, nd_axial_load
    
    def CheckViability(self, li_link, nd_axial_load):
        
        valid = False
        valid_graph = {}
        for key_graph, val_graph in self.graph.items():
            short_path = []
            valid_sub_graph = True
            for pos, li_nd in nd_axial_load.items():
                for nd in li_nd:
                    valid_iter = False
                    for bg in li_link:
                        if nx.has_path(val_graph, nd, bg):
                            short_path.append([p for p in nx.all_shortest_paths(val_graph, 
                                              source = nd, target = bg)])
                            valid_iter = True
                    if not valid_iter:
                        valid_sub_graph = False
            if valid_sub_graph:
                valid = True
                g = nx.DiGraph()
                for path in short_path:
                    for pt in path:
                        g.add_path(pt)
                valid_graph[val_graph] = g
#                plt.figure()
#                nx.draw_networkx(g, pos = self.positions)
        return valid
    
    def BearingAssemblyLoad(self, fa, fr):
        nb_radial_bearing = len([i for i in self.radial_load_linkage if i == True])
        fr_per_bearing = fr/(nb_radial_bearing*1.)
        
#        for key_graph, val_graph in self.graph.items():
#            for pos, li_nd in self.nd_axial_load.items():
#                for nd in li_nd:
#                    for bg in self.li_link:
#                        if nx.has_path(val_graph, nd, bg):
#                            short_path = [p for p in nx.all_shortest_paths(val_graph, 
#                                              source = nd, target = bg)]
#                            g = nx.DiGraph()
#                            for path in short_path:
#                                g.add_path(path)
#                            for node1, node2 in list(nx.dfs_edges(g, nd)):
#                                n1 = g.degree[node1]
#                                n2 = g.degree[node2]
#                                print(n1, n2)
#                            plt.figure()
#                            nx.draw_networkx(g)
#                            return g
            
#            g_nv = nx.DiGraph()
#            for (node1, node2) in self.nd_axial_load.values():
#                g1 = nx.dfs_edges(sub_gr, node1)
#                g2 = nx.dfs_edges(sub_gr, node2)
#                
#                g_nv.add_edges_from(list(g1))
#                g_nv.add_edges_from(list(g2))
#            plt.figure()
#            nx.draw_networkx(g_nv)
    
            
class CompositiveBearingAssembly:
    def __init__(self, list_bearing_assembly, list_position=None, pre_load=0, pos_x=None):
        
        self.list_bearing_assembly = list_bearing_assembly
        self.mass = self.MassCalculate()
        self.pos_x = pos_x
        
    def Update(self, pos_x, list_pos_unknown, list_load, d_shaft_min, axial_pos, 
               d_ext, length):
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
    
    def Plot(self):
        
        shaft = self.Shaft()
        contour_shaft = vm.Contour2D([shaft])
        f, a = contour_shaft.MPLPlot(style = '-k')
        
        for assembly_bg, pos in zip(self.list_bearing_assembly, self.pos_x):
            assembly_bg.Plot(pos, a)
            
    def Graph(self):
        for li_bg in self.list_bearing_assembly:
            G, positions, li_axial_link, nd_axial_load = li_bg.Graph()
            plt.figure()
            nx.draw_networkx(G, pos = positions)
    
    def ShaftLoad(self, pos1, pos2, list_pos_unknown, list_load, list_torque):
        
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
            fa = li_bg.BearingAssemblyLoad(fa = 0, fr = fr)
            li_fa.append(fa)
            li_fr.append(fr)
        

        return li_fa, li_fr
        

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