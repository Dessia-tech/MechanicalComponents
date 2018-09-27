import numpy as npy
#import math as mt
from scipy import interpolate
#import os
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math
from scipy.optimize import minimize,fsolve

from mechanical_components.bearings_snr import RadialRollerBearingSNR
from mechanical_components.catalogs.dico_bearings_ISO \
    import dico_rlts_iso,dico_roller_iso,dico_radial_clearance_iso,dico_rules


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
    def __init__(self, d, D, i, Z, Dw, alpha, Cr=None, C0r=None ,oil=oil_iso_vg_1500, 
                 material=material_iso, typ_contact=None):
        self.d = d
        self.D = D
        self.Dpw = (d + D)/2.
        self.i = i
        self.Z = Z
        self.Dw = Dw
        self.alpha = alpha
        self.oil = oil
        self.material = material
        self.typ_contact = typ_contact
        if Cr is not None:
            self.Cr = Cr
        if C0r is not None:
            self.C0r = C0r
    
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
        
        Pr = (Pr / total_cycles) ** (1/self.coeff_baselife)
#        Cr=self.BaseDynamicLoad()

        L10=(Cr/Pr)**(self.coeff_baselife)
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
            if (fr != 0.) or (fa != 0.):
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
        return total_cycles / nci_Lpi
    
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
    def __init__(self, d, D, i, Z, Dw, alpha, Cr=None, C0r=None ,oil=oil_iso_vg_1500, 
                 material=material_iso, typ_contact=None):
        RadialBearing.__init__(self, d, D, i, Z, Dw, alpha, Cr, C0r, oil, material, typ_contact)
        self.coeff_baselife = 3.
        
    def EquivalentStaticLoad(self, fr, fa=None):
        #Charge radiale statique équivalente
        X0 = 0.6
        Y0 = 0.5
        P0r = max(fr,X0*fr+Y0*fa)
        return P0r
    
    def EquivalentDynamicLoad(self, fr, fa = 0):
        alphap = fsolve((lambda alphap:npy.cos(5/180*npy.pi)/npy.cos(alphap) \
                    -(1+0.012534*(fa/(self.i*self.Z*((self.Dw*1e3)**2)*npy.sin(alphap)))**(2/3.))),self.alpha)[0]
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
        
class ConceptAngularBallBearing(RadialBearing):
    def __init__(self, d, D, i, Z, Dw, alpha, Cr=None, C0r=None ,oil=oil_iso_vg_1500, 
                 material=material_iso, typ_contact=None):
        RadialBearing.__init__(self, d, D, i, Z, Dw, alpha, Cr, C0r, oil, material, typ_contact)
        self.coeff_baselife = 3.
        
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
            
class ConceptSphericalBallBearing(RadialBearing):
    def __init__(self, d, D, i, Z, Dw, alpha, Cr=None, C0r=None ,oil=oil_iso_vg_1500, 
                 material=material_iso, typ_contact=None):
        RadialBearing.__init__(self, d, D, i, Z, Dw, alpha, Cr, C0r, oil, material, typ_contact)
        self.coeff_baselife = 3.
        
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
    def __init__(self, d, D, i, Z, Dw, alpha, Cr=None, C0r=None ,oil=oil_iso_vg_1500, 
                 material=material_iso, typ_contact=None):
        RadialBearing.__init__(self, d, D, i, Z, Dw, alpha, Cr, C0r, oil, material, typ_contact)
        self.coeff_baselife = 10/3.
        
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
        coeff0 = fsolve(f,ec*Cu/Pr)[0]
        coeff = min(coeff0, ec*Cu/Pr)
        a_iso = 0.1*(f(coeff)**(-9.185))
        return a_iso

class RadialRollerBearing(ConceptRadialRollerBearing):
    #Roulement à rouleaux
    def __init__(self, typ, B, d, D, d1, D1, Lw, Dw, r_roller, E, F, Z, i, 
                 alpha, bm=1.1, oil=oil_iso_vg_1500, material=material_iso):
        ConceptRadialRollerBearing.__init__(self, i, Z, Dw, alpha ,oil, material, typ_contact = 'linear_contact')
        self.typ = typ
        self.B = B
        self.d1 = d1
        self.D1 = D1
        self.E = E
        self.F = F
        #diametre rouleau moyen
        self.Dwe = Dw
        self.Lw = Lw
        self.r_roller = r_roller
        self.bm = bm
        self.Dpw,self.Lwe,self.jeu,self.ep = self.DefParam()
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
        Lwe = self.Lw-2*self.r_roller
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
    
    def InternalRingContour(self):
        if self.typ == 'NU':
            p = [vm.Point2D((0,self.d/2.))]
            p.append(vm.Point2D((-self.B/2.,self.d/2.)))
            p.append(vm.Point2D((-self.B/2.,self.F/2.)))
            p.append(vm.Point2D((self.B/2.,self.F/2.)))
            p.append(vm.Point2D((self.B/2.,self.d/2.)))
            p.append(p[0])
            # TODO handle radius of rollers
#            ref=vm.Contour2D(primitives2D.RoundedLines2D(p,{1:self.r_roller,2:self.r_roller,3:self.r_roller,4:self.r_roller},False).primitives)
            ref = vm.Contour2D(primitives2D.RoundedLines2D(p,{},False).primitives)
        elif self.typ == 'N' or self.typ == 'NF':
            p = [vm.Point2D((0,self.d/2.))]
            p.append(vm.Point2D((-self.B/2.,self.d/2.)))
            p.append(vm.Point2D((-self.B/2.,self.d1/2.)))
            p.append(vm.Point2D((-self.B/2.+self.ep,self.d1/2.)))
            p.append(vm.Point2D((-self.B/2.+self.ep,self.F/2.)))
            p.append(vm.Point2D((self.B/2.-self.ep,self.F/2.)))
            p.append(vm.Point2D((self.B/2.-self.ep,self.d1/2.)))
            p.append(vm.Point2D((self.B/2.,self.d1/2.)))
            p.append(vm.Point2D((self.B/2.,self.d/2.)))
            p.append(p[0])
#            ref=vm.Contour2D(primitives2D.RoundedLines2D(p,{1:self.r_roller,2:self.r_roller,3:self.r_roller,4:self.r_roller,5:self.r_roller,6:self.r_roller,7:self.r_roller,8:self.r_roller},False).primitives)
            ref = vm.Contour2D(primitives2D.RoundedLines2D(p,{},False).primitives)
        elif self.typ == 'NJ':
            p = [vm.Point2D((0,self.d/2.))]
            p.append(vm.Point2D((-self.B/2.,self.d/2.)))
            p.append(vm.Point2D((-self.B/2.,self.d1/2.)))
            p.append(vm.Point2D((-self.B/2.+self.ep,self.d1/2.)))
            p.append(vm.Point2D((-self.B/2.+self.ep,self.F/2.)))
            p.append(vm.Point2D((self.B/2.,self.F/2.)))
            p.append(vm.Point2D((self.B/2.,self.d/2.)))
            p.append(p[0])
#            ref=vm.Contour2D(primitives2D.RoundedLines2D(p,{1:self.r_roller,2:self.r_roller,3:self.r_roller,4:self.r_roller,5:self.r_roller,6:self.r_roller},False).primitives)
            ref = vm.Contour2D(primitives2D.RoundedLines2D(p,{},False).primitives)
        return ref
    
    def ExternalRingContour(self):
        if self.typ == 'N':
            p = [vm.Point2D((0,self.E/2.))]
            p.append(vm.Point2D((-self.B/2.,self.E/2.)))
            p.append(vm.Point2D((-self.B/2.,self.D/2.)))
            p.append(vm.Point2D((self.B/2.,self.D/2.)))
            p.append(vm.Point2D((self.B/2.,self.E/2.)))
            p.append(p[0])
            ref = vm.Contour2D(primitives2D.RoundedLines2D(p,{1:self.r_roller,2:self.r_roller,
                            3:self.r_roller,4:self.r_roller},False).primitives)
        elif self.typ == 'NU' or self.typ == 'NJ':
            p = [vm.Point2D((0,self.E/2.))]
            p.append(vm.Point2D((-self.B/2.+self.ep,self.E/2)))
            p.append(vm.Point2D((-self.B/2.+self.ep,self.D1/2)))
            p.append(vm.Point2D((-self.B/2.,self.D1/2.)))
            p.append(vm.Point2D((-self.B/2.,self.D/2.)))
            p.append(vm.Point2D((self.B/2.,self.D/2.)))
            p.append(vm.Point2D((self.B/2.,self.D1/2.)))
            p.append(vm.Point2D((self.B/2.-self.ep,self.D1/2.)))
            p.append(vm.Point2D((self.B/2.-self.ep,self.E/2.)))
            p.append(p[0])
            ref = vm.Contour2D(primitives2D.RoundedLines2D(p,{1:self.r_roller,2:self.r_roller,
                            3:self.r_roller,4:self.r_roller,5:self.r_roller,6:self.r_roller,7:self.r_roller,8:self.r_roller},False).primitives)
        elif self.typ == 'NF':
            p = [vm.Point2D((0,self.E/2))]
            p.append(vm.Point2D((-self.B/2+self.ep,self.E/2)))
            p.append(vm.Point2D((-self.B/2+self.ep,self.D1/2)))
            p.append(vm.Point2D((-self.B/2,self.D1/2)))
            p.append(vm.Point2D((-self.B/2,self.D/2)))
            p.append(vm.Point2D((self.B/2,self.D/2)))
            p.append(vm.Point2D((self.B/2,self.E/2)))
            p.append(p[0])
            ref = vm.Contour2D(primitives2D.RoundedLines2D(p,{1:self.r_roller,2:self.r_roller,
                            3:self.r_roller,4:self.r_roller,5:self.r_roller,6:self.r_roller},False).primitives)
        return ref
    
    def RollerContour(self):
        p = [vm.Point2D((0,0.))]
        p.append(vm.Point2D((-self.Lw/2,0.)))
        p.append(vm.Point2D((-self.Lw/2,self.Dw/2)))
        p.append(vm.Point2D((self.Lw/2,self.Dw/2)))
        p.append(vm.Point2D((self.Lw/2,0.)))
        p.append(p[0])
        ref = vm.Contour2D(primitives2D.RoundedLines2D(p,{2:self.r_roller,3:self.r_roller},
                                                       False).primitives)
        return ref
        
    def PlotData(self, x, heights, ys, zs, labels = True):
        transversal_plot_data = []
        axial_plot_data = []
        
        component_height = 0.5 * (self.D-self.d)
                
        # TODO Check rollers seem to be to low
        y = ys[0]
        z = zs[0]
        yroller = 0.25*(self.E+self.F)
        # interface of upper section
        axial_plot_data.append({'type' : 'rect',
                            'x' : x - 0.5 * self.B,
                            'y' : heights[0] + 0.5*self.d ,
                            'width' : self.B,
                            'height' : component_height,
                            'color' : (0, 0, 0),
                            'size' : 1,
                            'dash' : 'none'})
            
        # Roller of upper section
        axial_plot_data.append({'type' : 'rect',
                            'x' : x - 0.5*self.Lw,
                            'y' : heights[0] +  yroller -0.5*self.Dw,
                            'width' : self.Lw,
                            'height' : self.Dw,
                            'color' : (0, 0, 0),
                            'size' : 1,
                            'dash' : 'none'})


        # interface of lower section
        axial_plot_data.append({'type' : 'rect',
                            'x' : x - 0.5 * self.B,
                            'y' : heights[0] - 0.5*self.D ,
                            'width' : self.B,
                            'height' : component_height,
                            'color' : (0, 0, 0),
                            'size' : 1,
                            'dash' : 'none'})
            
        # Roller of upper section
        axial_plot_data.append({'type' : 'rect',
                            'x' : x - 0.5*self.Lw,
                            'y' : heights[0] - yroller - 0.5*self.Dw,
                            'width' : self.Lw,
                            'height' : self.Dw,
                            'color' : (0, 0, 0),
                            'size' : 1,
                            'dash' : 'none'})
                    
            
                    
        
        transversal_plot_data.append({'type' : 'circle',
                                  'cx' : y,
                                  'cy' : z,
                                  'r' : 0.5 * self.d,
                                  'color' : [0, 0, 0],
                                  'size' : 1,
                                  'group' : 3,
                                  'dash' : 'none',})
            
        transversal_plot_data.append({'type' : 'circle',
                                  'cx' : y,
                                  'cy' : z,
                                  'r' : 0.5 * self.D,
                                  'color' : [0, 0, 0],
                                  'size' : 1,
                                  'group' : 3,
                                  'dash' : 'none',})
    
        for i in range(self.Z):
            theta=2*npy.pi/self.Z*i
            transversal_plot_data.append({'type' : 'circle',
                                      'cx' : y + yroller * math.cos(theta),
                                      'cy' : z + yroller * math.sin(theta),
                                      'r' : 0.5 * self.Dw,
                                      'color' : [0, 0, 0],
                                      'size' : 1,
                                      'group' : 3,
                                      'dash' : 'none',})
                    

        return axial_plot_data, transversal_plot_data
        
    def VolumeModel(self, center = (0,0,0), axis = (1,0,0)):
        center = vm.Point3D(npy.round(center,6))
        x = vm.Vector3D(axis)
        x.vector = x.vector/x.Norm()
        
        y = x.RandomUnitNormalVector()
        y.vector = npy.round(y.vector,3)
        y.vector = y.vector/y.Norm() 
        
        z=vm.Vector3D(npy.cross(x.vector,y.vector))
        
        #bague interne
        IRC=self.InternalRingContour()        
        irc=primitives3D.RevolvedProfile(center,x, z,[IRC], center, x,angle=2*math.pi,name='irc')
        #bague externe
        ERC=self.ExternalRingContour()
        erc=primitives3D.RevolvedProfile(center,x, z, [ERC], center, x, angle=2*math.pi,name='erc')
        #roller
        ROL=self.RollerContour()
        radius=self.F/2.+self.jeu+self.Dw/2.
        rol=[]
        theta=2*npy.pi/self.Z
        for zi in range(int(self.Z)):
            center_roller = center + radius*math.cos(zi*theta) * y + radius*math.sin(zi*theta) * z
            rol.append(primitives3D.RevolvedProfile(center_roller, x, z, [ROL],
                                                    center_roller, x,
                                                    angle=2*math.pi,name='rol'))
        
        tot=[irc,erc]+rol
        model=vm.VolumeModel(tot)
        return model

    def FreeCADExport(self,file_path,export_types=['fcstd']):
        model = self.VolumeModel()
        model.FreeCADExport('python',file_path,'/usr/lib/freecad/lib',export_types)

class DrawnCupNeedleRollerBearing(RadialRollerBearing):
    #Douille à aiguilles
    def __init__(self, typ, B, d, D, d1, D1, Lw, Dw, r_roller, E, F, Z, i,
                 alpha,bm=1, weibull_e=9/8., weibull_c=31/3., weibull_h=7/3.,
                 B1=551.13373/0.483, mu_delta=0.83, c_gamma=0.05,
                 oil_name='iso_vg_100'):
        RadialRollerBearing.__init__(self, typ, B, d, D, d1, D1, Lw, Dw, r_roller, E,
                                     F, Z, i, alpha, bm, weibull_e, weibull_c,
                                     weibull_h, B1, mu_delta, c_gamma, oil_name)
        
class NeedleRollerBearing(RadialRollerBearing):
    #Cage à aiguilles
    def __init__(self, typ, B, d, D, d1, D1, Lw, Dw, r_roller, E, F, Z, i,
                 alpha, bm=1, weibull_e=9/8., weibull_c=31/3., weibull_h=7/3.,
                 B1=551.13373/0.483, mu_delta=0.83, c_gamma=0.05,
                 oil_name='iso_vg_100'):
        RadialRollerBearing.__init__(self, typ, B, d, D, d1, D1, Lw, Dw, r_roller, E,
                                     F, Z, i, alpha, bm, weibull_e, weibull_c,
                                     weibull_h, B1, mu_delta, c_gamma, oil_name)

class SphericalRollerBearing(RadialRollerBearing):
    #Roulement à rotule à rouleaux
    def __init__(self, typ, B, d, D, d1, D1, Lw, Dw, r_roller, E, F, Z, i,
                 alpha,bm=1.15, weibull_e=9/8., weibull_c=31/3., weibull_h=7/3.,
                 B1=551.13373/0.483, mu_delta=0.83, c_gamma=0.05,
                 oil_name='iso_vg_100'):
        RadialRollerBearing.__init__(self, typ, B, d, D, d1, D1, Lw, Dw, r_roller, E,
                                     F, Z, i, alpha, bm, weibull_e, weibull_c,
                                     weibull_h, B1, mu_delta, c_gamma,
                                     oil_name)