import numpy as npy
#import math as mt
from scipy import interpolate
#import os
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math


from scipy.optimize import minimize


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
        self.oil_kinematic_viscosity_curve=self.KinematicViscosity(oil_data)
        self.dict_oil_contamination=dict_oil_contamination
    
    def Dict(self):
        return self.__dict__
    
    def FunCoeff(self,x,data,type_x='Linear',type_y='Linear'):
        if type_x=='Log': 
            x=npy.log10(x)
        f = interpolate.interp1d(list(data[:,0]),list(data[:,1]), fill_value='extrapolate')
        sol=float(f(x))
        if type_y=='Log':
            sol=10**sol
        return sol
    
    def KinematicViscosity(self,oil_kinematic_viscosity):
        oil_kinematic_viscosity_curve={}
        for (key,val) in oil_kinematic_viscosity.items():
            val_np=npy.array(val)
            if key not in ['x','y']:
                oil_kinematic_viscosity_curve[key]={}
                A=(npy.log10(npy.log10(0.6+val_np[0,1]))-npy.log10(npy.log10(0.6+val_np[-1,1])))/(npy.log10(val_np[0,0])-npy.log10(val_np[-1,0]))
                B=npy.log10(npy.log10(0.6+val_np[0,1]))-A*npy.log10(val_np[0,0])
                oil_kinematic_viscosity_curve[key]['A']=A
                oil_kinematic_viscosity_curve[key]['B']=B
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
    def __init__(self,weibull_e=9/8,weibull_c=31/3,weibull_h=7/3,B1=551.13373/0.483,mu_delta=0.83,c_gamma=0.05):
        self.weibull_e=weibull_e
        self.weibull_c=weibull_c
        self.weibull_h=weibull_h
        self.B1=B1
        self.mu_delta=mu_delta
        self.c_gamma=c_gamma

    def Dict(self):
        return self.__dict__

material_iso=Material()
    
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
#class RadialBallBearing(persistent.Persistent):
    
class RadialRollerBearing:
    #Roulement à rouleaux
    def __init__(self, typ, B, d, D, d1, D1, Lw, Dw, r_roller, E, F, Z, i, 
                 alpha,bm=1.1, oil=oil_iso_vg_1500, material=material_iso):
        self.typ=typ
        self.B=B
        self.d=d
        self.D=D
        self.d1=d1
        self.D1=D1
        self.E=E
        self.F=F
        self.Dw=Dw
        #diametre rouleau moyen
        self.Dwe=Dw
        self.Lw=Lw
        self.Z=Z
        self.i=i
        self.r_roller=r_roller
        self.alpha=alpha
        self.bm=bm
        self.oil=oil
        self.material=material
        self.Dpw,self.Lwe,self.jeu,self.ep=self.DefParam()
        self.mass=self.Mass()
        
    def Update(self,d1,D1,E,F,Z):
        self.d1=d1
        self.D1=D1
        self.E=E
        self.F=F
        self.Z=Z
        self.Dpw,self.Lwe,self.jeu,self.ep=self.DefParam()
        self.mass=self.Mass()
        
    def DefParam(self):
        Dpw=(self.E+self.F)/2
        Lwe=self.Lw-2*self.r_roller
        jeu=(self.E-self.F-2*self.Dw)/4
        ep=(self.B-self.Lw-2*jeu)/2
        return Dpw,Lwe,jeu,ep
        
    def Mass(self):
        rho=7800
        m=self.Z*npy.pi*(self.Dw)**2/4*self.Lw*rho
        m+=(npy.pi*(self.D)**2/4-npy.pi*(self.E)**2/4)*self.B*rho
        m+=(npy.pi*(self.F)**2/4-npy.pi*(self.d)**2/4)*self.B*rho
        return m
    
    def BaseStaticLoad(self):
        #Charge radiale statique de base
        #besoin de convertir les dimensions en mm pour les formules ISO
        C0r=44*(1-(self.Dw*1e3*npy.cos(self.alpha))/(self.Dpw*1e3))*self.i*self.Z*self.Lwe*1e3*self.Dw*1e3*npy.cos(self.alpha)
        return C0r
    
    def EquivalentStaticLoad(self,Fr,Fa=None):
        #Charge radiale statique équivalente
        if self.alpha!=0:
            x0=0.5*self.i
            y0=0.22*1/npy.tan(self.alpha)*self.i
        else:
            x0,y0=1,0
        P0r=max(Fr,x0*Fr+y0*Fa)
        return P0r
        
    def BaseDynamicLoad(self):
        #Charge radiale dynamique de base
        mu=float((self.Dwe*1e3)*npy.cos(self.alpha)/(self.Dpw*1e3))
        fc=0.377*self.material.mu_delta*1/((2**((self.material.weibull_c+self.material.weibull_h-1)/(self.material.weibull_c-self.material.weibull_h+1)))*(0.5**(2*self.material.weibull_e/(self.material.weibull_c-self.material.weibull_h+1))))*self.material.B1*((1-mu)**((self.material.weibull_c+self.material.weibull_h-3)/(self.material.weibull_c-self.material.weibull_h+1))/((1+mu)**(2*self.material.weibull_e/(self.material.weibull_c-self.material.weibull_h+1))))*(mu**(2/(self.material.weibull_c-self.material.weibull_h+1)))*(1+(1.04*((1-mu)/(1+mu))**((self.material.weibull_c+self.material.weibull_h+2*self.material.weibull_e-3)/(self.material.weibull_c-self.material.weibull_h+1)))**((self.material.weibull_c-self.material.weibull_h+1)/2))**(-2/(self.material.weibull_c-self.material.weibull_h+1))
        Cr=fc*self.bm*self.i*((self.Lwe*1e3)*npy.cos(self.alpha))**((self.material.weibull_c-self.material.weibull_h-1)/(self.material.weibull_c-self.material.weibull_h+1))*self.Z**((self.material.weibull_c-self.material.weibull_h-2*self.material.weibull_e+1)/(self.material.weibull_c-self.material.weibull_h+1))*(self.Dwe*1e3)**((self.material.weibull_c-self.material.weibull_h-3)/(self.material.weibull_c-self.material.weibull_h+1))
        return Cr
    
    def EquivalentDynamicLoad(self,Fr,Fa=0):
        #Charge radiale dynamique équivalente
        e=1.5*npy.tan(self.alpha)
        if self.alpha!=0:
            if self.i==1:
                if Fa/Fr<=e:
                    Pr=Fr
                else:
                    Pr=0.4*Fr+0.4/(npy.tan(self.alpha))*Fa
            elif self.i==2:
                if Fa/Fr<=e:
                    Pr=Fr+0.45/(npy.tan(self.alpha))*Fa
                else:
                    Pr=0.67*Fr+0.67/(npy.tan(self.alpha))*Fa
        else:
            Pr=Fr
        return Pr
    
    def BaseLifeTime(self,Fr,Fa=0):
        # Durée de vie en millions de tour associée à une fiabilité de 90%
        Cr=self.BaseDynamicLoad()
        Pr=self.EquivalentDynamicLoad(Fr,Fa)
        L10=(Cr/Pr)**(10/3)
        return L10
        
    def AdjustedLifeTime(self,Fr,n,Fa=0,S=0.9,T=70):
        # Durée de vie corrigée en millions de tour associée à une fiabilité de S% pour un roulement tournant à la vitesse n (rad/s) et à la température de l'huile T
        a1=(1-self.material.c_gamma)*(npy.log(1/S)/npy.log(100/90))**(1/self.material.weibull_e)+self.material.c_gamma
        L10=self.BaseLifeTime(Fr,Fa)
        Pr=self.EquivalentDynamicLoad(Fr,Fa)
        C0r=self.BaseStaticLoad()
        # viscosité cinématique de référence
        if n<(1000*2*npy.pi/60):
            nu1=45000*(n*60/(2*npy.pi))**(-0.83)*(self.Dpw*1e3)**(-0.5)
        else:
            nu1=4500*(n*60/(2*npy.pi))**(-0.5)*(self.Dpw*1e3)**(-0.5)
        
        coeff_oil=self.oil.oil_kinematic_viscosity_curve
        nu=10**(10**(coeff_oil['A']*npy.log10(T)+coeff_oil['B']))-0.6
        kappa=nu/nu1
        #définition du paramètre de contamination
        ec=self.oil.OilParameterContamination(self.Dpw,3)
        #calcul de la limite de charge en fatigue
        if self.Dpw<0.1:
            Cu=C0r/8.2
        else:
            Cu=C0r/8.2*(100/(self.Dpw*1e3))**0.3
        #calcul du coefficient a_iso
        # TODO: Il y a des fois des warnings, il y a peut etre encore un problème!
        if kappa<0.4:
            a_iso=0.1*((1-(1.5859-1.3993/(kappa**0.054381))*((ec*Cu/Pr)**0.4))**(-9.185))
        elif kappa<1:
            a_iso=0.1*((1-(1.5859-1.2348/(kappa**0.19087))*((ec*Cu/Pr)**0.4))**(-9.185))
        else:
            kappa=min(kappa,4)
            a_iso=0.1*((1-(1.5859-1.2348/(kappa**0.071739))*((ec*Cu/Pr)**0.4))**(-9.185))
#            print(kappa,ec,Cu,Pr, a_iso)
        a_iso=min(50,a_iso)
        #calcul de la durée de vie corrigée
        Lnm=a1*a_iso*L10
        return Lnm
        
    def Dict(self):

        d={}
        for k,v in self.__dict__.items():
            tv=type(v)
            if tv==npy.int64:
                d[k]=int(v)
            elif tv==npy.float64:
                d[k]=float(v)
            else:
                d[k]=v

        d['oil'] = self.oil.Dict()
        d['material'] = self.material.Dict()
        return d
    
    def InternalRingContour(self):
        if self.typ=='NU':
            p=[vm.Point2D((0,self.d/2))]
            p.append(vm.Point2D((-self.B/2,self.d/2)))
            p.append(vm.Point2D((-self.B/2,self.F/2)))
            p.append(vm.Point2D((self.B/2,self.F/2)))
            p.append(vm.Point2D((self.B/2,self.d/2)))
            p.append(p[0])
            ref=vm.Contour2D(primitives2D.RoundedLines2D(p,{1:self.r_roller,2:self.r_roller,3:self.r_roller,4:self.r_roller},False).primitives)
        elif self.typ=='N' or self.typ=='NF':
            p=[vm.Point2D((0,self.d/2))]
            p.append(vm.Point2D((-self.B/2,self.d/2)))
            p.append(vm.Point2D((-self.B/2,self.d1/2)))
            p.append(vm.Point2D((-self.B/2+self.ep,self.d1/2)))
            p.append(vm.Point2D((-self.B/2+self.ep,self.F/2)))
            p.append(vm.Point2D((self.B/2-self.ep,self.F/2)))
            p.append(vm.Point2D((self.B/2-self.ep,self.d1/2)))
            p.append(vm.Point2D((self.B/2,self.d1/2)))
            p.append(vm.Point2D((self.B/2,self.d/2)))
            p.append(p[0])
            ref=vm.Contour2D(primitives2D.RoundedLines2D(p,{1:self.r_roller,2:self.r_roller,3:self.r_roller,4:self.r_roller,5:self.r_roller,6:self.r_roller,7:self.r_roller,8:self.r_roller},False).primitives)
        elif self.typ=='NJ':
            p=[vm.Point2D((0,self.d/2))]
            p.append(vm.Point2D((-self.B/2,self.d/2)))
            p.append(vm.Point2D((-self.B/2,self.d1/2)))
            p.append(vm.Point2D((-self.B/2+self.ep,self.d1/2)))
            p.append(vm.Point2D((-self.B/2+self.ep,self.F/2)))
            p.append(vm.Point2D((self.B/2,self.F/2)))
            p.append(vm.Point2D((self.B/2,self.d/2)))
            p.append(p[0])
            ref=vm.Contour2D(primitives2D.RoundedLines2D(p,{1:self.r_roller,2:self.r_roller,3:self.r_roller,4:self.r_roller,5:self.r_roller,6:self.r_roller},False).primitives)
        return ref
    
    def ExternalRingContour(self):
        if self.typ=='N':
            p=[vm.Point2D((0,self.E/2))]
            p.append(vm.Point2D((-self.B/2,self.E/2)))
            p.append(vm.Point2D((-self.B/2,self.D/2)))
            p.append(vm.Point2D((self.B/2,self.D/2)))
            p.append(vm.Point2D((self.B/2,self.E/2)))
            p.append(p[0])
            ref=vm.Contour2D(primitives2D.RoundedLines2D(p,{1:self.r_roller,2:self.r_roller,3:self.r_roller,4:self.r_roller},False).primitives)
        elif self.typ=='NU' or self.typ=='NJ':
            p=[vm.Point2D((0,self.E/2))]
            p.append(vm.Point2D((-self.B/2+self.ep,self.E/2)))
            p.append(vm.Point2D((-self.B/2+self.ep,self.D1/2)))
            p.append(vm.Point2D((-self.B/2,self.D1/2)))
            p.append(vm.Point2D((-self.B/2,self.D/2)))
            p.append(vm.Point2D((self.B/2,self.D/2)))
            p.append(vm.Point2D((self.B/2,self.D1/2)))
            p.append(vm.Point2D((self.B/2-self.ep,self.D1/2)))
            p.append(vm.Point2D((self.B/2-self.ep,self.E/2)))
            p.append(p[0])
            ref=vm.Contour2D(primitives2D.RoundedLines2D(p,{1:self.r_roller,2:self.r_roller,3:self.r_roller,4:self.r_roller,5:self.r_roller,6:self.r_roller,7:self.r_roller,8:self.r_roller},False).primitives)
        elif self.typ=='NF':
            p=[vm.Point2D((0,self.E/2))]
            p.append(vm.Point2D((-self.B/2+self.ep,self.E/2)))
            p.append(vm.Point2D((-self.B/2+self.ep,self.D1/2)))
            p.append(vm.Point2D((-self.B/2,self.D1/2)))
            p.append(vm.Point2D((-self.B/2,self.D/2)))
            p.append(vm.Point2D((self.B/2,self.D/2)))
            p.append(vm.Point2D((self.B/2,self.E/2)))
            p.append(p[0])
            ref=vm.Contour2D(primitives2D.RoundedLines2D(p,{1:self.r_roller,2:self.r_roller,3:self.r_roller,4:self.r_roller,5:self.r_roller,6:self.r_roller},False).primitives)
        return ref
    
    def RollerContour(self):
        p=[vm.Point2D((0,0.))]
        p.append(vm.Point2D((-self.Lw/2,0.)))
        p.append(vm.Point2D((-self.Lw/2,self.Dw/2)))
        p.append(vm.Point2D((self.Lw/2,self.Dw/2)))
        p.append(vm.Point2D((self.Lw/2,0.)))
        p.append(p[0])
        ref=vm.Contour2D(primitives2D.RoundedLines2D(p,{2:self.r_roller,3:self.r_roller},False).primitives)
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
        radius=self.F/2+self.jeu+self.Dw/2
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
                 alpha,bm=1, weibull_e=9/8, weibull_c=31/3, weibull_h=7/3,
                 B1=551.13373/0.483, mu_delta=0.83, c_gamma=0.05,
                 oil_name='iso_vg_100'):
        RadialRollerBearing.__init__(typ, B, d, D, d1, D1, Lw, Dw, r_roller, E,
                                     F, Z, i, alpha, bm, weibull_e, weibull_c,
                                     weibull_h, B1, mu_delta, c_gamma, oil_name)
        
class NeedleRollerBearing(RadialRollerBearing):
    #Cage à aiguilles
    def __init__(self, typ, B, d, D, d1, D1, Lw, Dw, r_roller, E, F, Z, i,
                 alpha, bm=1, weibull_e=9/8, weibull_c=31/3, weibull_h=7/3,
                 B1=551.13373/0.483, mu_delta=0.83, c_gamma=0.05,
                 oil_name='iso_vg_100'):
        RadialRollerBearing.__init__(typ, B, d, D, d1, D1, Lw, Dw, r_roller, E,
                                     F, Z, i, alpha, bm, weibull_e, weibull_c,
                                     weibull_h, B1, mu_delta, c_gamma, oil_name)

class SphericalRollerBearing(RadialRollerBearing):
    #Roulement à rotule à rouleaux
    def __init__(self, typ, B, d, D, d1, D1, Lw, Dw, r_roller, E, F, Z, i,
                 alpha,bm=1.15, weibull_e=9/8, weibull_c=31/3, weibull_h=7/3,
                 B1=551.13373/0.483, mu_delta=0.83, c_gamma=0.05,
                 oil_name='iso_vg_100'):
        RadialRollerBearing.__init__(typ, B, d, D, d1, D1, Lw, Dw, r_roller, E,
                                     F, Z, i, alpha, bm, weibull_e, weibull_c,
                                     weibull_h, B1, mu_delta, c_gamma,
                                     oil_name)
# TODO: Cumulative damages
class BearingCombination:
    """
    Objet avec 3 fonctions de selection des roulements cylindriques
     * Combinatoire sur les dimensions externe ISO
     * Combinatoire en prenant en compte les règles SKF
     * Estimation des durées de vie et charge dynamique et fonction de tri

    """
    def __init__(self):
        
        self.solutions=[]
       
    def OptimizerBearing(self, d, D, B, Fr, Fa, n, L10=None, C0r=None, Cr=None,
                         Lnm=None, grade=['Gr_gn'], S=0.9, T=40,
                         oil=oil_iso_vg_1500, material=material_iso,
                         nb_sol=1, maxi=None, mini=None, rsmin=None, typ='NF',
                         verbose = False):

        
        err_default=0.05
        def def_inter(data):
            if data==None:
                sol=[-npy.inf,npy.inf]
            else:
                if 'nom' in data.keys():
                    if 'err' in data.keys():
                        err=data['err']
                    else:
                        err=err_default
                    sol=[data['nom']*(1-err),data['nom']*(1+err)]
                elif 'min' in data.keys():
                    sol=[data['min'],data['max']]
            return sol
        
        
        #Approche combinatoire sur les dimensions externe d/D/B
        liste_keys=[]
        interval_d=def_inter(d)
        interval_D=def_inter(D)
        interval_B=def_inter(B)
        lk_d=[]
        for item_d in dico_rlts_iso.keys():
            if (item_d<=interval_d[1]) and (item_d>=interval_d[0]):
                lk_d.append([item_d])
        liste_keys=lk_d
        lk_D=[]
        for [item_d] in liste_keys:
            for item_D in dico_rlts_iso[item_d].keys():
                if (item_D<=interval_D[1]) and (item_D>=interval_D[0]):
                    lk_D.append([item_d,item_D])
        liste_keys=lk_D
        lk_B=[]
        for item_d,item_D in liste_keys:
            for item_B in dico_rlts_iso[item_d][item_D].keys():
                if (item_B<=interval_B[1]) and (item_B>=interval_B[0]):
                    lk_B.append([item_d,item_D,item_B])
        liste_keys=lk_B
        liste_sol_iso=[]
        for item_d,item_D,item_B in liste_keys:
            dk=[item_d,item_D,item_B]
            dk.append(list(dico_rlts_iso[item_d][item_D][item_B].keys())[0])
            dk.append(list(dico_rlts_iso[item_d][item_D][item_B].values())[0])
            liste_sol_iso.append(dk) #Liste dk contient [d,D,B,rsmin,serial]
        
        #Approche combinatoire sur les rouleaux
        def analyse_rule(Y,inter_min,inter_max,input_var,data_var):
            for (var_x,var_y,typ),[a,b] in dico_rules.items():
                if var_y==Y:
                    if var_x in input_var:
                        ind_x=input_var.index(var_x)
                        Var_x=data_var[ind_x]
                        if typ=='inf':
                            inter_min=max(inter_min,a*Var_x+b)
                        else:
                            inter_max=min(inter_max,a*Var_x+b)
            return inter_min,inter_max
                            
        liste_sol_roller_iso=[]
        for [d,D,B,rsmin,serial] in liste_sol_iso:
            
            #Estimation des plages admissibles de E et F
            E_inf,E_sup=-npy.inf,npy.inf
            F_inf,F_sup=-npy.inf,npy.inf
            E_inf,E_sup=analyse_rule('E',E_inf,E_sup,['d','D','B'],[d,D,B])
            F_inf,F_sup=analyse_rule('F',F_inf,F_sup,['d','D','B'],[d,D,B])
                            
            #Analyse borne sup/inf de Lw et Dw
            Lw_inf,Lw_sup=-npy.inf,npy.inf
            Dw_inf,Dw_sup=E_inf-F_sup,E_sup-F_inf
            Lw_inf,Lw_sup=analyse_rule('Lw',Lw_inf,Lw_sup,['d','D','B'],[d,D,B])
            Dw_inf,Dw_sup=analyse_rule('Dw',Dw_inf,Dw_sup,['d','D','B'],[d,D,B])
            
            #Choix des dimensions rouleaux compatible
            Dw_sup=min(Dw_sup,(D-d)/2)
            for Dw in dico_roller_iso.keys():
                if (Dw>=Dw_inf) and (Dw<=Dw_sup):
                    for Lw in dico_roller_iso[Dw].keys():
                        if (Lw>=Lw_inf) and (Lw<=Lw_sup):
                            E_inf,E_sup=analyse_rule('E',-npy.inf,npy.inf,['d','D','B','Dw','Lw'],[d,D,B,Dw,Lw])
                            F_inf,F_sup=analyse_rule('F',-npy.inf,npy.inf,['d','D','B','Dw','Lw'],[d,D,B,Dw,Lw])
                            if ((E_inf-F_sup)/2<=Dw) and ((E_sup-F_inf)/2>=Dw):
                                liste_sol_roller_iso.append([d,D,B,rsmin,serial,Dw,Lw,E_inf,E_sup,F_inf,F_sup])


        #Analyse de détail des roulements
        a_ED1_inf,b_ED1_inf=dico_rules[('E','D1','inf')]
        a_ED1_sup,b_ED1_sup=dico_rules[('E','D1','sup')]
        a_Fd1_inf,b_Fd1_inf=dico_rules[('F','d1','inf')]
        a_Fd1_sup,b_Fd1_sup=dico_rules[('F','d1','sup')]
        estim_masse=[]
        #Construction d'une fonctionnelle pour le tri afin d'engager une optimisation continue
        for [d,D,B,rsmin,serial,Dw,Lw,E_inf,E_sup,F_inf,F_sup] in liste_sol_roller_iso:
            E=(E_inf+E_sup)/2
            F=E-2*Dw
            Zmax=int(2*npy.pi/(2*npy.arcsin((Dw/2)/(F/2+Dw/2))))
            masse_elem=(npy.pi*D**2/4-npy.pi*E**2/4)*B
            masse_elem+=(npy.pi*F**2/4-npy.pi*d**2/4)*B
            masse_elem+=npy.pi*Dw**2/4*Lw*Zmax
            estim_masse.append(masse_elem)
        
        #Optimisation Minimize du détail du roulement (E,D1 et d1)
        fonct_sort=npy.argsort(estim_masse)
        liminf_lifetime=False
        # BUG: rien n'empeche i de dépasser du tableau!
        # Vérifier ma correction
        
        # TODO: permettre de ne donner qu'un borne min au L10/ Lnm
        
        i=0
        iter_max = len(liste_sol_roller_iso)
        while (not liminf_lifetime) and (i<iter_max):
            [d,D,B,rsmin,serial,Dw,Lw,E_inf,E_sup,F_inf,F_sup]=liste_sol_roller_iso[fonct_sort[i]]
            i+=1
            # Var X: (E,D1,d1)
            R1=RadialRollerBearing('N',B,d,D,0,0,Lw,Dw,rsmin,0,0,0,1,0)
            def fun(x):
                obj=0
                F=x[0]-2*Dw
                Zmax=int(2*npy.pi/(2*npy.arcsin((Dw/2)/(F/2+Dw/2))))
                R1.Update(x[2],x[1],x[0],F,Zmax)
                l10=R1.BaseLifeTime(Fr,Fa)
                obj+=(1/l10)**2
                obj+=(x[0]-x[1])**2 #minimisation de la hauteur de l'épaulement externe
                obj+=(F-x[2])**2 #minimisation de la hauteur de l'épaulement interne
                return obj
            def ineg(x):
                ine=[]
                D1_inf=a_ED1_inf*x[0]+b_ED1_inf
                D1_sup=a_ED1_sup*x[0]+b_ED1_sup
                ine.append(x[1]-D1_inf)
                ine.append(D1_sup-x[1])
                F=x[0]-2*Dw
                d1_inf=a_Fd1_inf*F+b_Fd1_inf
                d1_sup=a_Fd1_sup*F+b_Fd1_sup
                ine.append(x[2]-d1_inf)
                ine.append(d1_sup-x[2])
                ine.append(x[0]-x[1])
                ine.append(x[1]-x[2])
                ine.append((x[0]+F)/2-x[2])  #d1 inférieur au diamètre passant par l'axe des rouleaux
                return ine
            cons = ({'type': 'ineq','fun' : ineg})
            valid_optim=1
            while valid_optim==1:
                Bound=[[E_inf,E_sup],[d,D],[d,D]]
                x0=(npy.array(Bound)[:,1]-npy.array(Bound)[:,0])*npy.random.random(3)+npy.array(Bound)[:,0]
                res = minimize(fun,x0, method='SLSQP', bounds=Bound,constraints=cons)
                if (min(ineg(res.x))>=0):
                    valid_optim=0
                    
            #Validation finale vis à vis des brones du CDC
            if valid_optim==0:
                x_opt=res.x
                F=x_opt[0]-2*Dw
                Zmax=int(2*npy.pi/(2*npy.arcsin((Dw/2)/(F/2+Dw/2))))
                R1.Update(x_opt[2],x_opt[1],x_opt[0],F,Zmax)
                l10=R1.BaseLifeTime(Fr,Fa)
                
                liminf_lifetime=True
                if not L10 is None:
                    if (l10<L10['min']) or (l10>L10['max']):
                        liminf_lifetime=False
                if not Lnm is None:
                    lnm=R1.AdjustedLifeTime(Fr,n,Fa,S,T)
                    if (lnm<Lnm['min']) or (lnm>Lnm['max']):
                        liminf_lifetime=False
                if not C0r is None:
                    c0r=R1.BaseStaticLoad()
                    if (c0r<C0r['min']) or (c0r>C0r['max']):
                        liminf_lifetime=False
                if not Cr is None:
                    cr=R1.BaseDynamicLoad()
                    if (cr<Cr['min']) or (cr>Cr['max']):
                        liminf_lifetime=False
                if liminf_lifetime==True:
                    self.solutions.append(R1)
                    if len(self.solutions)<nb_sol:
                        liminf_lifetime=False
                    if verbose:
                        print('Bearing solution n°{} with L10:{},D:{},d:{},B:{},Dw:{},Z:{}'.format(len(self.solutions),l10,D,d,B,Dw,Zmax))