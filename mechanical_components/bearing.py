import numpy as npy
import math as mt
from scipy import interpolate
import os
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math
from scipy.linalg import norm
from scipy.optimize import minimize,fsolve
from scipy.interpolate import splprep, splev
import cma
#from sympy import *
import itertools
from jinja2 import Environment, PackageLoader, select_autoescape

import mechanical_components.LibSvgD3 as LibSvg
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D

import persistent
import pandas
from pandas.plotting import scatter_matrix
import matplotlib.pyplot as plt
#from dessia_common import ResultsDBClient
#import pyDOE
from operator import itemgetter

# =============================================================
# Object matériau (data huile) nécessaire pour le calcul de la durée de vie corrigée
# =============================================================

class OilData(persistent.Persistent):
    
    def __init__(self):
        
        data_array,dico_axis,dico_nom,type_axis=self.KinematicViscosityData()
        self.oil_kinematic_viscosity=self.GenereCoeff(data_array,dico_axis,dico_nom,type_axis)
        for (key,val) in self.oil_kinematic_viscosity.items():
            if key not in ['x','y']:
                if type_axis['x']=='Log':
                    self.oil_kinematic_viscosity[key][:,0]=10**val[:,0]
                if type_axis['y']=='Log':
                    self.oil_kinematic_viscosity[key][:,1]=10**val[:,1]
        self.KinematicViscosity()
        
    def GenereCoeff(self,data_array,dico_axis,dico_nom,type_axis):
        structure={}
        vect_x=npy.linspace(data_array[dico_nom['x']][0,0],data_array[dico_nom['x']][-1,0],len(data_array[dico_nom['x']][:,0]))
        vect_y=npy.linspace(data_array[dico_nom['y']][0,1],data_array[dico_nom['y']][-1,1],len(data_array[dico_nom['y']][:,1]))
        
        for i,j in type_axis.items():
            if j=='Log':
                axe_reel=[mt.log10(dico_axis[i][0]),mt.log10(dico_axis[i][-1])]
            if j=='Linear':
                axe_reel=dico_axis[i]
            if i=='x':
                ax,bx=self.AxisLinear(axe_reel,vect_x)
            if i=='y':
                ay,by=self.AxisLinear(axe_reel,vect_y)

        for key,data in dico_nom.items():
            if not key in ['x','y']:
                export=[]
                for item in data_array[data]:
                    data_x=item[0]*ax+bx
                    data_y=item[1]*ay+by
                    export.append([data_x,data_y])
                export=npy.array(export)
                structure[key]=export
        structure['x']=type_axis['x']
        structure['y']=type_axis['y']
        return structure
    
    def KinematicViscosityData(self):
        data_svg=['m 70.029703,610.96286 48.808577,0 48.80858,1.06106 46.68647,0 48.80858,-1.06106 48.80859,1.06106 47.74752,0 48.80858,0 47.74753,0 47.74752,-1.06106 45.62541,0','m 70.029703,512.28465 -2.122112,-168.70792 1.061056,-71.09076 0,-245.103963','M 200.5396,35.870456 340.59901,201.39521 505.06271,332.96617','M 174.0132,41.175736 273.75248,168.50247 456.25413,327.66088','m 160.21947,53.908409 88.06766,118.838281 154.91419,147.4868','M 147.4868,84.679036 257.83663,220.49422 384.10231,335.08828','M 128.38779,89.984317 249.34818,240.65428 362.88119,345.69884','M 117.77723,121.816 246.16502,273.54702 372.43069,379.65263','m 108.22772,151.52557 110.34984,129.44885 174.0132,142.18151','M 132.63201,227.92161 265.26403,361.61468 411.68977,468.78135','M 117.77723,256.57013 247.22607,382.8358 405.32343,495.30775','M 113.533,297.95131 264.20297,429.52227 450.94884,545.17739','m 106.10561,330.84405 147.4868,125.20462 114.59406,73.21287','m 120.9604,387.08003 128.38778,101.86138 160.21948,97.61716','m 89.128713,401.93481 96.556107,81.70132 158.09736,103.9835','m 94.433993,454.98762 102.922447,79.57921 93.37293,58.35808']
        dico_nom={'x':0,'y':1,'iso_vg_1500':2,'iso_vg_1000':3,'iso_vg_680':4,'iso_vg_460':5,'iso_vg_320':6,'iso_vg_220':7,'iso_vg_150':8,'iso_vg_100':9,'iso_vg_68':10,'iso_vg_46':11,'iso_vg_32':12,'iso_vg_22':13,'iso_vg_15':14,'iso_vg_10':15}
        dico_axis={'x':[20,30,40,50,60,70,80,90,100,110,120],'y':[10,1000]}
        type_axis={'x':'Linear','y':'Log'}
        data_array=self.SVG2Array(data_svg)
        return data_array,dico_axis,dico_nom,type_axis
    
    def KinematicViscosity(self):
        self.oil_kinematic_viscosity_curve={}
        for (key,val) in self.oil_kinematic_viscosity.items():
            if key not in ['x','y']:
                self.oil_kinematic_viscosity_curve[key]={}
                A=(npy.log10(npy.log10(0.6+val[0,1]))-npy.log10(npy.log10(0.6+val[-1,1])))/(npy.log10(val[0,0])-npy.log10(val[-1,0]))
                B=npy.log10(npy.log10(0.6+val[0,1]))-A*npy.log10(val[0,0])
                self.oil_kinematic_viscosity_curve[key]['A']=A
                self.oil_kinematic_viscosity_curve[key]['B']=B
                
    def OilContamination(self,Dpw,grade=3):
        #grade correspond au niveau de propreté de lhuile, 1 très propre et 7 pour une contamination extreme
        oil_contamination={1:'ultra clean',2:'high level clean',3:'normal clean',4:'slight contamination',5:'classical contamination',6:'serious contamination',7:'major contamination'}
        if Dpw<0.1:
            if grade==1:
                self.ec=1
            elif grade==2:
                self.ec=0.7
            elif grade==3:
                self.ec=0.55
            elif grade==4:
                self.ec=0.4
            elif grade==5:
                self.ec=0.2
            elif grade==6:
                self.ec=0.05
            elif grade==7:
                self.ec=0
        else:
            if grade==1:
                self.ec=1
            elif grade==2:
                self.ec=0.85
            elif grade==3:
                self.ec=0.7
            elif grade==4:
                self.ec=0.5
            elif grade==5:
                self.ec=0.3
            elif grade==6:
                self.ec=0.05
            elif grade==7:
                self.ec=0
                
    def FunCoeff(self,x,data,type_x='Linear',type_y='Linear'):
        if type_x=='Log': 
            x=npy.log10(x)
        f = interpolate.interp1d(list(data[:,0]),list(data[:,1]), fill_value='extrapolate')
        sol=float(f(x))
        if type_y=='Log':
            sol=10**sol
        return sol
    
    def AxisLinear(self,axe,vect):
        a=(axe[0]-axe[-1])/(vect[0]-vect[-1])
        b=axe[0]-a*vect[0]
        return a,b
    
    def SVG2Array(self,data):
        # en entrée une liste de data SVG positionnee en relatif et 
        # en sortie liste Array avec positionnement en absolu
        export={}
        for i,dat in enumerate(data):
            data_temp=[]
            sol=dat.split(' ')
            for item in sol[1:]:
                temp=item.split(',')
                data_temp.append([float(temp[0]),float(temp[1])])
            export_temp=[data_temp[0]]
            if sol[0]=='m':
                for item in data_temp[1::]:
                    export_temp.append([item[0]+export_temp[-1][0],item[1]+export_temp[-1][1]])
            if sol[0]=='M':
                export_temp=data_temp
            export[i]=npy.array(export_temp)
        return export
    
# =============================================================
# Object générique roulement cylindrique
# =============================================================

class RadialRollerBearing(persistent.Persistent):
    def __init__(self,typ,B,d,D,d1,D1,Lw,Dw,r_roller,E,F,Z,i,alpha,O1,weibull_e=9/8,weibull_c=31/3,weibull_h=7/3,B1=551.13373/0.483,mu_delta=0.83,bm=1.1,c_gamma=0.05,oil_name='iso_vg_100'):
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
        self.Dpw=(self.E+self.F)/2
        self.Lwe=self.Lw-2*self.r_roller
        self.weibull_e=weibull_e
        self.weibull_c=weibull_c
        self.weibull_h=weibull_h
        self.B1=B1
        self.mu_delta=mu_delta
        self.bm=bm
        self.c_gamma=c_gamma
        self.oil_name=oil_name
        self.O1=O1
        self.jeu=(self.E-self.F-2*self.Dw)/4
        self.ep=(self.B-self.Lw-2*self.jeu)/2
    def BaseStaticLoad(self):
        #le système d'unité en entrée est le SI
        self.C0r=44*(1-(self.Dw*1e3*npy.cos(self.alpha))/(self.Dpw*1e3))*self.i*self.Z*self.Lwe*1e3*self.Dw*1e3*npy.cos(self.alpha)
    def EquivalentStaticLoad(self,Fr,Fa=None):
        x0=0.5*self.i
        y0=0.22*1/npy.tan(self.alpha)*self.i
        self.P0r=max(Fr,x0*Fr+y0*Fa)
    def BaseDynamicLoad(self):
        mu=float((self.Dwe*1e3)*npy.cos(self.alpha)/(self.Dpw*1e3))
        delta=self.mu_delta/mu
        self.fc=0.377*self.mu_delta*1/((2**((self.weibull_c+self.weibull_h-1)/(self.weibull_c-self.weibull_h+1)))*(0.5**(2*self.weibull_e/(self.weibull_c-self.weibull_h+1))))*self.B1*((1-mu)**((self.weibull_c+self.weibull_h-3)/(self.weibull_c-self.weibull_h+1))/((1+mu)**(2*self.weibull_e/(self.weibull_c-self.weibull_h+1))))*(mu**(2/(self.weibull_c-self.weibull_h+1)))*(1+(1.04*((1-mu)/(1+mu))**((self.weibull_c+self.weibull_h+2*self.weibull_e-3)/(self.weibull_c-self.weibull_h+1)))**((self.weibull_c-self.weibull_h+1)/2))**(-2/(self.weibull_c-self.weibull_h+1))
        self.Cr=self.fc*self.bm*self.i*((self.Lwe*1e3)*npy.cos(self.alpha))**((self.weibull_c-self.weibull_h-1)/(self.weibull_c-self.weibull_h+1))*self.Z**((self.weibull_c-self.weibull_h-2*self.weibull_e+1)/(self.weibull_c-self.weibull_h+1))*(self.Dwe*1e3)**((self.weibull_c-self.weibull_h-3)/(self.weibull_c-self.weibull_h+1))
    def EquivalentDynamicLoad(self,Fr,Fa=0):
        e=1.5*npy.tan(self.alpha)
        if self.i==1:
            if Fa/Fr<=e:
                self.Pr=Fr
            else:
                self.Pr=0.4*Fr+0.4/(npy.tan(self.alpha))*Fa
        elif self.i==2:
            if Fa/Fr<=e:
                self.Pr=Fr+0.45/(npy.tan(self.alpha))*Fa
            else:
                self.Pr=0.67*Fr+0.67/(npy.tan(self.alpha))*Fa
    def BaseLifeTime(self,Fr,Fa=0):
        # Durée de vie en millions de tour associée à une fiabilité de 90%
        self.BaseDynamicLoad()
        self.EquivalentDynamicLoad(Fr,Fa)
        self.L10=(self.Cr/self.Pr)**(10/3)
    def AdjustedLifeTime(self,Fr,n,Fa=0,S=0.9,T=40):
        # Durée de vie corrigée en millions de tour associée à une fiabilité de S% pour un roulement tournant à la vitesse n (rad/s) et à la température de l'huile T
        self.a1=(1-self.c_gamma)*(npy.log(1/S)/npy.log(100/90))**(1/self.weibull_e)+self.c_gamma
        # viscosité cinématique de référence
        if n<(1000*2*npy.pi/60):
            nu1=45000*(n*60/(2*npy.pi))**(-0.83)*(self.Dpw*1e3)**(-0.5)
        else:
            nu1=4500*(n*60/(2*npy.pi))**(-0.5)*(self.Dpw*1e3)**(-0.5)
        
        coeff_oil=self.O1.oil_kinematic_viscosity_curve[self.oil_name]
        nu=10**(10**(coeff_oil['A']*npy.log10(T)+coeff_oil['B']))-0.6
        self.kappa=nu/nu1
        #définition du paramètre de contamination
        self.O1.OilContamination(self.Dpw)
        ec=self.O1.ec
        #calcul de la limite de charge en fatigue
        self.BaseStaticLoad()
        if self.Dpw<0.1:
            Cu=self.C0r/8.2
        else:
            Cu=self.C0r/8.2*(100/(self.Dpw*1e3))**0.3
        #calcul du coefficient a_iso
        self.EquivalentDynamicLoad(Fr,Fa)
        if self.kappa<0.4:
            self.a_iso=0.1*((1-(1.5859-1.3993/(self.kappa**0.054381))*((ec*Cu/self.Pr)**0.4))**(-9.185))
        elif self.kappa<1:
            self.a_iso=0.1*((1-(1.5859-1.2348/(self.kappa**0.19087))*((ec*Cu/self.Pr)**0.4))**(-9.185))
        else:
            self.a_iso=0.1*((1-(1.5859-1.2348/(self.kappa**0.071739))*((ec*Cu/self.Pr)**0.4))**(-9.185))
        self.a_iso=min(50,self.a_iso)
        #calcul de la durée de vie corrigée
        self.BaseLifeTime(Fr,Fa)
        self.Lnm=self.a1*self.a_iso*self.L10
        
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
        del d['O1']

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
        elif self.typ=='N':
            
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
        return ref
    
    def RollerContour(self):
        p=[vm.Point2D((0,-self.Dw/2))]
        p.append(vm.Point2D((-self.Lw/2,-self.Dw/2)))
        p.append(vm.Point2D((-self.Lw/2,self.Dw/2)))
        p.append(vm.Point2D((self.Lw/2,self.Dw/2)))
        p.append(vm.Point2D((self.Lw/2,-self.Dw/2)))
        p.append(p[0])
        ref=vm.Contour2D(primitives2D.RoundedLines2D(p,{1:self.r_roller,2:self.r_roller,3:self.r_roller,4:self.r_roller},False).primitives)
        return ref
        
    def FreeCADExport(self,file_path,export_types):
        #bague interne
        IRC=self.InternalRingContour()        
        irc=primitives3D.RevolvedProfile(vm.Point3D((0,0,0)),vm.Vector3D((0,0,1)),
                                           vm.Vector3D((0,1,0)),[IRC],vm.Vector3D((0,0,0)),
                                           vm.Vector3D((0,0,1)),angle=2*math.pi,name='irc')
        #bague externe
        ERC=self.ExternalRingContour()
        erc=primitives3D.RevolvedProfile(vm.Point3D((0,0,0)),vm.Vector3D((0,0,1)),
                                           vm.Vector3D((0,1,0)),[ERC],vm.Vector3D((0,0,0)),
                                           vm.Vector3D((0,0,1)),angle=2*math.pi,name='erc')
        #roller
        ROL=self.RollerContour()
        radius=self.F/2+self.jeu+self.Dw/2
        rol=[]
        theta=2*npy.pi/self.Z
        for z in range(int(self.Z)):
            rol.append(primitives3D.RevolvedProfile(vm.Point3D((radius*npy.sin(z*theta),radius*npy.cos(z*theta),0)),vm.Vector3D((0,0,1)),
                                               vm.Vector3D((radius*npy.sin(z*theta),radius*npy.cos(z*theta),0)),[ROL],vm.Vector3D((radius*npy.sin(z*theta),radius*npy.cos(z*theta),0)),
                                               vm.Vector3D((0,0,1)),angle=2*math.pi,name='rol'))
        
#        total=[ROL]
#        G1=vm.Contour2D(total)
#        G1.MPLPlot()
        
        tot=[irc,erc]+rol
        model=vm.VolumeModel(tot)
#        model=vm.VolumeModel([gear1,t1,gear2])
        model.FreeCADExport('python',file_path,'/usr/lib/freecad/lib',export_types)
                
# =============================================================
# Objet avec 3 fonctions de selection des roulements cylindriques
#   - Combinatoire sur les dimensions externe ISO
#   - Combinatoire en prenant en compte les règles SKF
#   - Estimation des durées de vie et charge dynamique et fonction de tri
# =============================================================
        
class BearingCombination():
    def __init__(self):
        self.tableau_serie=pandas.read_csv('../mechanical_components/catalogs/serie_rlts_iso.csv')
        self.roller=pandas.read_csv('../mechanical_components/catalogs/roller_iso.csv')
        self.radial_clearance=pandas.read_csv('../mechanical_components/catalogs/radial_clearance_iso.csv')
        df1=self.tableau_serie.copy()
        df2=self.roller.copy()
        df3=self.radial_clearance.copy()
        self.df=[df1,df2,df3]
        self.dic={}
        for i,it in enumerate(self.df):
            for d in it.columns:
                self.dic[d]=i
        
    def LoadSKFRules(self):
        self.rules_rlts_skf=pandas.read_csv('../mechanical_components/catalogs/rules_rlts_SKF.csv')
        
    def Analyze(self,limit,Fr,n,grade=['Gr_gn'],Fa=0):
        # Combinatoire sur les dimensions externe ISO
        self.Fr=Fr
        self.Fa=Fa
        self.n=n
        data_rlts=self.tableau_serie.copy()
        for (k1,v1) in limit.items():
            #choix des rlts
            data_rlts=data_rlts[(data_rlts[k1] >= v1[0]) & (data_rlts[k1] <= v1[1])]
        liste1=list(data_rlts.index)
        liste2=[]
        for i,index in enumerate(liste1):
            data_roll=self.roller.copy()
            #choix des rouleaux
            data_roll=data_roll[(data_rlts.D[index]-data_rlts.d[index])/2 > data_roll.Dw]
            data_roll=data_roll[data_rlts.B[index] > data_roll.Lw]
            for j in list(data_roll.index):
                liste2.append((index,j))
            
        liste3=[]
        data_clearance1=self.radial_clearance.copy()
        for i,ind in enumerate(liste2):
            a=data_rlts.d[ind[0]]
            data_clearance=data_clearance1
            for j,gr in enumerate(grade):
                liste_grade=['d_min','d_max']
                liste_grade.append(gr+'_min')
                liste_grade.append(gr+'_max')

                data_clearance=data_clearance[(a > data_clearance.d_min) & (a <= data_clearance.d_max)]
                
                for k in data_clearance.index:
                    liste3.append((ind[0],ind[1],k,gr))
#                matrix_concat=npy.concatenate((npy.tile(data_rlts.loc[ind,:].values,(len(data_clearance),1)),data_clearance.as_matrix()),axis=1)
#                tableau_serie_min_temp=pandas.DataFrame(matrix_concat,columns=list(data_rlts.columns)+list(data_clearance.columns))
        return liste3
    
    def AnalyseSKFRules(self,liste):
        # Combinatoire en prenant en compte les règles SKF
        df_rules=self.rules_rlts_skf.to_dict()
        liste_out=[]
        for item in liste:
            drap=1
            for k2,v2 in df_rules['type'].items():
                var_x=df_rules['x'][k2]
                var_y=df_rules['y'][k2]
                a=df_rules['a'][k2]
                b=df_rules['b'][k2]
                if (var_x in self.dic.keys()) & (var_y in self.dic.keys()):
                    typ=df_rules['type'][k2]
                    ind_x=self.dic[var_x]
                    d_x=self.df[ind_x][var_x][item[ind_x]]
                    ind_y=self.dic[var_y]
                    d_y=self.df[ind_y][var_y][item[ind_y]]
                    if typ=='inf':
                        if d_y<(a*d_x+b)*0.8:
                            drap=0
                    elif typ=='sup':
                        if d_y>(a*d_x+b)*1.2:
                            drap=0
            F_inter=self.AnalyseSKFInterRules(item,'F')
            Dw=self.AccesData('Dw',item)
            D=self.AccesData('D',item)
            if F_inter[0]>F_inter[1]:
                drap=0
            E_min=F_inter[0]+2*Dw
            if E_min>D:
                drap=0
            if drap==1:
                liste_out.append(item)
        return liste_out
    
    def AnalyseSKFValueRules(self,var_x,var_y,data_x):
        df_rules=self.rules_rlts_skf.to_dict()
        for k1,v1 in df_rules['type'].items():
            if (var_y==df_rules['y'][k1]) & (var_x==df_rules['x'][k1]):
                data_y=df_rules['a'][k1]*data_x+df_rules['b'][k1]
        return data_y
    
    def AnalyseSKFInterRules(self,item,var):
        # définition pour une variable donnée "var" de l'intervalle d'existance de cette variable 
        # pour les données du roulement défini par "item" (liste des adresses dans les catalogues ISO)
        df_rules=self.rules_rlts_skf.to_dict()
        borne_inf=-npy.inf
        borne_sup=npy.inf
        for k1,v1 in df_rules['type'].items():
            if var==df_rules['y'][k1]:
                if df_rules['type'][k1]=='inf':
                    if df_rules['x'][k1] in self.dic.keys():
                        varx=df_rules['x'][k1]
                        ind=self.dic[varx]
                        d1=self.df[ind][varx][item[ind]]
                        borne_inf=max(borne_inf,df_rules['a'][k1]*d1+df_rules['b'][k1])
                if df_rules['type'][k1]=='sup':
                    if df_rules['x'][k1] in self.dic.keys():
                        varx=df_rules['x'][k1]
                        ind=self.dic[varx]
                        d1=self.df[ind][varx][item[ind]]
                        borne_sup=min(borne_sup,df_rules['a'][k1]*d1+df_rules['b'][k1])
        return [borne_inf,borne_sup]
    
    def AccesData(self,var,item):
        return self.df[self.dic[var]][var][item[self.dic[var]]]
        
    def SortBearing(self,liste,const,S,T,oil_name,nb_sol,typ):
        # Estimation des durées de vie et charge dynamique et fonction de tri
        for ind,item in enumerate(liste):
            Dw=self.AccesData('Dw',item)
            Lw=self.AccesData('Lw',item)
            D=self.AccesData('D',item)
            d=self.AccesData('d',item)
            B=self.AccesData('B',item)
            F_inter=self.AnalyseSKFInterRules(item,'F')
            Gr_min=self.AccesData(item[-1]+'_min',item)
            Gr_max=self.AccesData(item[-1]+'_max',item)
            Gr=(Gr_min+Gr_max)/2
            rsmin=self.AccesData('rsmin',item)
            rsmax=self.AccesData('rsmax',item)
            r_roller=(rsmin+rsmax)/2
            
            data={'typ':typ,'B':B,'d':d,'D':D,'Lw':Lw,'Dw':Dw,'r_roller':r_roller,'i':1,'alpha':0,'O1':self.O1,'oil_name':oil_name}
            liste_F=npy.arange(F_inter[0],F_inter[1],(F_inter[1]-F_inter[0])/10)
#            sol=npy.array([[],[],[],[],[],[],[],[]])
            for f in liste_F:
                Zmax=int(2*npy.pi/(2*npy.arcsin((Dw/2)/(f/2+Dw/2))))
                E=f+2*Dw+Gr
                data['F']=f
                data['E']=E
                
                d1=self.AnalyseSKFValueRules('F','d1',f)
                D1=self.AnalyseSKFValueRules('E','D1',E)
                if typ=='NU':
                    data['d1']=f
                else:
                    data['d1']=d1
                if typ=='N':
                    data['D1']=E
                else:
                    data['D1']=D1
                    
                for Z in range(int(3/4*Zmax),Zmax+1):
                    data['Z']=Z
                    R1=RadialRollerBearing(**data)
                    R1.BaseDynamicLoad()
                    R1.BaseStaticLoad()
                    R1.BaseLifeTime(Fr=self.Fr)
                    R1.AdjustedLifeTime(Fr=self.Fr,n=self.n,Fa=self.Fa,S=S,T=T)
                    data_export=data.copy()
                    data_export['C0r']=R1.C0r
                    data_export['Cr']=R1.Cr
                    data_export['L10']=R1.L10
                    data_export['Lnm']=R1.Lnm
                    ll=[[ind]]
                    ll_pos=['B','d','D','d1','D1','Lw','Dw','r_roller','i','alpha','E','F','Z','C0r','Cr','L10','Lnm']
                    for key in ll_pos:
                        ll.append([data_export[key]])
                    try:
                        sol=npy.concatenate((sol,ll),axis=1)
                    except UnboundLocalError:
                        sol=npy.array(ll)
                    
            li=self._AnalyseOptim(sol,const,ll_pos,nb_sol=nb_sol+20)
            try:
                liste_sort=npy.concatenate((liste_sort,sol[:,li]),axis=1)
            except UnboundLocalError:
                liste_sort=npy.array(sol[:,li])
        li=self._AnalyseOptim(liste_sort,const,ll_pos,nb_sol=nb_sol+20)
        self.solution=[]
        for i in range(nb_sol):
            for k1 in data.keys():
                try:
                    data[k1]=liste_sort[ll_pos.index(k1)+1,li[i]]
                except ValueError:
                    pass
            R1=RadialRollerBearing(**data)
            R1.BaseDynamicLoad()
            R1.BaseStaticLoad()
            R1.BaseLifeTime(Fr=self.Fr)
            R1.AdjustedLifeTime(Fr=self.Fr,n=self.n,Fa=self.Fa,S=S,T=T)
            self.solution.append(R1)

    def _AnalyseOptim(self,matrix,const,ll,nb_sol=20):
        matrix1=matrix
        temp=matrix1.argsort()
        li_k2=[]
        for co in const:
            pos=ll.index(co['var'])
            if co['type']=='min':
                k2=list(temp[pos+1,0:min(nb_sol,temp.shape[1])])
                matrix1=matrix1[:,tuple(k2)]
                temp=matrix1.argsort()
            elif co['type']=='max':
                k2=list(temp[pos+1,temp.shape[1]-min(nb_sol,temp.shape[1]):-1])
                k2.reverse()
                matrix1=matrix1[:,tuple(k2)]
                temp=matrix1.argsort()
            elif co['type']=='bound_sup':
                val=co['val']
                k2=list(npy.where((matrix1[pos+1,:])<=val)[0])
                matrix1=matrix1[:,k2]
                temp=matrix1.argsort()
            elif co['type']=='bound_inf':
                val=co['val']
                k2=list(npy.where((matrix1[pos+1,:])>=val)[0])
                matrix1=matrix1[:,k2]
                temp=matrix1.argsort()
            li_k2.append(tuple(k2))
        tp=range(matrix.shape[1])
        for li in li_k2:
            tpo=[]
            for l in li:
                tpo.append(tp[l])
            tp=tpo
        return tp
    
    def OptimizerBearing(self,input_dict):
#        limit,grade,Fr,n=100,Fa=0,S=0.9,T=40,oil_name='iso_vg_100',nb_sol=10,maxi=None,mini=None
        maxi=None
        mini=None
        
        self.O1=OilData()
        self.LoadSKFRules()
        limit_ISO={}
        limit_sort=[]
        for it in input_dict:
            if 'min' in it.keys():
                if it['type'] in ['d','D','B','rsmin']:
                    limit_ISO[it['type']]=[it['min'],it['max']]
                elif it['type'] in ['L10','C0r','Cr','Lnm']:
                    limit_sort.append({'type':'bound_inf','var':it['type'],'val':it['min']})
                    limit_sort.append({'type':'bound_sup','var':it['type'],'val':it['max']})
            elif 'grade'==it['type']:
                grade=[it['nom']]
            elif 'Fr'==it['type']:
                Fr=it['nom']
            elif 'Fa'==it['type']:
                Fa=it['nom']
            elif 'n'==it['type']:
                n=it['nom']
            elif 'S'==it['type']:
                S=it['nom']
            elif 'T'==it['type']:
                T=it['nom']
            elif 'oil_name'==it['type']:
                oil_name=it['nom']
            elif 'nb_sol'==it['type']:
                nb_sol=it['nom']
            elif 'maxi'==it['type']:
                maxi=it['nom']
            elif 'mini'==it['type']:
                mini=it['nom']
            elif 'typ'==it['type']:
                typ=it['nom']

        liste_ind=self.Analyze(limit=limit_ISO,grade=grade,Fr=Fr,Fa=Fa,n=n)
        print('Nombre de Solution ISO: ',len(liste_ind))
        liste_ind=self.AnalyseSKFRules(liste_ind)
        print('Nombre de Solution avec règles SKF: ',len(liste_ind))
        
        limit_sort.extend([{'type':'min','var':'B','val':None},
                            {'type':'min','var':'D','val':None},
                            {'type':'max','var':'d','val':None}])
        if not maxi==None:
            limit_sort.append({'type':'max','var':maxi,'val':None})
        elif not mini==None:
            limit_sort.append({'type':'min','var':mini,'val':None})
            
        self.SortBearing(liste_ind,const=limit_sort,S=S,T=T,oil_name=oil_name,nb_sol=nb_sol,typ=typ)
                    
    
# =============================================================
# Test sur un objet RadialRollerBearing
# =============================================================

        



