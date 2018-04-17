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

import persistent
import pandas
from pandas.plotting import scatter_matrix
import matplotlib.pyplot as plt
#from dessia_common import ResultsDBClient
#import pyDOE

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

class RadialRollerBearing(persistent.Persistent):
    def __init__(self,B,d,D,Lw,Dw,r_roller,E,F,Z,i,alpha,weibull_e=9/8,weibull_c=31/3,weibull_h=7/3,B1=551.13373/0.483,mu_delta=0.83,bm=1.1,c_gamma=0.05,oil_name='iso_vg_100'):
        self.B=B
        self.d=d
        self.D=D
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
        O1=OilData()
        coeff_oil=O1.oil_kinematic_viscosity_curve[self.oil_name]
        nu=10**(10**(coeff_oil['A']*npy.log10(T)+coeff_oil['B']))-0.6
        self.kappa=nu/nu1
        #définition du paramètre de contamination
        O1.OilContamination(self.Dpw)
        ec=O1.ec
        #calcul de la limite de charge en fatigue
        self.BaseStaticLoad()
        if self.Dpw<0.1:
            Cu=self.C0r/8.2
        else:
            Cu=self.C0r/8.2*(100/(self.Dpw*1e3))**0.3
        print(Cu)
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

        
        
#R1=RadialRollerBearing(B=0.02,d=0.03,D=0.045,Lw=0.015,Dw=0.005,r_roller=0.001,E=0.04,F=0.036,Z=30,i=1,alpha=0)
#R1.BaseLifeTime(1000)
#R1.AdjustedLifeTime(1000,500,T=80,S=0.99)
                
        
class BearingCombination():
    def __init__(self,database):
        self.database=database
        
    def Analyze(self,limit):
        #analyse dimension generale du rlts
        data_rlts=self.database.tableau_serie.copy()
        for (k1,v1) in limit.items():
            #choix des rlts
            data_rlts=data_rlts[(data_rlts[k1] >= v1[0]) & (data_rlts[k1] <= v1[1])]
        for i,index in enumerate(data_rlts.index):
            data_roll=self.database.roller.copy()
            #choix des rouleaux
            data_roll=data_roll[(data_rlts.D[index]-data_rlts.d[index])/2 > data_roll.Dw]
            data_roll=data_roll[data_rlts.B[index] > data_roll.Lw]
            
            matrix_concat=npy.concatenate((npy.tile(data_rlts.loc[index,:].values,(len(data_roll),1)),data_roll.as_matrix()),axis=1)
            tableau_serie_min_temp=pandas.DataFrame(matrix_concat,columns=list(data_rlts.columns)+list(data_roll.columns))
            if i==0:
                tableau_serie_min=tableau_serie_min_temp
            if i >0:
                tableau_serie_min=pandas.concat([tableau_serie_min, tableau_serie_min_temp])
                
        return tableau_serie_min
        
#D1=DataBaseBearing()
#C1=BearingCombination(D1)
#tab_selection=C1.Analyze(limit={'d':[0.015,0.016],'D':[0.035,0.040],'B':[0.011,0.015]})
#print(len(tab_selection))
