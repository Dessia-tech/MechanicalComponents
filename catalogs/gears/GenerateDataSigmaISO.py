#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 11:24:40 2018

@author: Pierrem
"""
import numpy as npy
import math as mt
from scipy import interpolate
#import os
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math
from scipy.linalg import norm
from scipy.optimize import minimize,fsolve
from scipy.interpolate import splprep, splev
#import cma
#from sympy import *
import itertools
from jinja2 import Environment, PackageLoader, select_autoescape
import networkx as nx
#from interval import interval, inf, imath

import mechanical_components.LibSvgD3 as LibSvg
import powertransmission.tools as tools

import persistent


class Material():
    
    def __init__(self):
        
        data_array,dico_axis,dico_nom,type_axis=self.DataCoeffYBIso()
        self.data_coeff_YB_Iso=self.GenereCoeff(data_array,dico_axis,dico_nom,type_axis)
    
        data_array,dico_axis,dico_nom,type_axis=self.DataWholerCurve()
        self.data_wholer_curve=self.GenereCoeff(data_array,dico_axis,dico_nom,type_axis)
        
        data_array,dico_axis,dico_nom,type_axis=self.DataGearMaterial()
        self.data_gear_material=self.GenereCoeff(data_array,dico_axis,dico_nom,type_axis)
        
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
    
    def DataCoeffYBIso(self):
        data_svg=['56.097561,477.97196 81.707319,-1.21951 81.70732,1.21951 82.92682,0 79.2683,-1.21951 82.92683,0 81.70731,-1.21951 81.70732,0 81.70732,1.21951',
                               '57.317073,476.75245 -1.219512,-82.92683 0,-82.92683 1.219512,-84.14634 0,-82.92683 1.219512,-82.92683',
                               '56.097561,58.459766 76.829269,59.756094 71.95122,43.90244 89.02439,45.12195 78.04878,32.92683 75.60976,19.5122 74.39024,7.31707 90.2439,1.21951 97.56098,1.21952']
        dico_nom={'evol_coeff_yb_iso':2,'x':0,'y':1}
        dico_axis={'x':[0,5,10,15,20,25,30,35,40],'y':[0.5,0.6,0.7,0.8,0.9,1]}
        type_axis={'x':'Linear','y':'Linear'}
        data_array=self.SVG2Array(data_svg)
        return data_array,dico_axis,dico_nom,type_axis
    
    def DataWholerCurve(self):
        data_svg=['199.46895,423.19149 0,-42.67247 -0.99238,-41.68007 0.99238,-31.75625 -0.99238,-28.77911 0,-20.84004 0,-19.84765 0,-20.84004 0,-18.85528 0.99238,-17.86289 0,-14.88574 0,-16.87051 0,-12.90098 0,-9.92382 0,-11.9086 0,-11.90859 0,-10.916213','200.46133,423.19149 34.7334,-0.99239 20.84004,0.99239 14.88574,-0.99239 11.9086,0 36.29445,0.74708 36.12719,-0.80283 21.67632,0.80283 14.45088,0 11.23957,0 38.93709,-0.40141 36.1272,0.40141 21.67631,0 16.85936,-0.40141 11.23957,0 36.1272,-0.40142 37.33143,-0.40141 19.26784,0.40141 14.04947,0.80283 9.2325,0 36.93002,-0.40141','236,100.3622 63.42857,18.28572 61.14286,15.42857 93.71428,21.71429 57.71429,3.42857 65.71429,1.14285 65.71428,0.57143 38.28572,0','213.04348,126.27525 55.65217,34.78261 38.26087,22.60869 66.08696,6.95652 115.65217,-0.86956 103.47826,0 86.95652,-0.86957','221.14286,140.3622 60.57143,46.28572 34.85714,26.85714 43.42857,27.42857 13.14286,5.71429 78.28571,-0.57143 125.71429,1.14286 103.42857,-2.28572','234.28571,204.3622 50.28572,29.14286 37.71428,21.14286 41.14286,21.14286 31.42857,12 48,5.14285 65.71429,0 83.42857,-1.14285 88,0','237.3913,261.05786 112.17392,50.43478 83.47826,31.30435 69.56522,18.26087 84.34782,-0.86957 92.17391,-0.86956']
        dico_nom={'hardened_alloy_steel':2,'nitrided_alloy_steel':3,'through_hardened_steel':4,'surface_hardened_steel':5,'x':1,'y':0,'carbon_steel':6,'cast_iron':6,'bronze':6,'grey_iron':6}
        dico_axis={'x':[10000,100000000],'y':[20,100]}
        type_axis={'x':'Log','y':'Log'}
        data_array=self.SVG2Array(data_svg)
        return data_array,dico_axis,dico_nom,type_axis
    
    def DataGearMaterial(self):
        data_svg=['76.2931,579.71715 91.17956,0.9304 59.54583,0.9304 47.45059,-0.9304 134.90853,0 110.71804,-0.9304','76.2931,579.71715 0,-148.86459 -0.930404,-90.24915 0,-240.974553','350.76218,194.53003 53.96341,-28.84251 56.75463,-29.77292 40.93776,-22.32969','352.62299,209.41649 49.31139,-26.98171 39.07695,-19.53847 54.89382,-29.77292','319.12846,281.05757 71.64108,-46.52018 56.75462,-39.07696 66.98906,-45.58977','345.17976,370.37632 51.1722,-42.79857 53.96341,-42.79856','267.02585,377.81955 43.72897,-36.28574 38.14655,-33.49453 23.26009,-20.46888','163.75104,467.1383 37.21615,-30.70332 30.70332,-25.1209 16.74727,-17.67767','82.805926,535.98817 38.146554,-31.63372 20.46888,-15.81686 13.02565,-14.88646','103.27481,574.13472 37.21614,-33.49453 33.49453,-35.35534 25.1209,-26.0513 14.88646,-17.67767']
        dico_nom={'hardened_alloy_steel':2,'nitrided_alloy_steel':3,'through_hardened_steel':4,'surface_hardened_steel':5,'carbon_steel':6,'cast_iron':7,'bronze':8,'grey_iron':9,'x':0,'y':1}
        dico_axis={'x':[20,150],'y':[5,45]}
        type_axis={'x':'Log','y':'Log'}
        data_array=self.SVG2Array(data_svg)
        return data_array,dico_axis,dico_nom,type_axis
    
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
            for item in sol:
                temp=item.split(',')
                data_temp.append([float(temp[0]),float(temp[1])])
            export_temp=[data_temp[0]]
            for item in data_temp[1::]:
                export_temp.append([item[0]+export_temp[-1][0],item[1]+export_temp[-1][1]])
            export[i]=npy.array(export_temp)
        return export
    
M1=Material()
fichier=open('GenerateDataSigmaISO.txt','w')

fichier.write("#data_coeff_YB_Iso\n")
for k,v in M1.data_coeff_YB_Iso.items():
    if k not in ['x','y']:
        fichier.write(k+'={')
        if 'str' not in str(v.__class__):
            fichier.write('\'data\':[')
            for ligne in v[:-1]:
                fichier.write('[{},{}]'.format(*tuple(ligne)))
                fichier.write(',')
            fichier.write('[{},{}]'.format(*tuple(v[-1])))
            fichier.write(']')
        else:
            fichier.write('\'{}\''.format(v))
        fichier.write(',\'x\':\'{}\''.format(M1.data_coeff_YB_Iso['x']))
        fichier.write(',\'y\':\'{}\''.format(M1.data_coeff_YB_Iso['y']))
fichier.write('}\n')    

fichier.write("#data_wholer_curve\n")
for k,v in M1.data_wholer_curve.items():
    if k not in ['x','y']:
        fichier.write('wholer_'+k+'={')
        if 'str' not in str(v.__class__):
            fichier.write('\'data\':[')
            for ligne in v[:-1]:
                fichier.write('[{},{}]'.format(*tuple(ligne)))
                fichier.write(',')
            fichier.write('[{},{}]'.format(*tuple(v[-1])))
            fichier.write(']')
        else:
            fichier.write('\'{}\''.format(v))
        fichier.write(',\'x\':\'{}\''.format(M1.data_wholer_curve['x']))
        fichier.write(',\'y\':\'{}\''.format(M1.data_wholer_curve['y']))
        fichier.write('}\n') 

fichier.write("#data_gear_material\n")
for k,v in M1.data_gear_material.items():
    if k not in ['x','y']:
        fichier.write('sigma_'+k+'={')
        if 'str' not in str(v.__class__):
            fichier.write('\'data\':[')
            for ligne in v[:-1]:
                fichier.write('[{},{}]'.format(*tuple(ligne)))
                fichier.write(',')
            fichier.write('[{},{}]'.format(*tuple(v[-1])))
            fichier.write(']')
        else:
            fichier.write('\'{}\''.format(v))
        fichier.write(',\'x\':\'{}\''.format(M1.data_gear_material['x']))
        fichier.write(',\'y\':\'{}\''.format(M1.data_gear_material['y']))
        fichier.write('}\n') 

fichier.close()