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

def GenereCoeff(data_array,dico_axis,dico_nom,type_axis):
    structure={}
    vect_x=npy.linspace(data_array[dico_nom['x']][0,0],data_array[dico_nom['x']][-1,0],len(data_array[dico_nom['x']][:,0]))
    vect_y=npy.linspace(data_array[dico_nom['y']][0,1],data_array[dico_nom['y']][-1,1],len(data_array[dico_nom['y']][:,1]))
    
    for i,j in type_axis.items():
        if j=='Log':
            axe_reel=[mt.log10(dico_axis[i][0]),mt.log10(dico_axis[i][-1])]
        if j=='Linear':
            axe_reel=dico_axis[i]
        if i=='x':
            ax,bx=AxisLinear(axe_reel,vect_x)
        if i=='y':
            ay,by=AxisLinear(axe_reel,vect_y)
    
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
    
def DataSNR():
    data_svg=['m 203.14428,545.55955 53.9602,1.05804 55.01824,0 51.84411,-2.11608 53.9602,2.11608 53.9602,0 55.01824,0','m 203.14428,546.61759 1.05804,-105.80431 -1.05804,-105.80431 0,-104.74627 0,-104.74627','m 216.89884,185.82489 116.38474,209.49254 191.50581,62.42454','M 213.72471,306.44181 333.28358,446.1035 524.78939,489.48327','m 209.49254,422.82655 124.84909,74.06302 191.5058,20.10282']
    dico_nom={'x':0,'y':1,'punctual_load':2,'mixed_load':3,'constant_load':4}
    dico_axis={'x':[0,5,10,15,20,25,30],'y':[0,50,100,150,200]}
    type_axis={'x':'Linear','y':'Linear'}
    data_array=SVG2Array(data_svg)
    return data_array,dico_axis,dico_nom,type_axis
    
    
def FunCoeff(x,data,type_x='Linear',type_y='Linear'):
    if type_x=='Log': 
        x=npy.log10(x)
    f = interpolate.interp1d(list(data[:,0]),list(data[:,1]), fill_value='extrapolate')
    sol=float(f(x))
    if type_y=='Log':
        sol=10**sol
    return sol

def AxisLinear(axe,vect):
    a=(axe[0]-axe[-1])/(vect[0]-vect[-1])
    b=axe[0]-a*vect[0]
    return a,b

def SVG2Array(data):
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

data_array,dico_axis,dico_nom,type_axis=DataSNR()
data_snr=GenereCoeff(data_array,dico_axis,dico_nom,type_axis)

fichier=open('RadialRollerBearing_MaxAxialLoad.txt','w')
#
fichier.write("#RadialRollerBearing_MaxAxialLoad\n")
for k,v in data_snr.items():
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
        fichier.write(',\'x\':\'{}\''.format(data_snr['x']))
        fichier.write(',\'y\':\'{}\''.format(data_snr['y']))
        fichier.write('}\n')    

fichier.close()