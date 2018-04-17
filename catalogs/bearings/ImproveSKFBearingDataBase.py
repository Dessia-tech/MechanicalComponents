#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 11:47:48 2018

@author: dumouchel@dessia.tech
"""

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
import mechanical_components.bearing as bearing

import persistent
import pandas
from pandas.plotting import scatter_matrix

import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
from statsmodels.regression.quantile_regression import QuantReg

# =============================================================
# Object qui permet de completer la base de rlts SFK
# =============================================================

chemin_catalogs='../../mechanical_components/catalogs/'
class DataBaseBearing():
    def __init__(self):
        self.tableau_serie=pandas.read_csv(chemin_catalogs+'serie_rlts_iso.csv')
        self.roller=pandas.read_csv(chemin_catalogs+'roller_iso.csv')
        self.radial_clearance=pandas.read_csv(chemin_catalogs+'radial_clearance_iso.csv')
    def AddColumnDataBase(self,data,column=None):
        print(str(data.__class__))
        if data.__class__==float:
            self.database[column]=data
        elif 'pandas' in str(data.__class__):
            self.database=pandas.concat([self.database, data], axis=1,verify_integrity=True)
        return self.database
    def UpdateSKFDataBase(self,database=None):
        if 'pandas' in str(database.__class__):
            self.database=database
        self.database.to_csv(chemin_catalogs+'tableau_rlts_SKF.csv',index=False)
        return self.database
    def LoadSKFDataBase(self):
        self.database=pandas.read_csv(chemin_catalogs+'tableau_rlts_SKF.csv')
        return self.database
    def ImproveSKFDataBase(self):
        bearing_SKF=pandas.read_csv(chemin_catalogs+'tableau_rlts_SKF_brut.csv')
#        bearing_SKF['F']=pandas.Series(npy.NaN)
#        bearing_SKF['E']=pandas.Series(npy.NaN)
        bearing_SKF['type']=pandas.Series(npy.NaN)
        for ind in bearing_SKF.name.index:
            item=bearing_SKF.name[ind]
            if item[0]=='N':
                typ=item.split('_')[0]
                bearing_SKF.loc[ind,'type']=typ
            else:
                bearing_SKF.loc[ind,'type']='None'
        
        bearing_SKF_min=[]
        for index in (bearing_SKF.name).index:
            D=bearing_SKF['D'][index]
            d=bearing_SKF['d'][index]
            B=bearing_SKF['B'][index]
            item=bearing_SKF['type'][index]
            name=bearing_SKF['name'][index]
            
            if item in ['NU','N','NJ','NUP']:
                if item=='NU':
                    D1=bearing_SKF['D1'][index]
                    F=bearing_SKF['EF'][index]
                elif item=='N':
                    d1=bearing_SKF['d1'][index]
                    E=bearing_SKF['EF'][index]
                elif item in ['NJ','NUP']:
                    d1=bearing_SKF['d1'][index]
                    D1=bearing_SKF['D1'][index]
                    F=bearing_SKF['EF'][index]
                
                data_clearance=self.radial_clearance.copy()
                data_clearance=data_clearance[d >= data_clearance.d_min]
                data_clearance=data_clearance[d <  data_clearance.d_max]
                clearance=((data_clearance.Gr_gn_min+data_clearance.Gr_gn_max)/2).values[0]
                
                data_roll=self.roller.copy()
                if item in ['NU','NJ','NUP']:
                    data_roll=data_roll[(D-F)/2 > data_roll.Dw]
                    data_roll=data_roll[(D1-F)/2 < data_roll.Dw]
                elif item=='N':
                    data_roll=data_roll[(E-d)/2 > data_roll.Dw]
                    data_roll=data_roll[(E-d1)/2 < data_roll.Dw]
                data_roll=data_roll[B > data_roll.Lw]
                
                for ind in data_roll.index:
                    Lw=data_roll.Lw[ind]
                    Dw=data_roll.Dw[ind]
                    rsmin=data_roll.rsmin[ind]
                    rsmax=data_roll.rsmax[ind]
                    if item in ['NU','NJ','NUP']:
                        E=F+2*Dw+clearance
                    elif item=='N':
                        F=E-(2*Dw+clearance)
                    ratio_radial=((F-d)/2)/((D-E)/2)
                    ratio_axial=(B-Lw)/B
                    Zmax=int(2*npy.pi/(2*npy.arcsin((Dw/2)/(F/2+Dw/2))))
                    bearing_temp=bearing.RadialRollerBearing(D=D,d=d,B=B,Lw=Lw,Dw=Dw,r_roller=(rsmin+rsmax)/2,E=E,F=F,Z=Zmax,i=1,alpha=0)
                    bearing_temp.BaseStaticLoad()
                    bearing_temp.BaseDynamicLoad()
                    Z=Zmax
                    solution=[]
                    while bearing_temp.C0r>bearing_SKF['C0'][index]:
                        solution.append([ratio_radial,ratio_axial,Z+1,Lw,Dw,(rsmin+rsmax)/2,E,F,Zmax,bearing_temp.C0r,bearing_temp.Cr])
                        Z=Z-1
                        bearing_temp=bearing.RadialRollerBearing(D=D,d=d,B=B,Lw=Lw,Dw=Dw,r_roller=(rsmin+rsmax)/2,E=E,F=F,Z=Z,i=1,alpha=0)
                        bearing_temp.BaseStaticLoad()
                        bearing_temp.BaseDynamicLoad()
#                    if (bearing_temp.Cr>0.9*bearing_SKF['C'][index]) and (bearing_temp.Cr<1.1*bearing_SKF['C'][index]):
                        
                    
                    if not solution==[]:
                        post_solution=pandas.DataFrame(solution,columns=['ratio_r','ratio_a','Z','Lw','Dw','r_roller','E','F','Zmax','C0r','Cr'])
                        post_solution=post_solution[(post_solution.ratio_a<0.45) & (post_solution.ratio_a>0.1)]
                        post_solution=post_solution[(post_solution.ratio_r<1.4) & (post_solution.ratio_r>0.4)]
                        if len(post_solution)>0:
                            ind_max=post_solution.sort_values(by = ['ratio_r', 'ratio_a']).index[-1]
                            matrix_temp=npy.concatenate((bearing_SKF.loc[index,:].values,post_solution.loc[ind_max,:].values),axis=0)
                            try:
                                matrix_concat=npy.concatenate((matrix_concat,[matrix_temp]),axis=0)
                            except:
                                matrix_concat=[matrix_temp]
            try:
                print(item,len(matrix_concat),bearing_SKF.name[index])
            except:
                pass
        bearing_SKF_min=pandas.DataFrame(matrix_concat,columns=list(bearing_SKF.columns)+list(post_solution.columns))
        matrix_concat=[]
        return bearing_SKF_min
    
    def QuantileRegression(self,var1,var2,dataframe,out):
        
        data=pandas.DataFrame(dataframe, columns=[var1,var2])
        data.head()
        mod = smf.quantreg(var1+' ~ '+var2, data)
        res = mod.fit(q=.5)
        quantiles = npy.arange(.01, .99, .01)
        def fit_model(q):
            res = mod.fit(q=q)
            return [q, res.params['Intercept'], res.params[var2]] + \
                    res.conf_int().ix[var2].tolist()
            
        models = [fit_model(x) for x in quantiles]
        models = pandas.DataFrame(models, columns=['q', 'a', 'b','lb','ub'])
        
        ols = smf.ols(var1+' ~ '+var2, data).fit()
        ols_ci = ols.conf_int().ix[var2].tolist()
        ols = dict(a = ols.params['Intercept'],
                   b = ols.params[var2],
                   lb = ols_ci[0],
                   ub = ols_ci[1])
        
        if out==1:
            x = npy.arange(data[var2].min(), data[var2].max(),(data[var2].max()-data[var2].min())/250)
            get_y = lambda a, b:b * x + a
            
            fig, ax = plt.subplots(figsize=(8, 6))
            
            for i in [0,models.shape[0]-1]:
                y = get_y(models.a[i], models.b[i])
                ax.plot(x, y, linestyle='dotted', color='grey')
                
            y = get_y(ols['a'], ols['b'])
            
            ax.plot(x, y, color='red', label='OLS')
            ax.scatter(data[var2], data[var1], alpha=.1)
            #ax.set_xlim((0, 1))
            #ax.set_ylim((0, 1))
            legend = ax.legend()
            ax.set_xlabel(var2, fontsize=16)
            ax.set_ylabel(var1, fontsize=16);
        
        return [models.a[0], models.b[0]],[models.a[models.shape[0]-1], models.b[models.shape[0]-1]]
    
    def AddRule(self,curve,var1,var2,typ):
        if typ=='inf':
            df_temp = pandas.DataFrame([['inf', var2, var1,curve[1],curve[0]]],columns=['type','x', 'y','a','b'])
            df_rule = df_temp
            df_temp = pandas.DataFrame([['sup', var1, var2,1/curve[1],-curve[0]/curve[1]]],columns=['type','x', 'y','a','b'])
            df_rule = pandas.concat([df_rule, df_temp])
        elif typ=='sup':
            df_temp = pandas.DataFrame([['sup', var2, var1,curve[1],curve[0]]],columns=['type','x', 'y','a','b'])
            df_rule = df_temp
            df_temp = pandas.DataFrame([['inf', var1, var2,1/curve[1],-curve[0]/curve[1]]],columns=['type','x', 'y','a','b'])
            df_rule = pandas.concat([df_rule, df_temp])
        return df_rule
        
        
# =============================================================
# Génération de la base SKF modifiée à l'état brut
# =============================================================

D1=DataBaseBearing()
#bearing_SKF=D1.ImproveSKFDataBase()
#bearing_SKF.to_csv(chemin_catalogs+'tableau_rlts_SKF.csv',index=False)

# =============================================================
# Recherche du paramètre B1 de la norme ISO pour matcher les charges dynamiques SKF avec celles de la norme ISO (la norme propose B1=551.13373/0.483 et on obtient au final B1=6.765625*551.13373/0.483 afin que les chages dynamiques corrélent entre la démarche ISO et celle calculé par SKF)
# =============================================================

#bearing_SKF=D1.LoadSKFDataBase()
#coeff_min=6
#coeff_max=8
#def AnalyseCorr(coeff):
#    for ind in (bearing_SKF.name).index:
#        D=bearing_SKF['D'][ind]
#        d=bearing_SKF['d'][ind]
#        B=bearing_SKF['B'][ind]
#        E=bearing_SKF['E'][ind]
#        F=bearing_SKF['F'][ind]
#        Lw=bearing_SKF['Lw'][ind]
#        Dw=bearing_SKF['Dw'][ind]
#        Z=bearing_SKF['Z'][ind]
#        r_roller=bearing_SKF['r_roller'][ind]
#        bearing_temp=bearing.RadialRollerBearing(D=D,d=d,B=B,Lw=Lw,Dw=Dw,r_roller=r_roller,E=E,F=F,Z=Z,i=1,alpha=0,B1=coeff*551.13373/0.483)
#        bearing_temp.BaseStaticLoad()
#        bearing_temp.BaseDynamicLoad()
#        bearing_SKF.loc[ind,'Cr']=bearing_temp.Cr
#    C = bearing_SKF.C
#    Cr = bearing_SKF.Cr
#    Cr=sm.add_constant(Cr)
#    model = sm.OLS(C, Cr).fit()
#    return model.params.Cr
#
#drap=0
#while drap==0:
#    coeff=(coeff_min+coeff_max)/2
#    curve=AnalyseCorr(coeff)
#    if curve>1:
#        coeff_min=coeff
#    else:
#        coeff_max=coeff
#    print(curve,coeff)
#    if abs(curve-1)<0.001:
#        drap=1
#
#bearing_SKF=D1.UpdateSKFDataBase(bearing_SKF)
#coeff=6.765625
#B1=coeff*551.13373/0.483
#D1.AddColumnDataBase('B1',B1)
#bearing_SKF=D1.UpdateSKFDataBase(bearing_SKF)

# =============================================================
# Tracé de l'ensemble des charges statique et dynamique
# =============================================================

#bearing_SKF=D1.LoadSKFDataBase()
#adim1=((bearing_SKF.C-bearing_SKF.Cr)**2).max()
#adim2=((bearing_SKF.C0-bearing_SKF.C0r)**2).max()
#bearing_SKF['fonctionnel']=1/adim1*(bearing_SKF.C-bearing_SKF.Cr)**2+1/adim2*(bearing_SKF.C0-bearing_SKF.C0r)**2
#
#fig, ax = plt.subplots(figsize=(8, 6))
#ax.scatter(bearing_SKF.C, bearing_SKF.Cr,alpha=.1)
#ax.scatter(bearing_SKF.C0, bearing_SKF.C0r,alpha=.2)
#x = npy.arange(bearing_SKF.C.min(), bearing_SKF.C.max(),(bearing_SKF.C.max()-bearing_SKF.C.min())/250)
#get_y = lambda a, b: b + a * x
#y = get_y(1,0)
#C = bearing_SKF.C
#Cr = bearing_SKF.Cr
#Cr=sm.add_constant(Cr)
#model = sm.OLS(C, Cr).fit()
#Cr=get_y(model.params.Cr,model.params.const)
#ax.plot(x, y, linestyle='dotted', color='grey')
#ax.plot(x, Cr, linestyle='dotted', color='red')

# =============================================================
# Ajout de nouvelles colonnes
# =============================================================

#bearing_SKF=D1.LoadSKFDataBase()
#
#ep_ext=(bearing_SKF.D-bearing_SKF.E)/2
#ep_int=(bearing_SKF.F-bearing_SKF.d)/2
#ep=(bearing_SKF.D-bearing_SKF.d)/2
#ep_ax=bearing_SKF.B-bearing_SKF.Lw
#reg1=(ep_ext/ep).values
#reg2=(ep_int/ep).values
#reg3=((ep_ax)/bearing_SKF.B).values
#
#rules_rlts=pandas.DataFrame({'ep_ext':ep_ext,'ep_int':ep_int,
#                             'ep':ep,'ep_ax':ep_ax,'reg_C':bearing_SKF.C/bearing_SKF.Z/bearing_SKF.Lw})
##scatter_matrix(rules_rlts, alpha=0.5, figsize=(6, 6), diagonal='hist')
#D1.AddColumnDataBase(rules_rlts)
#bearing_SKF=D1.UpdateSKFDataBase()

# =============================================================
# Construction des règles métier SKF
# =============================================================

bearing_SKF=D1.LoadSKFDataBase()

#rule 1
var1='ep'
dat1=bearing_SKF.ep
var2='B'
dat2=bearing_SKF.B
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df_rule=D1.AddRule(curve1,var1,var2,'inf')
df=D1.AddRule(curve2,var1,var2,'sup')
df_rule = pandas.concat([df_rule, df])

#rule 2
var1='Lw'
dat1=bearing_SKF.Lw
var2='B'
dat2=bearing_SKF.B
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])
df=D1.AddRule(curve2,var1,var2,'sup')
df_rule = pandas.concat([df_rule, df])

#rule 3
var1='d'
dat1=bearing_SKF.d
var2='B'
dat2=bearing_SKF.B
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])
df=D1.AddRule(curve2,var1,var2,'sup')
df_rule = pandas.concat([df_rule, df])

#rule 4
var1='D'
dat1=bearing_SKF.D
var2='B'
dat2=bearing_SKF.B
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])
df=D1.AddRule(curve2,var1,var2,'sup')
df_rule = pandas.concat([df_rule, df])

#rule 5
var1='E'
dat1=bearing_SKF.E
var2='B'
dat2=bearing_SKF.B
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])
df=D1.AddRule(curve2,var1,var2,'sup')
df_rule = pandas.concat([df_rule, df])

#rule 6
var1='F'
dat1=bearing_SKF.F
var2='B'
dat2=bearing_SKF.B
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])
df=D1.AddRule(curve2,var1,var2,'sup')
df_rule = pandas.concat([df_rule, df])

#rule 7
var1='C0'
dat1=bearing_SKF.C0
var2='B'
dat2=bearing_SKF.B
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])

#rule 7
var1='C0r'
dat1=bearing_SKF.C0r
var2='B'
dat2=bearing_SKF.B
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])

#rule 8
var1='Z'
dat1=bearing_SKF.Z
var2='B'
dat2=bearing_SKF.B
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])

#rule 9
var1='ep_ax'
dat1=bearing_SKF.ep_ax
var2='B'
dat2=bearing_SKF.B
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])
df=D1.AddRule(curve2,var1,var2,'sup')
df_rule = pandas.concat([df_rule, df])

#rule 10
var1='Dw'
dat1=bearing_SKF.Dw
var2='B'
dat2=bearing_SKF.B
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])

#rule 11
var1='D'
dat1=bearing_SKF.D
var2='d'
dat2=bearing_SKF.d
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])

#rule 12
var1='Lw'
dat1=bearing_SKF.Lw
var2='d'
dat2=bearing_SKF.d
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])

#rule 13
var1='Dw'
dat1=bearing_SKF.Dw
var2='d'
dat2=bearing_SKF.d
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)

#rule 14
var1='E'
dat1=bearing_SKF.E
var2='d'
dat2=bearing_SKF.d
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])

#rule 15
var1='F'
dat1=bearing_SKF.F
var2='d'
dat2=bearing_SKF.d
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])

#rule 16
var1='C0'
dat1=bearing_SKF.C0
var2='d'
dat2=bearing_SKF.d
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)

#rule 17
var1='C0r'
dat1=bearing_SKF.C0r
var2='d'
dat2=bearing_SKF.d
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)

#rule 18
var1='Z'
dat1=bearing_SKF.Z
var2='d'
dat2=bearing_SKF.d
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)

#rule 19
var1='Lw'
dat1=bearing_SKF.Lw
var2='D'
dat2=bearing_SKF.D
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])

#rule 20
var1='Dw'
dat1=bearing_SKF.Dw
var2='D'
dat2=bearing_SKF.D
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)

#rule 21
var1='E'
dat1=bearing_SKF.E
var2='D'
dat2=bearing_SKF.D
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])
df=D1.AddRule(curve2,var1,var2,'sup')
df_rule = pandas.concat([df_rule, df])

#rule 22
var1='F'
dat1=bearing_SKF.F
var2='D'
dat2=bearing_SKF.D
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])
df=D1.AddRule(curve2,var1,var2,'sup')
df_rule = pandas.concat([df_rule, df])

#rule 23
var1='C0'
dat1=bearing_SKF.C0
var2='D'
dat2=bearing_SKF.D
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])

#rule 24
var1='C0r'
dat1=bearing_SKF.C0r
var2='D'
dat2=bearing_SKF.D
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])

#rule 25
var1='Dw'
dat1=bearing_SKF.Dw
var2='C0'
dat2=bearing_SKF.C0
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)

#rule 26
var1='Lw'
dat1=bearing_SKF.Lw
var2='C0'
dat2=bearing_SKF.C0
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)

#rule 27
var1='E'
dat1=bearing_SKF.E
var2='C0'
dat2=bearing_SKF.C0
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])
df=D1.AddRule(curve2,var1,var2,'sup')
df_rule = pandas.concat([df_rule, df])

#rule 28
var1='F'
dat1=bearing_SKF.F
var2='C0'
dat2=bearing_SKF.C0
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])
df=D1.AddRule(curve2,var1,var2,'sup')
df_rule = pandas.concat([df_rule, df])

#rule 29
var1='Dw'
dat1=bearing_SKF.Dw
var2='C0r'
dat2=bearing_SKF.C0r
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)

#rule 30
var1='Lw'
dat1=bearing_SKF.Lw
var2='C0r'
dat2=bearing_SKF.C0r
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)

#rule 31
var1='E'
dat1=bearing_SKF.E
var2='C0r'
dat2=bearing_SKF.C0r
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])

#rule 32
var1='F'
dat1=bearing_SKF.F
var2='C0r'
dat2=bearing_SKF.C0r
curve1,curve2=D1.QuantileRegression(var1,var2,pandas.DataFrame({var1:dat1,var2:dat2}),0)
df=D1.AddRule(curve1,var1,var2,'inf')
df_rule = pandas.concat([df_rule, df])
df=D1.AddRule(curve2,var1,var2,'sup')
df_rule = pandas.concat([df_rule, df])

df_rule.to_csv(chemin_catalogs+'rules_rlts_SKF.csv',index=False)