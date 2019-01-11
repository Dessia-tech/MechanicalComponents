#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 02:14:21 2018

@author: steven
"""

from mechanical_components.bearings import oil_iso_vg_1500, material_iso, dico_rlts_iso,dico_rules,dico_roller_iso
from mechanical_components.bearings import RadialBallBearing, AngularBallBearing, \
        SphericalBallBearing, RadialRollerBearing, TaperedRollerBearing, \
        BearingAssembly, BearingCombination, DetailedRadialRollerBearing, \
        BearingAssemblyOptimizationResults, LoadBearingAssemblyResults
import numpy as npy
from scipy.optimize import minimize
import pandas
#from pandas.plotting import scatter_matrix
import pkg_resources
from itertools import product

#import genmechanics
#import genmechanics.linkages as linkages
#import genmechanics.loads as loads
import math
import volmdlr as vm
#import volmdlr.primitives3D as primitives3D
#import volmdlr.primitives2D as primitives2D

# TODO: Cumulative damages
class RollerBearingContinuousOptimizer:

    """
    Objet avec 3 fonctions de selection des roulements cylindriques
     * Combinatoire sur les dimensions externe ISO
     * Combinatoire en prenant en compte les règles SKF
     * Estimation des durées de vie et charge dynamique et fonction de tri

    """
    def __init__(self):
        
        self.bearing_assemblies=[]
    
    def Optimization(self, d, D, B, Fr, Fa, N, t, T, L10=None, C0r=None, Cr=None,
                         Lnm=None, grade=['Gr_gn'], S=0.9,
                         oil=oil_iso_vg_1500, material=material_iso,
                         number_solutions=1, maxi=None, mini=None, rsmin=None, typ='NF',
                         verbose = False):

#        self.d = d
#        self.D = D
#        self.B = B
#        self.Fr = Fr
#        self.Fa = Fa
#        self.
        
        err_default=0.05
        def def_inter(data):
            if data==None:
                # TODO utiliser math.inf plutot
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
            Dw_sup=min(Dw_sup,(D-d)/2.)
            for Dw in dico_roller_iso.keys():
                if (Dw>=Dw_inf) and (Dw<=Dw_sup):
                    for Lw in dico_roller_iso[Dw].keys():
                        if (Lw>=Lw_inf) and (Lw<=Lw_sup):
                            E_inf,E_sup=analyse_rule('E',-npy.inf,npy.inf,['d','D','B','Dw','Lw'],[d,D,B,Dw,Lw])
                            F_inf,F_sup=analyse_rule('F',-npy.inf,npy.inf,['d','D','B','Dw','Lw'],[d,D,B,Dw,Lw])
                            if ((E_inf-F_sup)/2.<=Dw) and ((E_sup-F_inf)/2.>=Dw):
                                liste_sol_roller_iso.append([d,D,B,rsmin,serial,Dw,Lw,E_inf,E_sup,F_inf,F_sup])


        #Analyse de détail des roulements
        a_ED1_inf,b_ED1_inf=dico_rules[('E','D1','inf')]
        a_ED1_sup,b_ED1_sup=dico_rules[('E','D1','sup')]
        a_Fd1_inf,b_Fd1_inf=dico_rules[('F','d1','inf')]
        a_Fd1_sup,b_Fd1_sup=dico_rules[('F','d1','sup')]
        estim_masse=[]
        #Construction d'une fonctionnelle pour le tri afin d'engager une optimisation continue
        for [d,D,B,rsmin,serial,Dw,Lw,E_inf,E_sup,F_inf,F_sup] in liste_sol_roller_iso:
            E=(E_inf+E_sup)/2.
            F=E-2*Dw
            Zmax=int(2*npy.pi/(2.*math.arcsin((Dw/2)/(F/2.+Dw/2.))))
            masse_elem=(npy.pi*D**2/4.-npy.pi*E**2/4.)*B
            masse_elem+=(npy.pi*F**2/4.-npy.pi*d**2/4.)*B
            masse_elem+=npy.pi*Dw**2/4.*Lw*Zmax
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
            R1 = DetailedRadialRollerBearing(d, D, B, i = 1, Z = 0, Lw = Lw, Dw = Dw, radius = rsmin)
            def fun(x):
                obj=0
                F=x[0]-2*Dw
                Zmax=int(2*math.pi/(2*npy.arcsin((Dw/2.)/(F/2.+Dw/2.))))
                R1.Update(x[2],x[1],x[0],F,Zmax)
                l10=R1.BaseLifeTime(Fr, Fa, N, t, Cr = R1.BaseDynamicLoad())

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
                ine.append((x[0]+F)/2.-x[2])  #d1 inférieur au diamètre passant par l'axe des rouleaux
                return ine
            cons = ({'type': 'ineq','fun' : ineg})
            valid_optim=1
            while valid_optim==1:
                Bound=[[E_inf,E_sup],[d,D],[d,D]]
                x0=(npy.array(Bound)[:,1]-npy.array(Bound)[:,0])*npy.random.random(3)+npy.array(Bound)[:,0]
                res = minimize(fun,x0, method='SLSQP', bounds=Bound,constraints=cons)
                if (min(ineg(res.x))>=0):
                    valid_optim=0
                    
            #Validation finale vis à vis des bornes du CDC
            if valid_optim==0:
                x_opt=res.x
                F=x_opt[0]-2*Dw
                Zmax=int(2*math.pi/(2*npy.arcsin((Dw/2.)/(F/2.+Dw/2.))))
                R1.Update(x_opt[2],x_opt[1],x_opt[0],F,Zmax)
                l10=R1.BaseLifeTime(Fr, Fa, N, t, Cr = R1.BaseDynamicLoad())

                
                liminf_lifetime=True
                if L10 is not None:
                    if (l10<L10['min']) or (l10>L10['max']):
                        liminf_lifetime=False
                if Lnm is not None:
                    lnm=R1.AdjustedLifeTime(Fr, Fa, N, t, T, R1.BaseDynamicLoad(), R1.BaseStaticLoad(), S)

                    if (lnm<Lnm['min']) or (lnm>Lnm['max']):
                        liminf_lifetime=False
                if C0r is not None:
                    c0r=R1.BaseStaticLoad()
                    if (c0r<C0r['min']) or (c0r>C0r['max']):
                        liminf_lifetime=False
                if Cr is not None:
                    cr=R1.BaseDynamicLoad()
                    if (cr<Cr['min']) or (cr>Cr['max']):
                        liminf_lifetime=False
                if liminf_lifetime==True:
                    self.bearing_assemblies.append(R1)
                    if len(self.bearing_assemblies)<number_solutions:
                        liminf_lifetime=False
                    if verbose:
#                        print('Bearing solution n°{} with L10:{},D:{},d:{},B:{},Dw:{},Z:{}'.format(len(self.bearing_assemblies),l10,D,d,B,Dw,Zmax))
                        if L10 is not None:
                            print('L10: {}, specification: {}'.format(l10, L10))
                        if Lnm is not None:

                            print('Lnm: {}, specification: {}'.format(lnm, Lnm))
                            
with pkg_resources.resource_stream(pkg_resources.Requirement('mechanical_components'),'mechanical_components/catalogs/SNR.csv') as rlts_FNR:
    pandas_rlts_FNR = pandas.read_csv(rlts_FNR) 
            
pandas_sort = pandas_rlts_FNR[pandas_rlts_FNR['i'].notnull()]
pandas_sort = pandas_sort[pandas_sort['Z'].notnull()]
pandas_sort = pandas_sort[pandas_sort['Dw'].notnull()]
pandas_sort = pandas_sort[pandas_sort['mass'].notnull()]
#        pandas_sort = pandas_sort[pandas_sort['alpha'].notnull()]
pandas_sort = pandas_sort[pandas_sort['Cr'].notnull()]
pandas_sort = pandas_sort[pandas_sort['C0r'].notnull()]
base_bearing = pandas_sort

class DiscreteBearingCombinationOptimizer:
    
    def __init__(self, linkage, behavior_link, nb_rlts, d, D, length, number_solutions=10, 
                 number_load_case=1):
        
        axis_load = ['n', 'p'] # 'p' when the axial load applied on shaft is directed from left to right
        typ_rlts = {'rol_NU':0, 'rol_N':0, 
                    'rol_NJ_p':'p', 'rol_NJ_n':'n',
                    'rol_NF_p':'p', 'rol_NF_n':'n', 
                    'rol_NUP':'pn', 'ball':'pn', 
                    'ang_p':'p', 'ang_n':'n', 'tap_p':'p', 'tap_n':'n'}
        ring_mounting = ['n', 'p']
        
        self.bearing_combinations = []
        for li_rlts in product(typ_rlts.keys(), repeat = nb_rlts):
            if (linkage is 'ball_joint'):
                if (nb_rlts > 1):
                    continue
                elif li_rlts[0] not in ['ball', 'ang_p', 'ang_n']:
                    continue
            else:
                if (nb_rlts is 1) and (li_rlts[0] in ['ball', 'ang_p', 'ang_n']):
                    continue
            li_rlts_temp = []
            for bg in li_rlts:
                li_rlts_temp.append(typ_rlts[bg])
            if (('n' in li_rlts_temp) and ('p' in li_rlts_temp)) or ('pn' in li_rlts_temp):
                behavior_assembly = 'pn'
            elif 'p' in li_rlts_temp:
                behavior_assembly = 'p'
            elif 'n' in li_rlts_temp:
                behavior_assembly = 'n'
            else:
                behavior_assembly = 0
                
            if (behavior_link is 'pn') and (behavior_assembly is not 'pn'):
                continue
            elif (behavior_link is 'p') and (behavior_assembly in [0, 'n']):
                continue
            elif (behavior_link is 'n') and (behavior_assembly in [0, 'p']):
                continue
            elif (behavior_link is 0) and (behavior_assembly in ['p', 'n']):
                continue
            
            # selection combinatory
            # if one angular we don't want a tapered and so on
            if (('ang_p' in li_rlts) or ('ang_n' in li_rlts)) and (('tap_p' in li_rlts) or ('tap_n' in li_rlts)):
                continue
            #if 'n' and 'p' so we want two tappered or two angular
            if ('p' in li_rlts_temp) and ('n' in li_rlts_temp): 
                if (('tap_p' in li_rlts) and ('tap_n' in li_rlts)):
                    pass
                elif (('ang_p' in li_rlts) and ('ang_n' in li_rlts)):
                    pass
                else:
                    continue
                
            if behavior_link in ['pn', 'p']:
                if typ_rlts[li_rlts[0]] in ['n', 0]:
                    continue
                if typ_rlts[li_rlts[-1]] in ['n', 0]:
                    continue
                loc_axial_load = []
                ind_bg_axial = 0
                for num_bg, bg in enumerate(li_rlts[0:-1]):
                    if bg == li_rlts[num_bg + 1]:
                        ind_bg_axial += 1
                    else:
                        loc_axial_load.append(ind_bg_axial)
                        ind_bg_axial = 0
            
            if behavior_link == 'pn':
                mount = {'be':['p','n'], 'bi':['p','n']}
            elif behavior_link == 'p':
                mount = {'be':['p'], 'bi':['n']}
            elif behavior_link is 'n':
                mount = {'be':['n'], 'bi':['p']}
            elif behavior_link == 0:
                if behavior_assembly == 0:
                    mount = {'be':['p','n'], 'bi':['p','n']}
                elif behavior_assembly == 'pn':
                    mount = {'be':[], 'bi':['p','n']}
            else:
                print('Unknown behavior_link: ', behavior_link)
                raise NotImplementedError
            li_rlts_pandas = self.SearchCatalog(d, D, length, li_rlts, number_solutions = number_solutions)
            if li_rlts_pandas is not None:
                for li_bearing in li_rlts_pandas:
                    list_bearing = []
                    for index_bearing, code_bearing in zip(li_bearing, li_rlts):
                        list_bearing.append(self.GenerateBearing(index_bearing, code_bearing))
                    radial_load_linkage = [True]*len(list_bearing)
                    BA = BearingCombination(list_bearing, radial_load_linkage, connection_bi = mount['bi'], 
                                     connection_be = mount['be'], behavior_link = behavior_link,
                                     number_load_case = number_load_case)
                    if BA.check:
                        self.bearing_combinations.append(BA)
        
            
    def GenerateBearing(self, index, code_bearing):
        typ_rlt = base_bearing.loc[index,'typ_bearing']
        d = base_bearing.loc[index,'d']
        D = base_bearing.loc[index,'D']
        B = base_bearing.loc[index,'B']
        i = base_bearing.loc[index,'i']
        Z = base_bearing.loc[index,'Z']
        Dw = base_bearing.loc[index,'Dw']
        mass = base_bearing.loc[index,'mass']
        alp = base_bearing.loc[index,'alpha']
        if str(alp) == 'nan':
            alpha = 0
        else:
            alpha = alp
        Cr = base_bearing.loc[index,'Cr']
        C0r = base_bearing.loc[index,'C0r']
        if typ_rlt == 'radial_roller_bearing':
            if len(code_bearing.split('_')) == 3:
                typ = code_bearing.split('_')[1]
                direction = (code_bearing.split('_')[2] == 'n' and -1 or 1)
            elif len(code_bearing.split('_')) == 2:
                typ = code_bearing.split('_')[1]
                direction = 1
            return RadialRollerBearing(d, D, B, i, Z, Dw, alpha, Cr, C0r, 
                                              mass = mass, typ = typ, direction = direction)
        elif typ_rlt == 'radial_ball_bearing':
            return RadialBallBearing(d, D, B, i, Z, Dw, alpha, Cr, C0r, mass = mass)
        elif typ_rlt == 'angular_ball_bearing':
            direction = (code_bearing.split('_')[1] == 'n' and 1 or -1)
            return AngularBallBearing(d, D, B, i, Z, Dw, alpha,
                                             Cr, C0r, direction = direction, mass = mass)
        elif typ_rlt == 'spherical_ball_bearing':
            return SphericalBallBearing(d, D, B, i, Z, Dw, alpha, 
                                               Cr, C0r, mass = mass)
        elif typ_rlt == 'tapered_roller_bearing':
            direction = (code_bearing.split('_')[1] == 'n' and 1 or -1)
            return TaperedRollerBearing(d, D, B, i, Z, Dw, alpha, 
                                               Cr, C0r, mass = mass, direction = direction)
    
    def SearchCatalog(self, d, D, length, list_bearing=None, 
               constraint=[{'typ':'equal','var':'d'}, {'typ':'equal','var':'D'}],
               number_solutions=5,
               verbose=False):
        
        valid = True
        
        tab_rlts = base_bearing[d[1] >= base_bearing.d]
        tab_rlts = tab_rlts[d[0] <= tab_rlts.d]
        tab_rlts = tab_rlts[D[0] <= tab_rlts.D]
        tab_rlts = tab_rlts[D[1] >= tab_rlts.D]
        
        list_typ_bearing = []
        list_name_bearing = []
        list_index = []
        for bearing in list_bearing:
            name_bearing = bearing.split('_')[0]
            typ_bearing = None
            if name_bearing == 'rol':
                typ_bearing = bearing.split('_')[1]
                name_bearing = 'radial_roller_bearing'
            if name_bearing == 'ball':
                name_bearing = 'radial_ball_bearing'
            if name_bearing == 'ang':
                name_bearing = 'angular_ball_bearing'
            if name_bearing == 'tap':
                name_bearing = 'tapered_roller_bearing'
            list_name_bearing.append(name_bearing)
            list_typ_bearing.append(typ_bearing)
            if typ_bearing == None:
                tab_rlts_temp = tab_rlts[name_bearing == tab_rlts.typ_bearing]
            else:
                tab_rlts_temp = tab_rlts[(name_bearing == tab_rlts.typ_bearing) & (typ_bearing == tab_rlts.mounting)]
            if len(list(tab_rlts_temp.index)) == 0:
                valid = False
            list_index.extend(list(tab_rlts_temp.index))
        if valid:
            list_index = list(set(list_index))
    #        tab_rlts = tab_rlts.loc[tab_rlts['typ_bearing'].isin(list_typ_bearing)]
            tab_rlts = tab_rlts.loc[list_index, :]
            solu = pandas.crosstab([tab_rlts['typ_bearing']],[tab_rlts['d'],tab_rlts['D']])
            
            # crosstab approach is not perform on all type of roller_bearing but it's not a problem: differents typ (NU,N,NF,NUP) have same d, D and B
            crosstab_data = solu.values
            output_data = []
            for ind_column in range(crosstab_data.shape[1]):
                if 0 not in crosstab_data[:, ind_column]:
                    d = solu.columns.levels[0][solu.columns.labels[0][ind_column]]
                    D = solu.columns.levels[1][solu.columns.labels[1][ind_column]]
                    output_data.append({'d':d, 'D':D})
                    
            pd_bearing = []
            for output in output_data[0:5]:
                pandas_iter = tab_rlts[(tab_rlts.d == output['d']) & (tab_rlts.D == output['D'])]
                optim_a_bearing = []
                optim_b_bearing = []
                for typ_bear, name_bear in zip(list_typ_bearing,list_name_bearing):
                    if name_bear == 'radial_roller_bearing':
                        pandas_iter1 = pandas_iter[(pandas_iter.typ_bearing == name_bear) & (pandas_iter.mounting == typ_bear)]
                    else:
                        pandas_iter1 = pandas_iter[pandas_iter.typ_bearing == name_bear]
                    if not pandas_iter1.empty:
                        pandas_iter1 = pandas_iter1.sort_values(by = ['B']).loc[:,['d','D','B','Cr','C0r','Z','typ_bearing']]
                        optim_a_bearing.append(list(pandas_iter1.index)[0])
                        pandas_iter1 = pandas_iter1.sort_values(by = ['Cr'], ascending=False)
                        optim_b_bearing.append(list(pandas_iter1.index)[0])
                    
                for list_optim in [optim_a_bearing, optim_b_bearing]:
                    length_iter = 0
                    for index_bearing in list_optim:
                        li = tab_rlts.loc[index_bearing,'B']
                        length_iter += li
                    if (length_iter <= length[1]) and (len(list_optim) == len(list_bearing)):
                        pd_bearing.append(list_optim)
            if number_solutions is not None:
                return pd_bearing[0:min(number_solutions, len(pd_bearing))]
            else:
                return pd_bearing[0:]
        
class SortBearingCombinationOptimizer:
    def __init__(self, bearing_combinations, sort_arg={'min':'mass'}, number_solutions=10):
        
        bearing_combinations_sort = self.Sortarchitectures(bearing_combinations, sort_arg = {'min':'mass'})
        if number_solutions is not None:
            self.bearing_combinations = bearing_combinations_sort[0:min(len(bearing_combinations_sort), number_solutions)]
        else:
            self.bearing_combinations = bearing_combinations_sort[0:]
        self.check = True
        if self.bearing_combinations == []:
            self.check = False
            
    def Sortarchitectures(self, architectures, sort_arg):
        list_sort = []
        for sol in architectures:
            val = getattr(sol,sort_arg['min'])
            list_sort.append(val)
        list_sort_arg=list(npy.argsort(list_sort))
        return npy.array(architectures)[list_sort_arg]
    
class DiscreteBearingAssemblyOptimizer:
    def __init__(self, loads, speeds, operating_times,
                 inner_diameter=[0.02, 0.02],
                 axial_positions=[0, 0.1],
                 outer_diameter=[0.05, 0.05],
                 length=[0.04, 0.04],
                 linkage_types=[['all'], ['all']],
                 mounting_types=None,
                 number_bearings=[[1, 2], [1, 2]]):
        
        self.loads = loads
        self.speeds = speeds
        self.operating_times = operating_times
        self.inner_diameter = inner_diameter
        self.axial_positions = axial_positions
        self.outer_diameter = outer_diameter
        self.length = length
        self.linkage_types = linkage_types
        self.number_bearings = number_bearings
        
        nb_linkage = len(length)
        if mounting_types == None:
            mounting_types = []
            for typ_iter in product(['pn', 0, 'n', 'p'], repeat = nb_linkage):
                if typ_iter is not (0,)*nb_linkage:
                    if set(typ_iter) not in [set(('pn', 'p')), set(('pn', 'n')), set(('pn', 'pn')), set((0, 0))]:
                        mounting_types.append(typ_iter)
        self.mounting_types = mounting_types
        for mount in self.mounting_types:
            if set(mount) in [set(('pn', 'p')), set(('pn', 'n')), set(('pn', 'pn')), set((0, 0))]:
                raise KeyError('Modify the mounting_types data')
                
        self.bearing_combination_configurations = self.CombinationOptimize()

    def CombinationOptimize(self):
        list_linkage = self.linkage_types 
        for num_link,list_link in enumerate(self.linkage_types):
            if 'all' in list_link:
                list_linkage[num_link] = ['ball_joint', 'cylindric_joint']

        bearing_combination_configurations = []
        for li_mounting in self.mounting_types:
            for num_linkage, behavior_link in enumerate(li_mounting):
                for linkage, nb_rlts in product(list_linkage[num_linkage], 
                                                self.number_bearings[num_linkage]):
                    bearing_combination_configurations.append([num_linkage, linkage, behavior_link, nb_rlts])
        return bearing_combination_configurations
                    
    @classmethod
    def DefOptimizer(cls, loads, speeds, times,
                 inner_diameter, axial_positions, outer_diameter, length,
                 linkage_types, mounting_types, number_bearing,
                 sort, sort_arg, path, number_solutions):
        obj = cls(loads, speeds, times,
                 inner_diameter, axial_positions, outer_diameter, length,
                 linkage_types, mounting_types, number_bearing,
                 sort, sort_arg, path, number_solutions)
        return obj
        
class SortBearingAssemblyOptimizer:
    def __init__(self, bearing_assemblies, sort_arg = {'min':'mass'},
                 number_solutions=20):
        sort_bearing_assemblies = self.Sortarchitectures(bearing_assemblies, sort_arg)
        if number_solutions is not None:
            self.bearing_assemblies = sort_bearing_assemblies[0: min(len(sort_bearing_assemblies), number_solutions)]
        else:
            self.bearing_assemblies = sort_bearing_assemblies[0:]
        
    def Sortarchitectures(self, bearing_assemblies, sort_arg):
        list_sort = []
        for bearing_assembly in bearing_assemblies:
            val = getattr(bearing_assembly,sort_arg['min'])
            list_sort.append(val)
        list_sort_arg=list(npy.argsort(list_sort))
        return npy.array(bearing_assemblies)[list_sort_arg]
    
class ContinuousBearingAssemblyOptimizer:
    def __init__(self, bearing_assemblies, loads, speeds, operating_times,
                 inner_diameter=[0.02, 0.02],
                 axial_positions=[0, 0.1],
                 outer_diameter=[0.05, 0.05],
                 length=[0.04, 0.04],
                 number_solutions=5):
        
        self.inner_diameter = inner_diameter
        self.axial_positions = axial_positions
        self.outer_diameter = outer_diameter
        self.length = length
        load_bearing_assembly_results = LoadBearingAssemblyResults(loads, speeds,
                                                             operating_times)
        self.bearing_assemblies = self.Optimize(bearing_assemblies,
                                                load_bearing_assembly_results,
                                                number_solutions)
        
    def Optimize(self, bearing_assemblies, load_bearing_assembly_results,
                 number_solutions=5, Lnm_min=1e4, 
                 index_sol=None, verbose=False):
        
        if index_sol is not None:
            list_sol = []
            for num_sol, sol in enumerate(bearing_assemblies):
                if num_sol in index_sol:
                    list_sol.append(sol)
        else:
            list_sol = bearing_assemblies
        bearing_assemblies = []
        for bearing_assembly in list_sol:
            l1 = bearing_assembly.bearing_combinations[0].B
            l2 = bearing_assembly.bearing_combinations[1].B
            pos1_min = self.axial_positions[0] + l1/2.
            pos1_max = self.axial_positions[0] + self.length[0] - l1/2.
            pos2_min = self.axial_positions[1] + l2/2.
            pos2_max = self.axial_positions[1] + self.length[1] - l2/2.
            def fun(x):
                bearing_assembly.ShaftLoad([x[0], x[1]], load_bearing_assembly_results)
                ld1 = [bc.load_bearing_combination_results.radial_loads for bc in bearing_assembly.bearing_combinations]
                print(ld1)
                obj = 0
                for (l1) in ld1:
                    obj += 1/(l1[0])**2
#                Lnm1 = 0
#                for rlt1 in list_bearing1:
#                    Lnm1 += rlt1.AdjustedLifeTime(Fr = [fr1], Fa = [fa1], N = [300], t = [1e10], T = [60])
#                Lnm2 = 0
#                for rlt2 in list_bearing2:
#                    Lnm2 += rlt2.AdjustedLifeTime(Fr = [fr2], Fa = [fa2], N = [300], t = [1e10], T = [60])
                return obj
            def fineq(x):
                ineq = [0]
#                fa1, fr1, fa2, fr2 = bearing_assembly.ShaftLoad([x[0], x[1]], self.list_pos_unknown, 
#                                                    self.list_load, self.list_torque)
#                Lnm1 = 0
#                for rlt1 in list_bearing1:
#                    Lnm1 += rlt1.AdjustedLifeTime(Fr = [fr1], Fa = [fa1], N = [300], t = [1e10], T = [60])
#                Lnm2 = 0
#                for rlt2 in list_bearing2:
#                    Lnm2 += rlt2.AdjustedLifeTime(Fr = [fr2], Fa = [fa2], N = [300], t = [1e10], T = [60])
#                ineq.append(Lnm1 - Lnm_min)
#                ineq.append(Lnm2 - Lnm_min)
                return ineq
            Bound = [[pos1_min, pos1_max], [pos2_min, pos2_max]]
            sol_fun = npy.inf
            for p1, p2 in product(Bound[0],Bound[1]):
                cons = {'type': 'ineq','fun' : fineq}
                res = minimize(fun, [p1, p2], method='SLSQP', bounds=Bound, constraints = cons)
                if fun(res.x) < sol_fun:
                    sol_fun = fun(res.x)
                    sol_x = res.x
                    status = res.status
            for itera in range(0,5):
                x0 = (npy.array(Bound)[:,1]-npy.array(Bound)[:,0])*npy.random.random(2)+npy.array(Bound)[:,0]
                cons = {'type': 'ineq','fun' : fineq}
                res = minimize(fun, x0, method='SLSQP', bounds=Bound, constraints = cons)
                if fun(res.x) < sol_fun:
                    sol_fun = fun(res.x)
                    sol_x = res.x
                    status = res.status
            
            if number_solutions is not None:
                if len(bearing_assemblies) <= number_solutions:       
                    if status >= 0:
                        bearing_assembly.Update(sol_x, self.inner_diameter, self.axial_positions, 
                                            self.outer_diameter, self.length)
                        
                        if bearing_assembly not in bearing_assemblies:
                            bearing_assemblies.append(bearing_assembly)
                            
                else:
                    break
            else:
                if status >= 0:
                    bearing_assembly.Update(sol_x, self.inner_diameter, self.axial_positions, 
                                        self.outer_diameter, self.length)
                    
                    if bearing_assembly not in bearing_assemblies:
                        bearing_assemblies.append(bearing_assembly)
        return bearing_assemblies
    
class BearingAssemblyOptimizer:
    def __init__(self, loads, speeds, operating_times,
                 inner_diameter=[0.02, 0.02],
                 axial_positions=[0, 0.1],
                 outer_diameter=[0.05, 0.05],
                 length=[0.04, 0.04],
                 linkage_types=[['all'], ['all']],
                 mounting_types=None,
                 number_bearings=[[1, 2], [1, 2]],
                 sort={'typ':'Lnm', 'min':1e4, 'max':1e10},
                 sort_arg = {'min':'mass'},
                 number_solutions=[5, 20, 10, 10]):
        self.loads = loads
        self.speeds = speeds
        self.operating_times = operating_times
        self.inner_diameter = inner_diameter
        self.axial_positions = axial_positions
        self.outer_diameter = outer_diameter
        self.length = length
        self.linkage_types = linkage_types
        self.number_bearings = number_bearings
        self.sort = sort
        self.sort_arg = sort_arg
        self.number_solutions = number_solutions
        
        DBA = DiscreteBearingAssemblyOptimizer(loads, speeds, operating_times,
                                         inner_diameter=inner_diameter,
                                         axial_positions=axial_positions,
                                         outer_diameter=outer_diameter,
                                         length=length,
                                         linkage_types=linkage_types,
                                         mounting_types=mounting_types,
                                         number_bearings=number_bearings)
        
        sol_BC = {}
        for num_linkage, linkage, behavior_link, nb_rlts in DBA.bearing_combination_configurations:
            if num_linkage not in sol_BC.keys():
                sol_BC[num_linkage] = []
                
            DBC = DiscreteBearingCombinationOptimizer(linkage, 
                                                      behavior_link, 
                                                      nb_rlts, 
                                                      [self.inner_diameter[num_linkage], self.outer_diameter[num_linkage]], 
                                                      [self.inner_diameter[num_linkage], self.outer_diameter[num_linkage]], 
                                                      [0, self.length[num_linkage]], 
                                                      number_solutions=self.number_solutions[3], 
                                                      number_load_case=len(self.loads))
            
            SBC = SortBearingCombinationOptimizer(DBC.bearing_combinations, 
                                                  sort_arg={'min':'mass'}, 
                                                  number_solutions=self.number_solutions[2])
            
            if SBC.check is True:
                sol_BC[num_linkage].extend(SBC.bearing_combinations)
        analyze_architectures = True
        for num_linkage, bearing_combinations in sol_BC.items():
            if sol_BC[num_linkage] == []:
                analyze_architectures = False
        bearing_assemblies = []
        if analyze_architectures:
            for composite_solution in product(*sol_BC.values()):
                bearing_assemblies.append(BearingAssembly(composite_solution))
                    
        SBA = SortBearingAssemblyOptimizer(bearing_assemblies, sort_arg,
                                           number_solutions[1])
        
        CBA = ContinuousBearingAssemblyOptimizer(SBA.bearing_assemblies, loads, speeds, 
                                                 operating_times,
                                                 inner_diameter=inner_diameter,
                                                 axial_positions=axial_positions,
                                                 outer_diameter=outer_diameter,
                                                 length=length,
                                                 number_solutions=number_solutions[0])
        
        self.bearing_assemblies = CBA.bearing_assemblies
