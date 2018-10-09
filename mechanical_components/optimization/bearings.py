#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 02:14:21 2018

@author: steven
"""

from mechanical_components.bearings import oil_iso_vg_1500, material_iso, dico_rlts_iso,dico_rules,dico_roller_iso
from mechanical_components.bearings import ConceptRadialBallBearing, ConceptAngularBallBearing, \
        ConceptSphericalBallBearing, ConceptRadialRollerBearing, ConceptTaperedRollerBearing, \
        CompositiveBearingAssembly, BearingAssembly
import numpy as npy
from scipy.optimize import minimize
import pandas
from pandas.plotting import scatter_matrix
import pkg_resources
from itertools import product

import genmechanics
import genmechanics.linkages as linkages
import genmechanics.loads as loads

import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D

# TODO: Cumulative damages
class BearingContinuousOptimizer:
    """
    Objet avec 3 fonctions de selection des roulements cylindriques
     * Combinatoire sur les dimensions externe ISO
     * Combinatoire en prenant en compte les règles SKF
     * Estimation des durées de vie et charge dynamique et fonction de tri

    """
    def __init__(self):
        
        self.solutions=[]
    
    def Optimization(self, d, D, B, Fr, Fa, N, t, T, L10=None, C0r=None, Cr=None,
                         Lnm=None, grade=['Gr_gn'], S=0.9,
                         oil=oil_iso_vg_1500, material=material_iso,
                         nb_sol=1, maxi=None, mini=None, rsmin=None, typ='NF',
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
                    self.solutions.append(R1)
                    if len(self.solutions)<nb_sol:
                        liminf_lifetime=False
                    if verbose:
#                        print('Bearing solution n°{} with L10:{},D:{},d:{},B:{},Dw:{},Z:{}'.format(len(self.solutions),l10,D,d,B,Dw,Zmax))
                        if L10 is not None:
                            print('L10: {}, specification: {}'.format(l10, L10))
                        if Lnm is not None:
                            print('Lnm: {}, specification: {}'.format(lnm, Lnm))
                            
with pkg_resources.resource_stream(pkg_resources.Requirement('mechanical_components'),'mechanical_components/catalogs/tableau_rlts_SNR.csv') as rlts_FNR:
    pandas_rlts_FNR = pandas.read_csv(rlts_FNR) 
            
pandas_sort = pandas_rlts_FNR[pandas_rlts_FNR['i'].notnull()]
pandas_sort = pandas_sort[pandas_sort['Z'].notnull()]
pandas_sort = pandas_sort[pandas_sort['Dw'].notnull()]
pandas_sort = pandas_sort[pandas_sort['mass'].notnull()]
#        pandas_sort = pandas_sort[pandas_sort['alpha'].notnull()]
pandas_sort = pandas_sort[pandas_sort['Cr'].notnull()]
pandas_sort = pandas_sort[pandas_sort['C0r'].notnull()]
base_bearing = pandas_sort


class BearingAssemblyOptimizer:
    
    def __init__(self, linkage, behavior_link, nb_rlts, d, D, length, nb_sol=[10, 10], 
                 sort_arg = {'min':'mass'}):
        
        axis_load = ['n', 'p'] # 'p' when the axial load applied on shaft is directed from left to right
        typ_rlts = {'rol_NU':0, 'rol_N':0, 
                    'rol_NJ_p':'p', 'rol_NJ_n':'n',
                    'rol_NF_p':'p', 'rol_NF_n':'n', 
                    'rol_NUP':'pn', 'ball':'pn', 
                    'ang_p':'p', 'ang_n':'n', 'tap_p':'p', 'tap_n':'n'}
        ring_mounting = ['n', 'p']
        
        solutions = []
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
            if (('ang_p' in li_rlts) or ('ang_n' in li_rlts)) and \
                        (('tap_p' in li_rlts) or ('tap_n' in li_rlts)):
                continue

            if behavior_link is 'pn':
                mount = {'be':['p','n'], 'bi':['p','n']}
            elif behavior_link is 'p':
                mount = {'be':['p'], 'bi':['n']}
            elif behavior_link is 'n':
                mount = {'be':['n'], 'bi':['p']}
            elif behavior_link is 0:
                if behavior_assembly is 0:
                    mount = {'be':['p','n'], 'bi':['p','n']}
                elif behavior_assembly is 'pn':
                    mount = {'be':[], 'bi':['p','n']}
            li_rlts_pandas = self.SearchCatalog(d, D, length, li_rlts, nb_sol = nb_sol[1])
            if li_rlts_pandas is not None:
                for li_bearing in li_rlts_pandas:
                    list_bearing = []
                    for index_bearing, code_bearing in zip(li_bearing, li_rlts):
                        list_bearing.append(self.GenereBearing(index_bearing, code_bearing))
                    BA = BearingAssembly(list_bearing, connection_bi = mount['bi'], 
                                     connection_be = mount['be'])
                    solutions.append(BA)
        solutions_sort = self.SortSolutions(solutions, sort_arg = {'min':'mass'})
        self.solutions = solutions_sort[0:min(len(solutions_sort), nb_sol[0])]
        self.check = True
        if self.solutions == []:
            self.check = False
            
    def GenereBearing(self, index, code_bearing):
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
            return ConceptRadialRollerBearing(d, D, B, i, Z, Dw, alpha, Cr, C0r, 
                                              mass = mass, typ = typ, direction = direction)
        elif typ_rlt == 'radial_ball_bearing':
            return ConceptRadialBallBearing(d, D, B, i, Z, Dw, alpha, Cr, C0r, mass = mass)
        elif typ_rlt == 'angular_ball_bearing':
            direction = (code_bearing.split('_')[1] == 'n' and 1 or -1)
            return ConceptAngularBallBearing(d, D, B, i, Z, Dw, alpha,
                                             Cr, C0r, direction = direction, mass = mass)
        elif typ_rlt == 'spherical_ball_bearing':
            return ConceptSphericalBallBearing(d, D, B, i, Z, Dw, alpha, 
                                               Cr, C0r, mass = mass)
        elif typ_rlt == 'tapered_roller_bearing':
            direction = (code_bearing.split('_')[1] == 'n' and 1 or -1)
            return ConceptTaperedRollerBearing(d, D, B, i, Z, Dw, alpha, 
                                               Cr, C0r, mass = mass, direction = direction)
    
    def SearchCatalog(self, d, D, length, list_bearing=None, 
               constraint=[{'typ':'equal','var':'d'}, {'typ':'equal','var':'D'}],
               nb_sol=5,
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
                    if length_iter <= length[1]:
                        pd_bearing.append(list_optim)
            return pd_bearing[0:min(nb_sol, len(pd_bearing))]
        
    def SortSolutions(self, solutions, sort_arg):
        list_sort = []
        for sol in solutions:
            val = getattr(sol,sort_arg['min'])
            list_sort.append(val)
        list_sort_arg=list(npy.argsort(list_sort))
        return npy.array(solutions)[list_sort_arg]
        

class CompositeBearingAssemblyOptimizer:
    def __init__(self, list_pos_unknown, list_load, list_torque, list_speed, list_time,
                 d_shaft_min=0.02, axial_pos=[0, 0.1], d_ext=[0.05, 0.05], length=[0.04, 0.04],
                 typ_linkage=[['all'], ['all']],
                 typ_mounting=None, number_bearing=[[1, 2], [1, 2]],
                 sort={'typ':'Lnm', 'min':1e4, 'max':npy.inf},
                 sort_arg = {'min':'mass'},
                 path='Export_Solution.txt', nb_sol=[20, 10, 10]):
        
        nb_linkage = len(length)
        self.path = path
        self.list_pos_unknown = list_pos_unknown
        self.list_load = list_load
        self.list_torque = list_torque
        self.list_speed = list_speed
        self.list_time = list_time
        self.d_shaft_min = d_shaft_min
        self.axial_pos = axial_pos
        self.d_ext = d_ext
        self.length = length
        self.typ_linkage = typ_linkage
        if typ_mounting == None:
            typ_mounting = []
            for typ_iter in product(['pn', 0, 'n', 'p'], repeat = nb_linkage):
                if typ_iter is not (0,)*nb_linkage:
                    if 'pn' in typ_iter:
                        if ('p' not in typ_iter) and ('n' not in typ_iter):
                            typ_mounting.append(typ_iter)
                    else:
                        typ_mounting.append(typ_iter)
        self.typ_mounting = typ_mounting
        for mount in self.typ_mounting:
            if set(mount) in [set(('pn', 'p')), set(('pn', 'n')), set(('pn', 'pn')), set((0, 0))]:
                raise KeyError('Modify the typ_linkage data')
                
        self.sort = sort
        list_linkage = self.typ_linkage 
        for num_link,list_link in enumerate(self.typ_linkage):
            if 'all' in list_link:
                list_linkage[num_link] = ['ball_joint', 'cylindric_joint']

        solutions = []
        for li_mounting in self.typ_mounting:
            sol_BA = {}
            for num_linkage, behavior_link in enumerate(li_mounting):
                sol_BA[num_linkage] = []
                for linkage, nb_rlts in product(list_linkage[num_linkage], 
                                                number_bearing[num_linkage]):
                    BA = BearingAssemblyOptimizer(linkage, behavior_link, nb_rlts, [self.d_shaft_min[num_linkage], 
                                                self.d_ext[num_linkage]], [self.d_shaft_min[num_linkage], 
                                                self.d_ext[num_linkage]], [0, self.length[num_linkage]],
                                                nb_sol = nb_sol[1::])
                    if BA.check is True:
                        sol_BA[num_linkage].extend(BA.solutions)
            for num_linkage, behavior_link in enumerate(li_mounting):
                analyze_solutions = True
                if sol_BA[num_linkage] == []:
                    analyze_solutions = False
            if analyze_solutions:
                for composite_solution in product(*sol_BA.values()):
                    solutions.append(CompositiveBearingAssembly(composite_solution))
        
        solutions = self.SortSolutions(solutions, sort_arg)
        self.solutions = solutions[0: min(len(solutions), nb_sol[0])]
        
    def SortSolutions(self, solutions, sort_arg):
        list_sort = []
        for composite_bg in solutions:
            val = getattr(composite_bg,sort_arg['min'])
            list_sort.append(val)
        list_sort_arg=list(npy.argsort(list_sort))
        return npy.array(solutions)[list_sort_arg]
        
    def Optimize(self, nb_sol=5, Lnm_min=1e4, index_sol=None, verbose=False):
        if index_sol is not None:
            list_sol = []
            for num_sol, sol in enumerate(self.solutions):
                if num_sol in index_sol:
                    list_sol.append(sol)
        else:
            list_sol = self.solutions
        self.solutions = []
        for composite_bg in list_sol:
            l1 = composite_bg.list_bearing_assembly[0].B
            l2 = composite_bg.list_bearing_assembly[1].B
            pos1_min = self.axial_pos[0] + l1/2
            pos1_max = self.axial_pos[0] + self.length[0] - l1/2
            pos2_min = self.axial_pos[1] + l2/2
            pos2_max = self.axial_pos[1] + self.length[1] - l2/2
            def fun(x):
                fa1, fr1, fa2, fr2 = composite_bg.ShaftLoad(x[0], x[1], self.list_pos_unknown, 
                                                    self.list_load, self.list_torque)
#                Lnm1 = 0
#                for rlt1 in list_bearing1:
#                    Lnm1 += rlt1.AdjustedLifeTime(Fr = [fr1], Fa = [fa1], N = [300], t = [1e10], T = [60])
#                Lnm2 = 0
#                for rlt2 in list_bearing2:
#                    Lnm2 += rlt2.AdjustedLifeTime(Fr = [fr2], Fa = [fa2], N = [300], t = [1e10], T = [60])
                obj = 1/(fr1)**2
                return obj
            def fineq(x):
                ineq = [0]
#                fa1, fr1, fa2, fr2 = composite_bg.ShaftLoad(x[0], x[1], self.list_pos_unknown, 
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
            
            if len(self.solutions) <= nb_sol:       
                if status >= 0:
                    composite_bg.Update(pos_x = sol_x)
                    self.solutions.append(composite_bg)
            else:
                break
    
        
class ShaftOptimizer:
    def __init__(self, list_pos_unknown, list_load, list_torque, list_speed, list_time,
                 d_shaft_min=0.02, axial_pos=[0, 0.1], d_ext=[0.05, 0.05], length=[0.04, 0.04],
                 typ_linkage=[['all'], ['all']],
                 typ_mounting=None, sort={'typ':'Lnm', 'min':1e4, 'max':npy.inf},
                 path='Export_Solution.txt'):
        file = open(path, 'w')
        file.write('''$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  \n''')
        file.write('''$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  \n''')
        file.write(''' _______                                 ______   ______        \n''')
        file.write('''|       \                               |      \ /      \       \n''')
        file.write('''| $$$$$$$\  ______    _______   _______  \$$$$$$|  $$$$$$\      \n''')
        file.write('''| $$  | $$ /      \  /       \ /       \  | $$  | $$__| $$      \n''')
        file.write('''| $$  | $$|  $$$$$$\|  $$$$$$$|  $$$$$$$  | $$  | $$    $$      \n''')
        file.write('''| $$  | $$| $$    $$ \$$    \  \$$    \   | $$  | $$$$$$$$      \n''')
        file.write('''| $$__/ $$| $$$$$$$$ _\$$$$$$\ _\$$$$$$\ _| $$_ | $$  | $$      \n''')
        file.write('''| $$    $$ \$$     \|       $$|       $$|   $$ \| $$  | $$      \n''')
        file.write(''' \$$$$$$$   \$$$$$$$ \$$$$$$$  \$$$$$$$  \$$$$$$ \$$   \$$      \n''')
        file.write(''' \n''')
        file.write('''$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  \n''')
        file.write('''$$$$$$$$$$$$$$$$$$     Bearings Results   $$$$$$$$$$$$$$$$$$$$  \n''')
        file.write('''$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  \n''')
        file.write(''' \n''')

        file.close()
        self.path = path
        
        self.list_pos_unknown = list_pos_unknown
        self.list_load = list_load
        self.list_torque = list_torque
        self.list_speed = list_speed
        self.list_time = list_time
        self.d_shaft_min = d_shaft_min
        self.axial_pos = axial_pos
        self.d_ext = d_ext
        self.length = length
        self.typ_linkage = typ_linkage
        if typ_mounting == None:
            typ_mounting = ['X', 'O', 'fix_sliding', 'sliding_fix']
        self.typ_mounting = typ_mounting
        self.sort = sort
        
        list_linkage = self.typ_linkage 
        for i,list_link in enumerate(self.typ_linkage):
            if 'all' in list_link:
                list_linkage[i] = ['linear_annular', 'ball_joint', 'pivot', 'cylindric_joint']
        
        list_mounting = []
        for typ_rlts1, typ_rlts2 in product(list_linkage[0], list_linkage[1]):
            if [typ_rlts1, typ_rlts2] == ['ball_joint', 'ball_joint']:
                if 'X' in self.typ_mounting:
                    list_mounting.append(['ball_joint', 'ball_joint', 'X'])
                if 'O' in self.typ_mounting:
                    list_mounting.append(['ball_joint', 'ball_joint', 'O'])
            if [typ_rlts1, typ_rlts2] == ['pivot', 'pivot']:
                if 'X' in self.typ_mounting:
                    list_mounting.append(['pivot', 'pivot', 'X'])
                if 'O' in self.typ_mounting:
                    list_mounting.append(['pivot', 'pivot', 'O'])
            if (typ_rlts1 in ['ball_joint', 'pivot']) and (typ_rlts2 in ['linear_annular', 'cylindric_joint']):
                if 'fix_sliding' in self.typ_mounting:
                    list_mounting.append([typ_rlts1, typ_rlts2, 'fix_sliding'])
            if (typ_rlts2 in ['ball_joint', 'pivot']) and (typ_rlts1 in ['linear_annular', 'cylindric_joint']):
                if 'sliding_fix' in self.typ_mounting:
                    list_mounting.append([typ_rlts1, typ_rlts2, 'sliding_fix'])

        axis_load = ['n', 'p'] # 'p' when the axial load applied on shaft is directed from left to right
        typ_rlts = {'rol_NU':0, 'rol_N':0, 'rol_NJ_p':'p', 'rol_NJ_n':'n',
                    'rol_NF_p':'p', 'rol_NF_n':'n', 'rol_NUP':'pn', 'ball':'pn', 
                    'ang_p':'p', 'ang_n':'n', 'tap_p':'p', 'tap_n':'n'}
        ring_mounting = ['n', 'p']
        solution_mounting = []
        compt = 0
        for linkage1, linkage2, typ_link in list_mounting:
            for nb_rlts1 in [1, 2]:
                if (linkage1 in ['ball_joint', 'linear_annular']) and (nb_rlts1 > 1):
                    continue
                for nb_rlts2 in [1, 2]:
                    if (linkage2 in ['ball_joint', 'linear_annular']) and (nb_rlts2 > 1):
                        continue
                    for li_rlts1 in product(typ_rlts.keys(), repeat = nb_rlts1):
                        if (linkage1 in ['ball_joint', 'linear_annular']) and (li_rlts1[0] not in ['ball', 'ang_p', 'ang_n']):
                            continue
                        li1 = []
                        for r1 in li_rlts1:
                            li1.append(typ_rlts[r1])
                        if (('n' in li1) and ('p' in li1)) or ('pn' in li1):
                            cond1 = 'pn'
                        elif 'p' in li1:
                            cond1 = 'p'
                        elif 'n' in li1:
                            cond1 = 'n'
                        else:
                            cond1 = 0
                            
                        if typ_link == 'fix_sliding':
                            if cond1 != 'pn':
                                continue
                            else:
                                mount = {1:{'be':['p','n'], 'bi':['p','n']}, 2:{'be':[], 'bi':[]}}
                        if typ_link == 'sliding_fix':
                            if (cond1 == 'n') or (cond1 == 'p'):
                                continue
                            else:
                                if cond1 == 0:
                                    mount = {1:{'be':['p','n'], 'bi':['p','n']}, 2:{'be':[], 'bi':[]}}
                                else:
                                    if (li1[0] == 0) or (li1[-1] == 0):
                                        continue
                                    else:
                                        mount = {1:{'be':[], 'bi':['p','n']}, 2:{'be':[], 'bi':[]}}
                        for li_rlts2 in product(typ_rlts.keys(), repeat = nb_rlts2):
                            if (linkage2 in ['ball_joint', 'linear_annular']) and (li_rlts2[0] not in ['ball', 'ang_p', 'ang_n']):
                                continue
                            # selection combinatory
                            all_rlts = li_rlts1 + li_rlts2
                            if ('ang_p' in all_rlts) or ('ang_n' in all_rlts):
                                if ('tap_p' in all_rlts) or ('tap_n' in all_rlts):
                                    continue
                            compt_rol = 0
                            for rl in li_rlts1:
                                if 'rol' in rl:
                                    compt_rol += 1
                            if compt_rol > 1:
                                continue
                            compt_rol = 0
                            for rl in li_rlts2:
                                if 'rol' in rl:
                                    compt_rol += 1
                            if compt_rol > 1:
                                continue
                                    
                            li2 = []
                            for r2 in li_rlts2:
                                li2.append(typ_rlts[r2])
                            if (('n' in li2) and ('p' in li2)) or ('pn' in li2):
                                cond2 = 'pn'
                            elif 'p' in li2:
                                cond2 = 'p'
                            elif 'n' in li2:
                                cond2 = 'n'
                            else:
                                cond2 = 0
                                
                            if typ_link == 'sliding_fix':
                                if cond2 != 'pn':
                                    continue
                                else:
                                    mount[2]['be']=['p','n']
                                    mount[2]['bi']=['p','n']
                            if typ_link == 'fix_sliding':
                                if (cond2 == 'n') or (cond2 == 'p'):
                                    continue
                                else:
                                    if cond2 == 0:
                                        mount[2]['be']=['p','n']
                                        mount[2]['bi']=['p','n']
                                    else:
                                        if (li2[0] == 0) or (li2[-1] == 0):
                                            continue
                                        else:
                                            mount[2]['be']=[]
                                            mount[2]['bi']=['p','n']
                            if typ_link == 'X':
                                if (cond1 != 'n') or (cond2 != 'p'):
                                    continue
                                else:
                                    mount = {1:{'be':['n'], 'bi':['p']}, 2:{'be':['p'], 'bi':['n']}}
                                    if li1[0] == 0:
                                        mount[1]['bi'].append('n')
                                    if li1[-1] == 0:
                                        mount[1]['be'].append('p')
                                    if li2[0] == 0:
                                        mount[2]['be'].append('n')
                                    if li2[-1] == 0:
                                        mount[2]['bi'].append('p')
                            if typ_link == 'O':
                                if (cond1 != 'p') or (cond2 != 'n'):
                                    continue
                                else:
                                    mount = {1:{'be':['p'], 'bi':['n']}, 2:{'be':['n'], 'bi':['p']}}
                                    if li1[0] == 0:
                                        mount[1]['be'].append('n')
                                    if li1[-1] == 0:
                                        mount[1]['bi'].append('p')
                                    if li2[0] == 0:
                                        mount[2]['bi'].append('n')
                                    if li2[-1] == 0:
                                        mount[2]['be'].append('p')
                            compt += 1
                            mount[1]['list_rlts'] = li_rlts1
                            mount[2]['list_rlts'] = li_rlts2
                            mount['typ'] = typ_link
                            solution_mounting.append(mount)

        self.solution_mounting = solution_mounting
                    
    
    
    def Plot(self, solution):
        
        length_bg1 = 0
        for bg in solution['bearing1']:
            length_bg1 += bg.B
        length_bg2 = 0
        for bg in solution['bearing2']:
            length_bg2 += bg.B
            
        def BearingFixation(be, bi, D, d, B, pos, sign):
            ep = min(0.1*D, 0.1*d)
            De = D + ep
            Dg, Dd = D, D
            if 'n' in be:
                Dg = Dg - ep
            if 'p' in be:
                Dd = Dd -ep
            di = d - ep
            dg, dd = d, d
            if 'n' in bi:
                dg = dg + ep
            if 'p' in bi:
                dd = dd + ep
            be = vm.Polygon2D([vm.Point2D((-B/2., sign*D/2.)), vm.Point2D((-B/2., sign*Dg/2.)),
                               vm.Point2D((-B/2. - ep, sign*Dg/2.)), vm.Point2D((-B/2. - ep, sign*De/2.)),
                               vm.Point2D((B/2. + ep, sign*De/2.)), vm.Point2D((B/2. + ep, sign*Dd/2.)),
                               vm.Point2D((B/2., sign*Dd/2.)), vm.Point2D((B/2., sign*D/2.))])
            bi = vm.Polygon2D([vm.Point2D((-B/2., sign*d/2.)), vm.Point2D((-B/2., sign*dg/2.)),
                               vm.Point2D((-B/2. - ep, sign*dg/2.)), vm.Point2D((-B/2. - ep, sign*di/2.)),
                               vm.Point2D((B/2. + ep, sign*di/2.)), vm.Point2D((B/2. + ep, sign*dd/2.)),
                               vm.Point2D((B/2., sign*dd/2.)), vm.Point2D((B/2., sign*d/2.))])
            be.Translation((pos, 0))
            bi.Translation((pos, 0))
            return vm.Contour2D([be, bi])
        
        bearing1_1_fixation = BearingFixation(solution['mounting'][1]['be'], solution['mounting'][1]['bi'], solution['bearing1'][0].D,
                                            solution['bearing1'][0].d, length_bg1, solution['pos_bearing1'], sign = 1)
        bearing1_2_fixation = BearingFixation(solution['mounting'][1]['be'], solution['mounting'][1]['bi'], solution['bearing1'][0].D,
                                            solution['bearing1'][0].d, length_bg1, solution['pos_bearing1'], sign = -1)
        bearing2_1_fixation = BearingFixation(solution['mounting'][2]['be'], solution['mounting'][2]['bi'], solution['bearing2'][0].D,
                                            solution['bearing2'][0].d, length_bg2, solution['pos_bearing2'], sign = 1)
        bearing2_2_fixation = BearingFixation(solution['mounting'][2]['be'], solution['mounting'][2]['bi'], solution['bearing2'][0].D,
                                            solution['bearing2'][0].d, length_bg2, solution['pos_bearing2'], sign = -1)
            
        pos1 = solution['pos_bearing1'] - length_bg1/2.
        pos2 = solution['pos_bearing2'] - length_bg2/2.
        list_bg = []
        for bg in solution['bearing1']:
            bg1 = bg.Plot(pos = pos1 + bg.B/2.)
            pos1 += bg.B
            list_bg.append(bg1)
        for bg in solution['bearing2']:
            bg2 = bg.Plot(pos = pos2 + bg.B/2.)
            pos2 += bg.B
            list_bg.append(bg2)
            
        poly_shaft = vm.Polygon2D([vm.Point2D((self.axial_pos[0] - 5e-3, -self.d_shaft_min[0]/2.)),
                      vm.Point2D((self.axial_pos[0] - 5e-3, self.d_shaft_min[0]/2.)),
                      vm.Point2D(((self.axial_pos[0] + self.length[0] + self.axial_pos[1])/2., self.d_shaft_min[0]/2.)),
                      vm.Point2D(((self.axial_pos[0] + self.length[0] + self.axial_pos[1])/2., self.d_shaft_min[1]/2.)),
                      vm.Point2D((self.axial_pos[1] + self.length[1] + 5e-3, self.d_shaft_min[1]/2.)),
                      vm.Point2D((self.axial_pos[1] + self.length[1] + 5e-3, -self.d_shaft_min[1]/2.)),
                      vm.Point2D(((self.axial_pos[0] + self.length[0] + self.axial_pos[1])/2., -self.d_shaft_min[1]/2.)),
                      vm.Point2D(((self.axial_pos[0] + self.length[0] + self.axial_pos[1])/2., -self.d_shaft_min[0]/2.)),])
    
        poly_rlt1 = vm.Polygon2D([vm.Point2D((self.axial_pos[0], self.d_shaft_min[0]/2.)),
                      vm.Point2D((self.axial_pos[0], self.d_ext[0]/2.)),
                      vm.Point2D((self.axial_pos[0] + self.length[0], self.d_ext[0]/2.)),
                      vm.Point2D((self.axial_pos[0] + self.length[0], self.d_shaft_min[0]/2.))])
        poly_rlt1b = poly_rlt1.Translation((0, -self.d_ext[0]/2. - self.d_shaft_min[0]/2.), True)
    
        poly_rlt2 = vm.Polygon2D([vm.Point2D((self.axial_pos[1], self.d_shaft_min[1]/2.)),
                      vm.Point2D((self.axial_pos[1], self.d_ext[1]/2.)),
                      vm.Point2D((self.axial_pos[1] + self.length[1], self.d_ext[1]/2.)),
                      vm.Point2D((self.axial_pos[1] + self.length[1], self.d_shaft_min[1]/2.))])
        poly_rlt2b = poly_rlt2.Translation((0, -self.d_ext[1]/2. - self.d_shaft_min[1]/2.), True)

        list_bg.append(poly_shaft)
        contour_bg = vm.Contour2D(list_bg)
        contour_area = vm.Contour2D([poly_rlt1, poly_rlt2, poly_rlt1b, poly_rlt2b])
        f, a = contour_area.MPLPlot(style = '-r')
        
        bearing1_1_fixation.MPLPlot(a,'-g')
        bearing1_2_fixation.MPLPlot(a,'-g')
        bearing2_1_fixation.MPLPlot(a,'-g')
        bearing2_2_fixation.MPLPlot(a,'-g')
        
        total_length = (solution['pos_bearing1'] - length_bg1/2. + solution['pos_bearing2'] + length_bg2/2.)/2.
        fa1, fr1, fa2, fr2 = solution['load']
        ratio_long = 0.2*total_length/max([fa1, fr1, fa2, fr2])
        
        if fa1 !=0:
            a.arrow(solution['pos_bearing1'],0,ratio_long*fa1,0)
        if fa2 != 0:
            a.arrow(solution['pos_bearing2'],0,ratio_long*fa2,0)
        a.arrow(solution['pos_bearing1'],0,0,ratio_long*fr1)
        a.arrow(solution['pos_bearing2'],0,0,ratio_long*fr2)
        
        a.arrow(self.list_pos_unknown[0][0],(self.list_pos_unknown[0][1]**2 + self.list_pos_unknown[0][2]**2)**0.5,ratio_long*self.list_load[0][0],ratio_long*(self.list_load[0][1]**2 + self.list_load[0][2]**2)**0.5)
    
        contour_bg.MPLPlot(a)
        
    def Export(self, solution, num_sol):
        R1 = BearingCatalogOptimizer()
        file = open(self.path, 'a')
        file.write('''$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  \n''')
        file.write('''$$$$$$$$$$$$$$$$$     Optimum number {}   $$$$$$$$$$$$$$$$$$$$  \n'''.format(num_sol))
        file.write('''$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  \n''')
        file.write(''' \n''')
        file.write('''position; name bearing; type bearing, d, D, B, i, Z, Dw, alpha  \n''')
        for num_ligne in solution['bearing1_catalog']:
            name_bearing = R1.pd_FNR.loc[num_ligne,'name_bearing']
            typ_rlt2 = R1.pd_FNR.loc[num_ligne,'typ_bearing']
            d = R1.pd_FNR.loc[num_ligne,'d']
            D = R1.pd_FNR.loc[num_ligne,'D']
            B = R1.pd_FNR.loc[num_ligne,'B']
            i = R1.pd_FNR.loc[num_ligne,'i']
            Z = R1.pd_FNR.loc[num_ligne,'Z']
            Dw = R1.pd_FNR.loc[num_ligne,'Dw']
            alpha = R1.pd_FNR.loc[num_ligne,'alpha']
            file.write('''{}; {}; {}; {}; {}; {}; {}; {}; {}; {} \n'''.format('left', name_bearing, typ_rlt2, d, D, B, i, Z, Dw, alpha))
        for num_ligne in solution['bearing2_catalog']:
            name_bearing = R1.pd_FNR.loc[num_ligne,'name_bearing']
            typ_rlt2 = R1.pd_FNR.loc[num_ligne,'typ_bearing']
            d = R1.pd_FNR.loc[num_ligne,'d']
            D = R1.pd_FNR.loc[num_ligne,'D']
            B = R1.pd_FNR.loc[num_ligne,'B']
            i = R1.pd_FNR.loc[num_ligne,'i']
            Z = R1.pd_FNR.loc[num_ligne,'Z']
            Dw = R1.pd_FNR.loc[num_ligne,'Dw']
            alpha = R1.pd_FNR.loc[num_ligne,'alpha']
            file.write('''{}; {}; {}; {}; {}; {}; {}; {}; {}; {} \n'''.format('right', name_bearing, typ_rlt2, d, D, B, i, Z, Dw, alpha))
        file.close()
        

#print(S1.list_solutions)
    
#B1 = BearingCatalogOptimizer()
#df=B1.Optimization([0.02,0.03],[0.035,0.045],[0.01,0.02],100,100,1,1,1)
#print(df.loc[:,['d','D','B']])