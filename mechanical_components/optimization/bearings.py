#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 02:14:21 2018

@author: steven
"""

from mechanical_components.bearings import oil_iso_vg_1500, material_iso, iso_bearings, bearing_rules, iso_rollers

from mechanical_components.bearings import RadialBallBearing, AngularBallBearing, \
        SphericalBallBearing, RadialRollerBearing, TaperedRollerBearing, \
        NUPRadialRollerBearing, NRadialRollerBearing, \
        NFRadialRollerBearing, NURadialRollerBearing, \
        BearingAssembly, BearingCombination, DetailedRadialRollerBearing, \
        BearingAssemblySimulationResult, \
        BearingCombinationSimulationResult, BearingSimulationResult,\
        BearingAssemblySimulation, bearing_classes, dict_bearing_classes, \
        ConceptualBearingCombination, BearingCatalogue
        
import numpy as npy

npy.seterr(divide='raise', over='ignore', under='ignore', invalid='ignore')

from scipy.optimize import minimize, fsolve
import pandas
from copy import deepcopy
#from pandas.plotting import scatter_matrix
import pkg_resources
from itertools import product

from dectree import DecisionTree

#import genmechanics
#import genmechanics.linkages as linkages
#import genmechanics.loads as loads
import math
#import volmdlr as vm
#import volmdlr.primitives3D as primitives3D
#import volmdlr.primitives2D as primitives2D

from mechanical_components.tools import StringifyDictKeys

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
                sol=[-math.inf,math.inf]
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
        for item_d in iso_bearings.keys():
            if (item_d<=interval_d[1]) and (item_d>=interval_d[0]):
                lk_d.append([item_d])
        liste_keys=lk_d
        lk_D=[]
        for [item_d] in liste_keys:
            for item_D in iso_bearings[item_d].keys():
                if (item_D<=interval_D[1]) and (item_D>=interval_D[0]):
                    lk_D.append([item_d,item_D])
        liste_keys=lk_D
        lk_B=[]
        for item_d,item_D in liste_keys:
            for item_B in iso_bearings[item_d][item_D].keys():
                if (item_B<=interval_B[1]) and (item_B>=interval_B[0]):
                    lk_B.append([item_d,item_D,item_B])
        liste_keys=lk_B
        liste_sol_iso=[]
        for item_d,item_D,item_B in liste_keys:
            dk=[item_d,item_D,item_B]
            dk.append(list(iso_bearings[item_d][item_D][item_B].keys())[0])
            dk.append(list(iso_bearings[item_d][item_D][item_B].values())[0])
            liste_sol_iso.append(dk) #Liste dk contient [d,D,B,rsmin,serial]
        
        #Approche combinatoire sur les rouleaux
        def analyse_rule(Y,inter_min,inter_max,input_var,data_var):
            for (var_x,var_y,typ),[a,b] in bearing_rules.items():
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
            E_inf,E_sup=-math.inf,math.inf
            F_inf,F_sup=-math.inf,math.inf
            E_inf,E_sup=analyse_rule('E',E_inf,E_sup,['d','D','B'],[d,D,B])
            F_inf,F_sup=analyse_rule('F',F_inf,F_sup,['d','D','B'],[d,D,B])
                            
            #Analyse borne sup/inf de Lw et Dw
            Lw_inf,Lw_sup=-math.inf,math.inf
            Dw_inf,Dw_sup=E_inf-F_sup,E_sup-F_inf
            Lw_inf,Lw_sup=analyse_rule('Lw',Lw_inf,Lw_sup,['d','D','B'],[d,D,B])
            Dw_inf,Dw_sup=analyse_rule('Dw',Dw_inf,Dw_sup,['d','D','B'],[d,D,B])
            
            #Choix des dimensions rouleaux compatible
            Dw_sup=min(Dw_sup,(D-d)/2.)
            for Dw in iso_rollers.keys():
                if (Dw>=Dw_inf) and (Dw<=Dw_sup):
                    for Lw in iso_rollers[Dw].keys():
                        if (Lw>=Lw_inf) and (Lw<=Lw_sup):
                            E_inf,E_sup=analyse_rule('E',-math.inf,math.inf,['d','D','B','Dw','Lw'],[d,D,B,Dw,Lw])
                            F_inf,F_sup=analyse_rule('F',-math.inf,math.inf,['d','D','B','Dw','Lw'],[d,D,B,Dw,Lw])
                            if ((E_inf-F_sup)/2.<=Dw) and ((E_sup-F_inf)/2.>=Dw):
                                liste_sol_roller_iso.append([d,D,B,rsmin,serial,Dw,Lw,E_inf,E_sup,F_inf,F_sup])


        #Analyse de détail des roulements
        a_ED1_inf,b_ED1_inf=bearing_rules[('E','D1','inf')]
        a_ED1_sup,b_ED1_sup=bearing_rules[('E','D1','sup')]
        a_Fd1_inf,b_Fd1_inf=bearing_rules[('F','d1','inf')]
        a_Fd1_sup,b_Fd1_sup=bearing_rules[('F','d1','sup')]
        estim_masse=[]
        #Construction d'une fonctionnelle pour le tri afin d'engager une optimisation continue
        for [d,D,B,rsmin,serial,Dw,Lw,E_inf,E_sup,F_inf,F_sup] in liste_sol_roller_iso:
            E=(E_inf+E_sup)/2.
            F=E-2*Dw
            Zmax=int(2*math.pi/(2.*math.arcsin((Dw/2)/(F/2.+Dw/2.))))
            masse_elem=(math.pi*D**2/4.-math.pi*E**2/4.)*B
            masse_elem+=(math.pi*F**2/4.-math.pi*d**2/4.)*B
            masse_elem+=math.pi*Dw**2/4.*Lw*Zmax
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
                            
with pkg_resources.resource_stream(pkg_resources.Requirement('mechanical_components'),
                                   'mechanical_components/catalogs/SNR.csv') as rlts_FNR:
    pandas_rlts_FNR = pandas.read_csv(rlts_FNR) 
            
pandas_sort = pandas_rlts_FNR[pandas_rlts_FNR['i'].notnull()]
pandas_sort = pandas_sort[pandas_sort['Z'].notnull()]
pandas_sort = pandas_sort[pandas_sort['Dw'].notnull()]
pandas_sort = pandas_sort[pandas_sort['mass'].notnull()]
#        pandas_sort = pandas_sort[pandas_sort['alpha'].notnull()]
pandas_sort = pandas_sort[pandas_sort['Cr'].notnull()]
pandas_sort = pandas_sort[pandas_sort['C0r'].notnull()]
base_bearing = pandas_sort


class BearingCombinationOptimizer:
    
    def __init__(self, linkage, mounting, d, D, length,
                 bearing_classes):
        
        self.bearing_classes = bearing_classes
        self.linkage = linkage
        self.mounting = mounting
        self.d = d
        self.D = D
        self.length = length
        
        
    def ConceptualBearingCombinations(self, max_bearings=3):
        configurations = []
        
        dt = DecisionTree()
        nclasses = len(self.bearing_classes)
        while not dt.finished:
            # Constructing BearingCombination
            valid = True
            bearings = []
            directions = []
            for depth, node in enumerate(dt.current_node):
                if depth % 2 == 0:
                    bearing = self.bearing_classes[node]
                    bearings.append(bearing)
                else:
                    directions.append([1, -1][node])

            if ((dt.current_depth + 1) % 2 != 0) and dt.current_depth > 1:
                # Instanciating 
                cbc = ConceptualBearingCombination(bearings, directions, self.mounting)
                valid = cbc.CheckKinematic()
            
            # Testing
            if valid:
            
                # Counting possibilities
                if dt.current_depth % 2 != 0:
                    # counting sides
                    if bearing.symmetric:
                        dt.SetCurrentNodeNumberPossibilities(1)
                    else:
                        dt.SetCurrentNodeNumberPossibilities(2)
                else:
                    if dt.current_depth // 2 == max_bearings:
                        dt.SetCurrentNodeNumberPossibilities(0)    
                        configurations.append(cbc)
                    else:
                        dt.SetCurrentNodeNumberPossibilities(nclasses)
            else:
                dt.SetCurrentNodeNumberPossibilities(0) 
            
            dt.NextNode(valid)
                            
        return configurations
    
class BearingAssemblyOptimizer:
    dessia_db_attributes = [{'name':'bearing_assembly_simulations',
                             'class':'mechanical_components.bearings.BearingAssemblySimulation',
                             'type':'list'}]
    
    def __init__(self, loads, speeds, operating_times,
                 inner_diameters,
                 outer_diameters,
                 axial_positions,
                 lengths,
                 linkage_types=[['ball_joint', 'cylindric_joint'],
                                ['ball_joint', 'cylindric_joint']],
                 mounting_types=list(product(['both', 'right', 'left', 'free'],
                                             ['both', 'right', 'left', 'free'])),
                 number_bearings=[[1, 2], [1, 2]],
                 bearing_classes=bearing_classes,
                 bearing_assembly_simulations=None):
        self.loads = loads
        self.speeds = speeds
        self.operating_times = operating_times
        self.inner_diameters = inner_diameters
        self.axial_positions = axial_positions
        self.outer_diameters = outer_diameters
        self.lengths = lengths
        self.linkage_types = linkage_types
        self.mounting_types = mounting_types
        self.number_bearings = number_bearings
        self.bearing_classes = bearing_classes
        self.bearing_assembly_simulations = bearing_assembly_simulations
        self.bearing_catalogue = BearingCatalogue()
        
    def Configurations(self):
        configurations = []

        # test axial load sign
        list_load_cases = []
        for load_cases in self.loads:
            axial_load = 0
            for pos, ld, tq in load_cases:    
                axial_load += ld[0]
            if axial_load > 0:
                list_load_cases.append('right')
            elif axial_load < 0:
                list_load_cases.append('left')
        
        for mounting_type_left, mounting_type_right in self.mounting_types:
            valid_load_case = True
            if 'both' not in set((mounting_type_left, mounting_type_right)):
                for load_case in list_load_cases:
                    if load_case not in (mounting_type_left, mounting_type_right):
                        valid_load_case = False
            if valid_load_case:
                for linkage_left, linkage_right in product(*self.linkage_types):
                    configurations.append(((mounting_type_left, linkage_left),
                                           (mounting_type_right, linkage_right)))
        return configurations


    def ConceptualBearingCombinations(self, max_bearings=[2, 2]):
        bearing_combinations_possibilities = {}
        for left, right in self.Configurations():
            if not (left, right) in bearing_combinations_possibilities:
                DBC_l = BearingCombinationOptimizer(left[1], 
                                                    left[0], 
                                                    self.inner_diameters[0],
                                                    self.outer_diameters[0], 
                                                    self.lengths[0], 
                                                    bearing_classes = self.bearing_classes
                                                    )
                DBC_r = BearingCombinationOptimizer(right[1], 
                                                      right[0], 
                                                      self.inner_diameters[1],
                                                      self.outer_diameters[1], 
                                                      self.lengths[1], 
                                                      bearing_classes = self.bearing_classes
                                                     )
                bearing_combinations_possibilities[(left, right)] = \
                    (DBC_l.ConceptualBearingCombinations(max_bearings[0]), DBC_r.ConceptualBearingCombinations(max_bearings[1]))
        self.bearing_combinations_possibilities = bearing_combinations_possibilities
        
    def SelectBestBearingCombinations(self, first_bearing_possibilies, conceptual_bearing_combination, 
                                      radial_elementary_load, L10_objective):
        bearing_possibilies = []
        for i, bearing in enumerate(first_bearing_possibilies):
            list_L10_max = []
            bearings = []
            
            d = bearing.d
            D = bearing.D
            L10 = bearing.EstimateBaseLifeTime(Fr = radial_elementary_load,
                                                            N = self.speeds, 
                                                            t = self.operating_times, Cr = bearing.Cr)
            list_L10_max.append(L10)
            bearings.append(bearing)
            if len(conceptual_bearing_combination.bearing_classes) > 1:
                for conceptual_bearing in conceptual_bearing_combination.bearing_classes[1:]:
                    try:
                        bearing_possibilities = self.bearing_catalogue.SearchDictBearingCatalogue(conceptual_bearing, d , D)
                        list_L10_max.append(bearing_possibilities[-1].EstimateBaseLifeTime(Fr = radial_elementary_load,
                                                                N = self.speeds, 
                                                                t = self.operating_times, Cr = bearing_possibilities[-1].Cr))
                        bearings.append(bearing_possibilities[-1])
                    except:
                        break
            
            if len(list_L10_max) == len(conceptual_bearing_combination.bearing_classes):
                L10_max = BearingAssembly.EstimateBaseLifeTime(list_L10_max)
                if L10_max > L10_objective:
                    bearing_possibilies.append(bearing)
        return bearing_possibilies, bearings
        
    def SelectBearingCombinations(self, bearing_combinations_possibility, radial_load, 
                                  mounting, L10_objective,
                                  max_bearing_assemblies=10):
        
        radial_load_left, radial_load_right = radial_load
        left, right = mounting
        
        compt_continue = 0
        for conceptual_bearing_combination_left, conceptual_bearing_combination_right in \
                product(*bearing_combinations_possibility):
                    
            dt = DecisionTree()
            
            first_bearing_left_possibilies = self.bearing_catalogue\
                .SearchBearingCatalogue(conceptual_bearing_combination_left.bearing_classes[0], self.inner_diameters[0], self.outer_diameters[0])
            first_bearing_right_possibilies = self.bearing_catalogue\
                .SearchBearingCatalogue(conceptual_bearing_combination_right.bearing_classes[0], self.inner_diameters[1], self.outer_diameters[1])
            nb_bearings_left = len(conceptual_bearing_combination_left.bearing_classes)
            nb_bearings_right = len(conceptual_bearing_combination_right.bearing_classes)
            nb_bearings = nb_bearings_left + nb_bearings_right
            
            # selection of the best first_bearing_left_possibilies
            bearing_left_possibilies, bearings_left_max = self.SelectBestBearingCombinations(first_bearing_left_possibilies,
                                               conceptual_bearing_combination_left,
                                               [r/(1.*nb_bearings_left) for r in radial_load_left],
                                               L10_objective/2.)
                    
            # selection of the best first_bearing_right_possibilies
            bearing_right_possibilies, bearings_right_max = self.SelectBestBearingCombinations(first_bearing_right_possibilies,
                                               conceptual_bearing_combination_right,
                                               [r/(1.*nb_bearings_right) for r in radial_load_right],
                                               L10_objective/2.)
            
            bc_left = conceptual_bearing_combination_left.BearingCombination(bearings_left_max)
            bc_right = conceptual_bearing_combination_right.BearingCombination(bearings_right_max)
            bearing_assembly = BearingAssembly([bc_left, bc_right])
            bc_results = []
            for bearing_combination in bearing_assembly.bearing_combinations:
                li_bg_results = []
                for bearing in bearing_combination.bearings:
                    li_bg_results.append(BearingSimulationResult())
                bc_results.append(BearingCombinationSimulationResult(li_bg_results))
            bearing_assembly_simulation_result = BearingAssemblySimulationResult(bc_results, 
                                                            self.loads, self.speeds, self.operating_times)
            pos1_min = self.axial_positions[0]
            pos1_max = self.axial_positions[0] + self.lengths[0]
            pos2_min = self.axial_positions[1]
            pos2_max = self.axial_positions[1] + self.lengths[1]
            bearing_assembly.ShaftLoad([(pos1_min + pos1_max)/2., (pos2_min + pos2_max)/2.], 
                                        bearing_assembly_simulation_result)
            L10 = bearing_assembly_simulation_result.L10
            if L10 < L10_objective:
                compt_continue += 1
                if compt_continue > 10:
                    break
                continue
                        
            list_bearing_possibilities = [0]*nb_bearings
            list_Cr = []

            compt_same_configuration = 0
            while not dt.finished:
                # Constructing BearingCombination
                valid = True
                bearings = []
                B_left = 0
                B_right = 0
                Cr_current_node = []
                
                for depth, node in enumerate(dt.current_node):
                    if depth == 0:
                        bearing = bearing_left_possibilies[node]
                        bearings.append(bearing)
                        d = bearing.d
                        d_left = d
                        D = bearing.D
                        B_left += bearing.B
                        
                    elif depth < nb_bearings_left:
                        bearing_classe = conceptual_bearing_combination_left.bearing_classes[depth]
                        bearing_possibilities = list_bearing_possibilities[depth]
                        bearings.append(bearing_possibilities[node])
                        B_left += bearing.B
                        
                    elif depth == nb_bearings_left:
                        bearing = bearing_right_possibilies[node]
                        bearings.append(bearing)
                        d = bearing.d
                        d_right = d
                        D = bearing.D
                        B_right += bearing.B
                        if d_right < d_left:
                            valid = False
                            break
                        
                    elif depth < nb_bearings:
                        bearing_classe = conceptual_bearing_combination_right.bearing_classes[depth - nb_bearings_left]
                        bearing_possibilities = list_bearing_possibilities[depth]
                        bearings.append(bearing_possibilities[node])
                        B_right += bearing.B
                    
                    Cr_current_node.append(bearings[-1].Cr)
                    
                    if len(list_Cr) > 0:
                        if Cr_current_node[-1] < min([cr[depth] for cr in list_Cr]):
                            valid = False
                            break
                        
                if (dt.current_depth == nb_bearings) and valid:
                        
                    if (B_left > self.lengths[0]) or (B_right > self.lengths[1]):
                        valid = False
    
                if (dt.current_depth == nb_bearings) and valid:
                    
                    bc_left = conceptual_bearing_combination_left.BearingCombination(bearings[0: nb_bearings_left])
                    bc_right = conceptual_bearing_combination_right.BearingCombination(bearings[nb_bearings_left:])
                    bearing_assembly = BearingAssembly([bc_left, bc_right])
                    bc_results = []
                    for bearing_combination in bearing_assembly.bearing_combinations:
                        li_bg_results = []
                        for bearing in bearing_combination.bearings:
                            li_bg_results.append(BearingSimulationResult())
                        bc_results.append(BearingCombinationSimulationResult(li_bg_results))
                    bearing_assembly_simulation_result = BearingAssemblySimulationResult(bc_results, 
                                                                    self.loads, self.speeds, self.operating_times)
                    pos1_min = self.axial_positions[0]
                    pos1_max = self.axial_positions[0] + self.lengths[0]
                    pos2_min = self.axial_positions[1]
                    pos2_max = self.axial_positions[1] + self.lengths[1]
                    bearing_assembly.ShaftLoad([(pos1_min + pos1_max)/2., (pos2_min + pos2_max)/2.], 
                                                bearing_assembly_simulation_result)
                    L10 = bearing_assembly_simulation_result.L10
                    
                    if L10 < L10_objective:
#                        print(L10, L10_objective, Cr_current_node, dt.current_node)
                        valid = False
                        
                        def funct(alpha):
                            bgs = [deepcopy(bg) for bg in bearings]
                            for bg in bgs:
                                bg.Cr = bg.Cr*alpha[0]
                            bc_left = conceptual_bearing_combination_left.BearingCombination(bgs[0: nb_bearings_left])
                            bc_right = conceptual_bearing_combination_right.BearingCombination(bgs[nb_bearings_left:])
                            bearing_assembly = BearingAssembly([bc_left, bc_right])
                            bc_results = []
                            for bearing_combination in bearing_assembly.bearing_combinations:
                                li_bg_results = []
                                for bearing in bearing_combination.bearings:
                                    li_bg_results.append(BearingSimulationResult())
                                bc_results.append(BearingCombinationSimulationResult(li_bg_results))
                            bearing_assembly_simulation_result = BearingAssemblySimulationResult(bc_results, 
                                                                            self.loads, self.speeds, self.operating_times)
                            pos1_min = self.axial_positions[0]
                            pos1_max = self.axial_positions[0] + self.lengths[0]
                            pos2_min = self.axial_positions[1]
                            pos2_max = self.axial_positions[1] + self.lengths[1]
                            bearing_assembly.ShaftLoad([(pos1_min + pos1_max)/2., (pos2_min + pos2_max)/2.], 
                                                        bearing_assembly_simulation_result)
                            L10 = bearing_assembly_simulation_result.L10
                            
                            return (L10 - L10_objective)**2
                        
                        valid_fsolve = False
                        for i in range(5):
                            cond_init = npy.random.random(1)*1e4
                            sol = fsolve(funct, cond_init[0])[0]
                            if funct([sol]) < 1e-4:
                                coefficient_Cr = sol
                                valid_fsolve = True
#                                print('analyse coeff ', coefficient_Cr, funct([coefficient_Cr]))
                                break
                        if not valid_fsolve:
                            break
                        

                        
                        Cr_current_node_max = []
                        for bg in bearings:
                            Cr_current_node_max.append(bg.Cr*coefficient_Cr)
                        valid_Cr = True
                        for Crs in list_Cr:
                            if min(npy.array(Crs) - npy.array(Cr_current_node_max)) >= 0:
                                valid_Cr = False
                        if valid_Cr:
                            list_Cr.append(Cr_current_node_max)
                            
                # Testing
                if valid:
                    # Counting possibilities
                    if dt.current_depth == 0:
                        dt.SetCurrentNodeNumberPossibilities(len(bearing_left_possibilies))
                        list_bearing_possibilities[dt.current_depth] = bearing_left_possibilies
                    elif dt.current_depth < nb_bearings_left:
                        bearing_classe = conceptual_bearing_combination_left.bearing_classes[dt.current_depth]
                        bearing = bearing_left_possibilies[dt.current_node[0]]
                        d = bearing.d
                        D = bearing.D
                        try:
                            bearing_possibilities = self.bearing_catalogue.SearchDictBearingCatalogue(bearing_classe, d , D)
                            dt.SetCurrentNodeNumberPossibilities(len(bearing_possibilities))
                            list_bearing_possibilities[dt.current_depth] = bearing_possibilities
                        except:
                            valid = False
                            dt.SetCurrentNodeNumberPossibilities(0)
                    elif dt.current_depth == nb_bearings_left:
                        dt.SetCurrentNodeNumberPossibilities(len(bearing_right_possibilies))
                        list_bearing_possibilities[dt.current_depth] = bearing_right_possibilies
                    elif dt.current_depth < nb_bearings:
                        bearing_classe = conceptual_bearing_combination_right.bearing_classes[dt.current_depth - nb_bearings_left]
                        bearing = bearing_right_possibilies[dt.current_node[nb_bearings_left]]
                        d = bearing.d
                        D = bearing.D
                        try:
                            bearing_possibilities = self.bearing_catalogue.SearchDictBearingCatalogue(bearing_classe, d , D)
                            dt.SetCurrentNodeNumberPossibilities(len(bearing_possibilities))
                            list_bearing_possibilities[dt.current_depth] = bearing_possibilities
                        except:
                            valid = False
                            dt.SetCurrentNodeNumberPossibilities(0)
                    elif dt.current_depth == nb_bearings:
                        dt.SetCurrentNodeNumberPossibilities(0) 
                        bc_left = conceptual_bearing_combination_left.BearingCombination(bearings[0: nb_bearings_left])
                        bc_right = conceptual_bearing_combination_right.BearingCombination(bearings[nb_bearings_left:])
                        ba = BearingAssembly([bc_left, bc_right])
                        compt_same_configuration += 1
                        yield ba
                        if compt_same_configuration > 0:
                            break
                    else:
                        dt.SetCurrentNodeNumberPossibilities(0)      
                else:
                    dt.SetCurrentNodeNumberPossibilities(0)
                
                dt.NextNode(valid)
    
    def Optimize(self, max_solutions=10):
        
        L10_objective = 0
        for speed, time in zip(self.speeds, self.operating_times):
            L10_objective += speed/(2*math.pi)*time
        L10_objective = L10_objective/1e6
        print('the L10 objective is {}'.format(L10_objective))
        
        pos1_min = self.axial_positions[0]
        pos1_max = self.axial_positions[0] + self.lengths[0]
        pos2_min = self.axial_positions[1]
        pos2_max = self.axial_positions[1] + self.lengths[1]
        load = BearingAssembly.QuickShaftLoad(((pos1_min + pos1_max)/2., (pos2_min + pos2_max)/2.), self.loads)
        radial_load_left = []
        radial_load_right = []
        for rl in load:
            radial_load_left.append((rl[0]**2 + (0)**2)**0.5)
            radial_load_right.append((rl[1]**2 + (0)**2)**0.5)
            
            
        bearing_assembly_simulations = []
        sort_bearing_assembly_simulations = []
        combination_number_bearings = list(product(*self.number_bearings))
        combination_number_bearings = [combination_number_bearings[j] for j in npy.argsort([sum(i) for i in combination_number_bearings])]
        for max_bearings_left, max_bearings_right in combination_number_bearings:
            print('number of bearings analyzed: {} left and {} right'.format(max_bearings_left, max_bearings_right))
            self.ConceptualBearingCombinations(max_bearings = (max_bearings_left, max_bearings_right))
            for (left, right), bearing_combinations_possibility in self.bearing_combinations_possibilities.items():
                print((left, right))
                bearing_assemblies = self.SelectBearingCombinations(bearing_combinations_possibility, 
                                                                    radial_load = (radial_load_left, radial_load_right), 
                                                                    mounting = (left, right),
                                                                    L10_objective = L10_objective,
                                                                    max_bearing_assemblies = max_solutions)

                for i_bearing_assembly, bearing_assembly in enumerate(bearing_assemblies):
#                    print('try continuous optim ', i_bearing_assembly)
                    bearing_assembly_simulation = self.ContinuousOptimize(bearing_assembly)
                    if bearing_assembly_simulation != False:
                        L10 = bearing_assembly_simulation.bearing_assembly_simulation_result.L10
                        print(L10, L10_objective)
                        if L10 >= L10_objective:
                            bearing_assembly_simulations.append(bearing_assembly_simulation)
                            sort_bearing_assembly_simulations.append(L10)
                            print('solution with L10 {}, nb solutions {}'.format(L10, len(bearing_assembly_simulations)))
#                            break   
                    if len(bearing_assembly_simulations) > max_solutions:
                        break
                if len(bearing_assembly_simulations) > max_solutions:
                    break
            if len(bearing_assembly_simulations) > max_solutions:
                break
            
        self.bearing_assembly_simulations = [bearing_assembly_simulations[i] for i in npy.argsort(sort_bearing_assembly_simulations)]
        print('Number of solutions: {}'.format(len(self.bearing_assembly_simulations)))
        
    def ContinuousOptimize(self, bearing_assembly):
        
        bc_results = []
        for bearing_combination in bearing_assembly.bearing_combinations:
            li_bg_results = []
            for bearing in bearing_combination.bearings:
                li_bg_results.append(BearingSimulationResult())
            bc_results.append(BearingCombinationSimulationResult(li_bg_results))
        bearing_assembly_simulation_result = BearingAssemblySimulationResult(bc_results, 
                                                        self.loads, self.speeds, self.operating_times)
            
        l1 = bearing_assembly.bearing_combinations[0].B
        l2 = bearing_assembly.bearing_combinations[1].B
        pos1_min = self.axial_positions[0] + l1/2.
        pos1_max = self.axial_positions[0] + self.lengths[0] - l1/2.
        pos2_min = self.axial_positions[1] + l2/2.
        pos2_max = self.axial_positions[1] + self.lengths[1] - l2/2.
        pos1_moy = (pos1_min + pos1_max)/2.
        pos2_moy = (pos2_min + pos2_max)/2.
        def fun(x):
            
            obj = 0
            try:
                bearing_assembly.ShaftLoad([x[0], x[1]], bearing_assembly_simulation_result)
                L10 = bearing_assembly_simulation_result.L10
                obj += 1/(L10)**2
                return obj
            except:
                L10 = bearing_assembly_simulation_result.L10
                bearing_assembly.ShaftLoad([x[0], x[1]], bearing_assembly_simulation_result)
            
        def fineq(x):
            
            
            ineq = [0]

            return ineq
        Bound = [[pos1_min, pos1_max], [pos2_min, pos2_max]]
        sol_fun = math.inf
        for p1, p2 in product(Bound[0] + [pos1_moy],Bound[1] + [pos2_moy]):
#            cons = {'type': 'ineq','fun' : fineq}
            res = minimize(fun, [p1, p2], method='SLSQP', bounds=Bound)
            if fun(res.x) < sol_fun:
                sol_fun = fun(res.x)
                sol_x = res.x
                status = res.status
        for itera in range(0,5):
            # TODO: optimize this!
            x0 = (npy.array(Bound)[:,1]-npy.array(Bound)[:,0])*npy.random.random(2)+npy.array(Bound)[:,0]
#            cons = {'type': 'ineq','fun' : fineq}
            res = minimize(fun, x0, method='SLSQP', bounds=Bound)
            if fun(res.x) < sol_fun:
                sol_fun = fun(res.x)
                sol_x = res.x
                status = res.status
           
        if status >= 0:
            bearing_assembly.Update(sol_x, self.inner_diameters, self.axial_positions, 
                                self.outer_diameters, self.lengths)
            bearing_assembly.ShaftLoad(sol_x, bearing_assembly_simulation_result)
            
            bearing_assembly_simulation = BearingAssemblySimulation(bearing_assembly, bearing_assembly_simulation_result)
            return bearing_assembly_simulation
        else:
            return False

        
    def Dict(self, subobjects_id = {}, stringify_keys=True):
        """
        Export dictionary
        """
        d = {}
        for k,v in self.__dict__.items():
            tv = type(v)
            if tv == npy.int64:
                d[k] = int(v)
            elif tv==npy.float64:
                d[k]=round(float(v), 5)
            else:
                d[k] = v
        if self.bearing_assembly_simulations is not None:
            bar_dict = []
            for bar in self.bearing_assembly_simulations:
                if bar in subobjects_id:
                    bar_dict.append(subobjects_id[bar])
                else:                
                    bar_dict.append(bar.Dict())
        else:
            bar_dict = None
        d['bearing_assembly_simulations'] = bar_dict
        
        bc = []
        for bearing_classe in self.bearing_classes:
            bc.append(str(bearing_classe))
        d['bearing_classes'] = bc
                
        if stringify_keys:
            return StringifyDictKeys(d)

        return d
    
    @classmethod
    def DictToObject(cls, d):
        
        if 'bearing_assembly_simulations' in d:
            if d['bearing_assembly_simulations'] is None:
                li_bar = None
            else:
                li_bar = []
                for bar in d['bearing_assembly_simulations']:
                    li_bar.append(BearingAssemblySimulation.DictToObject(bar))
        else:
            li_bar = None
            
        bearing_classes = []
        for bearing_classe in d['bearing_classes']:
            bearing_classes.append(dict_bearing_classes[bearing_classe])
            
        obj = cls(loads = d['loads'], 
                  speeds = d['speeds'], 
                  operating_times = d['operating_times'],
                 inner_diameters = d['inner_diameters'],
                 axial_positions = d['axial_positions'],
                 outer_diameters = d['outer_diameters'],
                 lengths = d['lengths'],
                 linkage_types = d['linkage_types'],
                 mounting_types = d['mounting_types'],
                 number_bearings = d['number_bearings'],
                 bearing_classes = bearing_classes,
                 bearing_assembly_simulations = li_bar)
        return obj
