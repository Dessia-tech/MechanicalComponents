#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# cython: language_level=3
"""
Created on Wed Jul 10 15:49:05 2019

@author: steven
"""

import math
from itertools import product
from dectree import DecisionTree
from copy import deepcopy

import numpy as npy
npy.seterr(divide='raise', over='ignore', under='ignore', invalid='ignore')

from scipy.optimize import minimize, fsolve

from dessia_common import DessiaObject, dict_merge, Evolution, get_python_class_from_class_name
from typing import TypeVar, List



from mechanical_components.bearings import RadialBallBearing, AngularBallBearing, \
        SphericalBallBearing, \
        BearingAssembly, \
        BearingAssemblySimulationResult, BearingCombination,\
        BearingCombinationSimulationResult, BearingCombinationSimulation, BearingSimulationResult,\
        BearingAssemblySimulation, \
        ConceptualBearingCombination, \
        strength_bearing_classes, \
        BearingL10Error, CatalogSearchError, Linkage, Mounting, \
        CombinationMounting, SelectionLinkage
        
import mechanical_components
        

#import pandas

#from pandas.plotting import scatter_matrix
#import pkg_resources



class AxialPositionConvergenceError(Exception):
    def __init__(self):
        super().__init__('Fail in axial position optimization')


class BearingCombinationOptimizer(DessiaObject):
    def Configurations(self):
        
        configurations = []
        list_load_cases = []
        for axial_load in self.axial_loads:
            if axial_load > 0:
                list_load_cases.append(Mounting(right=True))
            elif axial_load < 0:
                list_load_cases.append(Mounting(left=True))
        
        for mounting_type in self.mounting_types:
            valid_load_case = True
            if len(list_load_cases) > 0:
                if mounting_type not in list_load_cases:
                    if not mounting_type.both:
                        valid_load_case = False
                elif len(list_load_cases) > 1:
                    valid_load_case = False
            if valid_load_case:
                for linkage in self.linkage_types:
                    configurations.append((mounting_type, linkage))
        return configurations
    
    def ConceptualBearingCombinations(self, max_bearings=2):
        bearing_combinations_possibilities = {}
        for mounting_type, linkage in self.Configurations():
            if not (mounting_type, linkage) in bearing_combinations_possibilities:
                DBC = mechanical_components.optimization.bearings.ConceptualBearingCombinationOptimizer(linkage, 
                                                    mounting_type, 
                                                    self.inner_diameter,
                                                    self.outer_diameter, 
                                                    self.length, 
                                                    bearing_classes = self.bearing_classes)
                bearing_combinations_possibilities[(mounting_type, linkage)] = \
                    DBC.ConceptualBearingCombinations(max_bearings)
        self.bearing_combinations_possibilities = bearing_combinations_possibilities
        
    def SelectBestBearingCombinations(self, first_bearing_possibilies, conceptual_bearing_combination, 
                                      L10_objective,max_speed =None):
        dt = DecisionTree()
        nb_bearings = len(conceptual_bearing_combination.bearing_classes)
        
        d_pre = 0
        bearing_possibilies = []
        list_next_bearings = [0]*nb_bearings
        current_Cr_pre = []
        while not dt.finished:
            valid = True
            bearings = []
            current_Cr = []
            for depth, node in enumerate(dt.current_node):
                if depth == 0:
                    bearings.append(first_bearing_possibilies[node])
                    d = bearings[-1].d
                    D = bearings[-1].D
                    B = bearings[-1].B
                    if (dt.current_depth == 1) and (d == d_pre):
                        valid = False
                        break
                    d_pre = d
                elif depth < nb_bearings:
                    next_bearings = list_next_bearings[depth]
                    bearings.append(next_bearings[node])
                    B += bearings[-1].B
                    
                current_Cr.append(bearings[-1].Cr)
                
                if B > self.length:
                    valid = False
                    break
                if bearings[-1].speed_limit:
                   
                    if max_speed>bearings[-1].speed_limit:
                        valid=False
                        
                        break
               
                if (bearings[-1].Cr < 20) or (d > D) or (d == 0) or (D == 0):
                    valid = False
                    
                    break
                
            if (dt.current_depth == nb_bearings) and valid:
                if current_Cr == current_Cr_pre:
                    valid = False
                current_Cr_pre = current_Cr
            
            # Testing
            if valid:
                # Counting possibilities
                if dt.current_depth == 0:
                    dt.SetCurrentNodeNumberPossibilities(len(first_bearing_possibilies))
                elif dt.current_depth < nb_bearings:
                    conceptual_bearing = conceptual_bearing_combination.bearing_classes[dt.current_depth]
                    list_next_bearings[dt.current_depth] = self.catalog.next_bearing_catalog(conceptual_bearing, d , D)
                    if list_next_bearings[dt.current_depth] is not False: 
                        dt.SetCurrentNodeNumberPossibilities(len(list_next_bearings[dt.current_depth]))
                    else:
                        dt.SetCurrentNodeNumberPossibilities(0)
                elif dt.current_depth == nb_bearings:
                    list_L10 = []
                    for bearing in bearings:
                        list_L10.append(bearing.estimate_base_life_time(Fr = self.radial_loads,
                                                                N = self.speeds, 
                                                                t = self.operating_times, Cr = bearing.Cr))
                    
                    L10 = BearingCombination.estimate_base_life_time(list_L10)
                    if L10 > L10_objective:
                        if bearings[0] not in bearing_possibilies:
                            bearing_possibilies.append(bearings[0])
                           
                    dt.SetCurrentNodeNumberPossibilities(0)
                else:
                    dt.SetCurrentNodeNumberPossibilities(0)
            else:
                dt.SetCurrentNodeNumberPossibilities(0)
                    
            dt.NextNode(valid)
        
        return bearing_possibilies
        
    def SelectBearingCombinations(self, bearing_combinations_possibility, L10_objective, max_speed =None):
        
        select_configurations = []
        li_quote = []
        bearing_combinations = []
        for conceptual_bearing_combination in bearing_combinations_possibility:
            quote = 0
            for bearing_classe in conceptual_bearing_combination.bearing_classes:
                quote += strength_bearing_classes[str(bearing_classe)]
            li_quote.append(quote)
            bearing_combinations.append(conceptual_bearing_combination)
        
        for conceptual_bearing_combination in [bearing_combinations[i] for i in npy.argsort(li_quote)]:
                    
            first_bearing_possibilies = self.catalog\
                .search_bearing_catalog(conceptual_bearing_combination.bearing_classes[0],
                                        self.inner_diameter, self.outer_diameter)
            
            if len(first_bearing_possibilies) == 0:
                continue
            
            # selection of the best first_bearing_left_possibilies
            bearing_possibilies = self.SelectBestBearingCombinations(first_bearing_possibilies,
                                               conceptual_bearing_combination,
                                               L10_objective,max_speed = max_speed)
            
            if len(bearing_possibilies) > 0:
                select_configurations.append([conceptual_bearing_combination, bearing_possibilies])

        return select_configurations
    
    def AnalyzeBearingCombinations(self, select_configurations, L10_objective,
                                  max_bearing_combinations=10):
        
        conceptual_bearing_combination, bearing_possibilies = select_configurations
       
        nb_bearings = len(conceptual_bearing_combination.bearing_classes)

        list_bearing_possibilities = [0]*nb_bearings
        list_Cr = []

        compt_same_configuration = 0
        compt_nb_eval_L10 = 0
        Cr_current_node_m = []
        
        d_pre = 0
        
        dt = DecisionTree()
        while not dt.finished:
            # Constructing BearingCombination
            valid = True
            bearings = []
            B = 0
            Cr_current_node = []
            
            for depth, node in enumerate(dt.current_node):
                
                if depth == 0:
                    
                    bearing = bearing_possibilies[node]
                    bearings.append(bearing)
                    d = bearing.d
                    if (d == d_pre) and (dt.current_depth == 1):
                        valid = False
                        
                        break
                    d_pre = d
                    D = bearing.D
                    B += bearing.B                    
                    
                elif depth < nb_bearings:
                    bearing_classe = conceptual_bearing_combination.bearing_classes[depth]
                    bearing_possibilities = list_bearing_possibilities[depth]
                    
                    bearings.append(bearing_possibilities[node])
                    B += bearing.B
                    
                if (d > D):
                    valid = False
                  
                    break
                if d == 0:
                    
                    valid = False
                  
                    break
                if D == 0:
                    valid = False
               
                    break
                
                if B > self.length:
                    
                    valid = False
                  
                    break
                
                Cr_current_node.append(bearings[-1].Cr)
                
#                if len(Cr_current_node_m) > 0:
#                    if Cr_current_node[-1] == Cr_current_node_m[depth]:
#                        valid = False
#                        break
                
                # a supprimer suite MAJ base rlts
                if Cr_current_node[-1] < 20:
                    valid = False
                    
                    break
                
                if len(list_Cr) > 0:
                    if Cr_current_node[-1] < min([cr[depth] for cr in list_Cr]):
                        
                        valid = False
                       
                        break
                
                    
            if (dt.current_depth == nb_bearings) and valid:
               
                if Cr_current_node_m == Cr_current_node:
                    valid = False
                    
                Cr_current_node_m = Cr_current_node

            if (dt.current_depth == nb_bearings) and valid:
                
                bc = conceptual_bearing_combination.bearing_combination(bearings)
                li_bg_results = []
                for bearing in bearings:
                    li_bg_results.append(BearingSimulationResult())
                bearing_combination_simulation_result = BearingCombinationSimulationResult(li_bg_results, self.axial_loads,
                                                                self.radial_loads, self.speeds, self.operating_times)
                # bc.plot()
                
                check = bc.base_life_time(bearing_combination_simulation_result)
                
                if check == False:
                    break
                
                L10 = bearing_combination_simulation_result.L10
                
                
                if L10 is False:
                    break
                if L10 < L10_objective:
                   
                    valid = False
                
                    compt_nb_eval_L10 += 1
                    if compt_nb_eval_L10 > 5:
                        break
                    
                    def funct(alpha):
                        bgs = [deepcopy(bg) for bg in bearings]
                        for bg in bgs:
                            bg.Cr = bg.Cr*alpha[0]
                        bc = conceptual_bearing_combination.bearing_combination(bgs)
                        li_bg_results = []
                        for bearing in bearings:
                            li_bg_results.append(BearingSimulationResult())
                        bearing_combination_simulation_result = BearingCombinationSimulationResult(li_bg_results, self.axial_loads,
                                                                        self.radial_loads, self.speeds, self.operating_times)
        
                        check = bc.base_life_time(bearing_combination_simulation_result)
                        if not check:
                            return False
                        
                        L10 = bearing_combination_simulation_result.L10
                        
                        if L10 is not False:
                            return (L10 - L10_objective)**2
                        else:
                            return False
                    
                    valid_fsolve = False
                    for i in range(2):
                        cond_init = npy.random.random(1)*3
                        sol = fsolve(funct, cond_init[0])[0]
                        if funct([sol]) < 1e-4:
                            coefficient_Cr = sol
                            valid_fsolve = True
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
                    
                    dt.SetCurrentNodeNumberPossibilities(len(bearing_possibilies))
                    list_bearing_possibilities[dt.current_depth] = bearing_possibilies
                elif dt.current_depth < nb_bearings:
                    
                    bearing_classe = conceptual_bearing_combination.bearing_classes[dt.current_depth]
                    bearing = bearing_possibilies[dt.current_node[0]]
                    d = bearing.d
                    D = bearing.D
                    try:
                        
                        bearing_possibilities = self.catalog.next_bearing_catalog(bearing_classe, d , D)
                        dt.SetCurrentNodeNumberPossibilities(len(bearing_possibilities))
                        list_bearing_possibilities[dt.current_depth] = bearing_possibilities
                       
                    except CatalogSearchError:
                        valid = False
                       
                        dt.SetCurrentNodeNumberPossibilities(0)
                elif dt.current_depth == nb_bearings:
                    
                    dt.SetCurrentNodeNumberPossibilities(0) 
                    bc = conceptual_bearing_combination.bearing_combination(bearings)
                    compt_same_configuration += 1
                    
                    yield bc
                    
                    if compt_same_configuration > 0:
                        break
                else:
                    dt.SetCurrentNodeNumberPossibilities(0)      
            else:
                dt.SetCurrentNodeNumberPossibilities(0)
            
            dt.NextNode(valid)
        
    def optimize(self, max_solutions:int=10,max_speed=None)->None:
        
        L10_objective = 0
        for speed, time in zip(self.speeds, self.operating_times):
            L10_objective += speed/(2*math.pi)*time
        L10_objective = L10_objective/1e6
        
#        pos_min = -self.length/2.
#        pos_max = self.length/2.

        bearing_combination_simulations = []
        sort_bearing_combination_simulations = []
        for max_bearings in self.number_bearings:
            self.ConceptualBearingCombinations(max_bearings = (max_bearings))
            for (mounting_type, linkage), bearing_combinations_possibility in self.bearing_combinations_possibilities.items():
                bearing_combination_configurations = self.SelectBearingCombinations(bearing_combinations_possibility, 
                                                                                    L10_objective = L10_objective,max_speed=max_speed)
                try:
                    li_bearing_assembly_configurations.extend(bearing_combination_configurations)
                except NameError:
                    li_bearing_assembly_configurations = bearing_combination_configurations
                    
            for bearing_combination_configurations in li_bearing_assembly_configurations:
                bearing_combinations = self.AnalyzeBearingCombinations(bearing_combination_configurations, 
                                                                     L10_objective = L10_objective,
                                                                     max_bearing_combinations = max_solutions)
                
                
                

                for i_bearing_combination, bearing_combination in enumerate(bearing_combinations):
                    
                    li_bg_results = []
                    for bearing in bearing_combination.bearings:
                        li_bg_results.append(BearingSimulationResult())
                    bearing_combination_simulation_result = BearingCombinationSimulationResult(li_bg_results, self.axial_loads,
                                                                    self.radial_loads, self.speeds, self.operating_times)
                    check = bearing_combination.base_life_time(bearing_combination_simulation_result)
                  
                    
                    if check != False:
                        
                        bearing_combination.update(-self.length/2., self.inner_diameter, self.outer_diameter, self.length)
                        bearing_combination_simulation = BearingCombinationSimulation(bearing_combination, bearing_combination_simulation_result)
                        L10 = bearing_combination_simulation_result.L10
                        if L10 >= L10_objective:
                            bearing_combination_simulations.append(bearing_combination_simulation)
                            sort_bearing_combination_simulations.append(L10)
                            
#                            break   
                    if len(bearing_combination_simulations) > max_solutions:
                        break
                if len(bearing_combination_simulations) > max_solutions:
                    break
            if len(bearing_combination_simulations) > max_solutions:
                break
            
        self.bearing_combination_simulations = [bearing_combination_simulations[i] for i in npy.argsort(sort_bearing_combination_simulations)]

        

class ConceptualBearingCombinationOptimizer(DessiaObject):
    
    def CheckLinkage(self, bearings):
        check = False
        if (len(bearings) > 1) and (self.linkage.cylindric_joint):
            check = True
        if (len(bearings) == 1) and (self.linkage.cylindric_joint):
            if bearings[0] not in [RadialBallBearing, 
                                   AngularBallBearing,
                                   SphericalBallBearing]:
                check = True
        if (len(bearings) == 1) and (self.linkage.ball_joint):
            if bearings[0] in [RadialBallBearing, 
                               AngularBallBearing,
                               SphericalBallBearing]:
                check = True
        return check
            
        
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

                    if type(self.bearing_classes[node]) == str:
                        bearing = get_python_class_from_class_name(self.bearing_classes[node])
                        bearings.append(bearing)
                    else:
                        bearing = self.bearing_classes[node]
                        bearings.append(bearing)

                else:
                    
                    directions.append([-1, 1][node])

            if dt.current_depth // 2 == max_bearings:
                # Instanciating 
                valid = self.CheckLinkage(bearings)
                
                if valid:
                    cbc = ConceptualBearingCombination(bearings, directions, self.mounting)
                    valid = cbc.check_kinematic()  
                   
            
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

class BearingAssemblyOptimizer(DessiaObject):

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
        

        for mounting_type in self.mounting_types:
            mounting_type_left, mounting_type_right = mounting_type.mountings
            valid_load_case = True
            if not mounting_type_left.both and not mounting_type_right.both:
                for load_case in list_load_cases:
                 
                    mouting_verification=[]
                    for mouting in [mounting_type_left,mounting_type_right]:
                        if mouting.right:
                            mouting_verification.append('right')
                        if mouting.left:
                            mouting_verification.append('left')
                    
                    
                    if load_case not in mouting_verification:
                        valid_load_case = False
            
            
            if mounting_type_left.free and mounting_type_right.free:
                valid_load_case = False
            # if 'both' in set((mounting_type_left, mounting_type_right)):
            
            if mounting_type_left.both or mounting_type_right.both:
                if not mounting_type_left.free and not mounting_type_right.free:
                # if 'free' not in set((mounting_type_left, mounting_type_right)):
                    valid_load_case = False
                    
            
            if valid_load_case:
                for linkage_left, linkage_right in product(*[lt.linkages for lt in self.linkage_types]):
                    configurations.append(((mounting_type_left, linkage_left),
                                           (mounting_type_right, linkage_right)))
        return configurations

    def ConceptualBearingCombinations(self, max_bearings=[2, 2]):
        bearing_combinations_possibilities = {}
        
        for left, right in self.Configurations():
            
            if not (left, right) in bearing_combinations_possibilities:
                
                
                DBC_l = mechanical_components.optimization.bearings.ConceptualBearingCombinationOptimizer(left[1], 
                                                    left[0], 
                                                    self.inner_diameters[0],
                                                    self.outer_diameters[0], 
                                                    self.lengths[0], 
                                                    bearing_classes = self.bearing_classes
                                                    )
                DBC_r =  mechanical_components.optimization.bearings.ConceptualBearingCombinationOptimizer(right[1], 
                                                      right[0], 
                                                      self.inner_diameters[1],
                                                      self.outer_diameters[1], 
                                                      self.lengths[1], 
                                                      bearing_classes = self.bearing_classes
                                                     )
                
                combination_left = DBC_l.ConceptualBearingCombinations(max_bearings[0])
                combination_right = DBC_r.ConceptualBearingCombinations(max_bearings[1])
           
                if (len(combination_left) > 0) and (len(combination_right) > 0):
                    bearing_combinations_possibilities[(left, right)] = (combination_left, combination_right)

        return bearing_combinations_possibilities
        
    def SelectBestBearingCombinations(self, first_bearing_possibilies, conceptual_bearing_combination, 
                                      radial_elementary_load, L10_objective, length):
        dt = DecisionTree()
        nb_bearings = len(conceptual_bearing_combination.bearing_classes)
        
        best_L10 = 0
        d_pre = 0
        bearing_possibilies = []
        list_next_bearings = [0]*nb_bearings
#        current_Cr_pre = []
        while not dt.finished:
            valid = True
            bearings = []
            current_Cr = []
            for depth, node in enumerate(dt.current_node):
                if depth == 0:
                    bearings.append(first_bearing_possibilies[node])
                    d = bearings[-1].d
                    D = bearings[-1].D
                    B = bearings[-1].B
                    if (dt.current_depth == 1) and (d == d_pre):
                        valid = False
                        break
#                    d_pre = d
                elif depth < nb_bearings:
                    next_bearings = list_next_bearings[depth]
                    bearings.append(next_bearings[node])
                    B += bearings[-1].B
                    
                current_Cr.append(bearings[-1].Cr)
                
                if B > length:
                    valid = False
                    break
                
                if (bearings[-1].Cr < 1) or (d > D) or (d == 0) or (D == 0):
                    valid = False
                    break
                
#            if (dt.current_depth == nb_bearings) and valid:
#                if current_Cr == current_Cr_pre:
#                    valid = False
#                current_Cr_pre = current_Cr
            
            # Testing
            if valid:
                # Counting possibilities
                if dt.current_depth == 0:
                    dt.SetCurrentNodeNumberPossibilities(len(first_bearing_possibilies))
                elif dt.current_depth < nb_bearings:
                    conceptual_bearing = conceptual_bearing_combination.bearing_classes[dt.current_depth]
                    list_next_bearings[dt.current_depth] = self.catalog.next_bearing_catalog(conceptual_bearing, d , D)
                    if list_next_bearings[dt.current_depth] is not False: 
                        dt.SetCurrentNodeNumberPossibilities(len(list_next_bearings[dt.current_depth]))
                    else:
                        dt.SetCurrentNodeNumberPossibilities(0)
                elif dt.current_depth == nb_bearings:
                    list_L10 = []
                    for bearing in bearings:
                        list_L10.append(bearing.estimate_base_life_time(Fr = radial_elementary_load,
                                                                N = self.speeds, 
                                                                t = self.operating_times, Cr = bearing.Cr))
                        L10 = BearingAssembly.estimate_base_life_time(list_L10)
    
                        if L10 > L10_objective:
                            best_L10 = max(best_L10, L10)

                            if bearings[0] not in bearing_possibilies:
                                bearing_possibilies.append(bearings[0])

                    dt.SetCurrentNodeNumberPossibilities(0)
                else:
                    dt.SetCurrentNodeNumberPossibilities(0)
            else:
                dt.SetCurrentNodeNumberPossibilities(0)
                    
            dt.NextNode(valid)
        
        return best_L10, bearing_possibilies
        
    def SelectBearingCombinations(self, bearing_combinations_possibility, radial_load, 
                                  L10_objective):
        
        radial_load_left, radial_load_right = radial_load
        
        compt_continue = 0
        select_configurations = []
        select_configurations_L10 = []
        
        li_quote = []
        bearing_combinations = []
        for conceptual_bearing_combination_left, conceptual_bearing_combination_right in \
                product(*bearing_combinations_possibility):
            quote = 0
            for bearing_classe in conceptual_bearing_combination_left.bearing_classes:
                quote += strength_bearing_classes[str(bearing_classe)]
            for bearing_classe in conceptual_bearing_combination_right.bearing_classes:
                quote += strength_bearing_classes[str(bearing_classe)]
            li_quote.append(quote)
            bearing_combinations.append([conceptual_bearing_combination_left, conceptual_bearing_combination_right])
        
        for conceptual_bearing_combination_left, conceptual_bearing_combination_right in \
                [bearing_combinations[i] for i in npy.argsort(li_quote)]:
                    
            first_bearing_left_possibilies = self.catalog\
                .search_bearing_catalog(conceptual_bearing_combination_left.bearing_classes[0],
                                        self.inner_diameters[0], self.outer_diameters[0])
            first_bearing_right_possibilies = self.catalog\
                .search_bearing_catalog(conceptual_bearing_combination_right.bearing_classes[0],
                                        self.inner_diameters[1], self.outer_diameters[1])
            nb_bearings_left = len(conceptual_bearing_combination_left.bearing_classes)
            nb_bearings_right = len(conceptual_bearing_combination_right.bearing_classes)
            nb_bearings = nb_bearings_left + nb_bearings_right# TODO why unussed?
            
            if len(first_bearing_left_possibilies) == 0 or len(first_bearing_right_possibilies) == 0:
                continue
            
            # selection of the best first_bearing_left_possibilies
            best_L10_left, bearing_left_possibilies = self.SelectBestBearingCombinations(first_bearing_left_possibilies,
                                               conceptual_bearing_combination_left,
                                               [r/(1.*nb_bearings_left) for r in radial_load_left],
                                               L10_objective, self.lengths[0])
#            bearing_left_possibilies = first_bearing_left_possibilies
                    
            # selection of the best first_bearing_right_possibilies
            best_L10_right, bearing_right_possibilies = self.SelectBestBearingCombinations(first_bearing_right_possibilies,
                                               conceptual_bearing_combination_right,
                                               [r/(1.*nb_bearings_right) for r in radial_load_right],
                                               L10_objective, self.lengths[1])
            
            if (best_L10_left != 0) and (best_L10_right != 0):
                L10 = BearingAssembly.estimate_base_life_time([best_L10_left, best_L10_right])
            else:
                L10 = 0

            if L10 < L10_objective:
                compt_continue += 1

#                if compt_continue > 30:
#                    break
                continue
            else:
                select_configurations.append([conceptual_bearing_combination_left, conceptual_bearing_combination_right,
                                              bearing_left_possibilies, bearing_right_possibilies])
                select_configurations_L10.append(L10)
#            if len(select_configurations) > 10:
#                break
            
        return select_configurations, select_configurations_L10
        
    def AnalyzeBearingCombinations(self, select_configurations, L10_objective,
                                  max_bearing_assemblies=10):
        
        conceptual_bearing_combination_left, conceptual_bearing_combination_right,\
            bearing_left_possibilies, bearing_right_possibilies = select_configurations

        nb_bearings_left = len(conceptual_bearing_combination_left.bearing_classes)
        nb_bearings_right = len(conceptual_bearing_combination_right.bearing_classes)
        nb_bearings = nb_bearings_left + nb_bearings_right

        list_bearing_possibilities = [0]*nb_bearings
        list_Cr = []

        compt_same_configuration = 0
        compt_nb_eval_L10 = 0
#        Cr_current_node_m = []
#        
#        d_left_pre = 0
#        d_right_pre = 0
        
        dt = DecisionTree()

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
#                    if (d == d_left_pre) and (dt.current_depth == 1):
#                        valid = False
#                        break
#                    d_left_pre = d
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
#                    if (d == d_right_pre) and (dt.current_depth == 1):
#                        valid = False
#                        break
#                    d_right_pre = d
                    D = bearing.D
                    B_right += bearing.B
#                        if d_right < d_left:
#                            valid = False
#                            break
                    
                elif depth < nb_bearings:
                    bearing_classe = conceptual_bearing_combination_right.bearing_classes[depth - nb_bearings_left]
                    bearing_possibilities = list_bearing_possibilities[depth]
                    bearings.append(bearing_possibilities[node])
                    B_right += bearing.B
                    
                if (d > D):
                    valid = False
                    break
                if d == 0:
                    valid = False
                    break
                if D == 0:
                    valid = False
                    break
                
                if (B_left > self.lengths[0]) or (B_right > self.lengths[1]):
                    valid = False
                    break
                
                Cr_current_node.append(bearings[-1].Cr)
#                if len(Cr_current_node_m) > 0:
#                    if Cr_current_node[-1] == Cr_current_node_m[depth]:
#                        valid = False
#                        break
                
                # a supprimer suite MAJ base rlts
#                if Cr_current_node[-1] < 20:
#                    valid = False
#                    break
                
#                if len(list_Cr) > 0:
#                    if Cr_current_node[-1] < min([cr[depth] for cr in list_Cr]):
#                        valid = False
#                        break
                    
            if (dt.current_depth == nb_bearings) and valid:
                    
                if (B_left > self.lengths[0]) or (B_right > self.lengths[1]):
                    valid = False
                    
#                if Cr_current_node_m == Cr_current_node:
#                    valid = False
#                Cr_current_node_m = Cr_current_node
            
            if (dt.current_depth == nb_bearings) and valid:
                
                bc_left = conceptual_bearing_combination_left.bearing_combination(bearings[0: nb_bearings_left])
                bc_right = conceptual_bearing_combination_right.bearing_combination(bearings[nb_bearings_left:])
                bearing_assembly = BearingAssembly([bc_left, bc_right])
                bc_results = []
                for bearing_combination in bearing_assembly.bearing_combinations:
                    li_bg_results = []
                    for bearing in bearing_combination.bearings:
                        li_bg_results.append(BearingSimulationResult())
                        
                    bc_results.append(BearingCombinationSimulationResult(li_bg_results))
                bearing_assembly_simulation_result = BearingAssemblySimulationResult(bc_results, 
                                                                self.loads, self.speeds, self.operating_times)
                
#                pos1_min = self.axial_positions[0]
#                pos1_max = self.axial_positions[0] + self.lengths[0]
#                pos2_min = self.axial_positions[1]
#                pos2_max = self.axial_positions[1] + self.lengths[1]
                
                pos1_min = self.axial_positions[0] + B_left/2.
                pos1_max = max(pos1_min, self.axial_positions[0] + self.lengths[0] - B_left/2.)
                pos2_min = self.axial_positions[1] + B_right/2.
                pos2_max = max(pos2_min, self.axial_positions[1] + self.lengths[1] - B_right/2.)
                pos1_moy = (pos1_min + pos1_max)/2.
                pos2_moy = (pos2_min + pos2_max)/2.
                
                L10 = 0
                for pos1, pos2 in product([pos1_min, pos1_moy, pos1_max], [pos2_min, pos2_moy, pos2_max]):
                    try:
                        bearing_assembly.shaft_load([pos1, pos2], 
                                                    bearing_assembly_simulation_result)
                        
                        L10 = max(L10, bearing_assembly_simulation_result.L10)
                    except BearingL10Error:
                        pass
                
                if L10 == 0:
                    break
                
                if L10 < L10_objective:

                    valid = False
                
                    compt_nb_eval_L10 += 1
                    if compt_nb_eval_L10 > 5:
                        break
                    
                    def funct(alpha):
                        bgs = [deepcopy(bg) for bg in bearings]
                        for bg in bgs:
                            bg.Cr = bg.Cr*alpha[0]
                        bc_left = conceptual_bearing_combination_left.bearing_combination(bgs[0: nb_bearings_left])
                        bc_right = conceptual_bearing_combination_right.bearing_combination(bgs[nb_bearings_left:])
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
                        try:
                            bearing_assembly.shaft_load([(pos1_min + pos1_max)/2., (pos2_min + pos2_max)/2.], 
                                                        bearing_assembly_simulation_result)
                            L10 = bearing_assembly_simulation_result.L10
                            return (L10 - L10_objective)**2
                        except BearingL10Error:
                            return False
                    
                    valid_fsolve = False
                    for i in range(5):
                        cond_init = npy.random.random(1)*3
                        sol = fsolve(funct, cond_init[0])[0]
                        if funct([sol]) < 1e-4:
                            coefficient_Cr = sol
                            valid_fsolve = True

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
                        bearing_possibilities = self.catalog.next_bearing_catalog(bearing_classe, d , D)
                        
                        
                        dt.SetCurrentNodeNumberPossibilities(len(bearing_possibilities))
                        list_bearing_possibilities[dt.current_depth] = bearing_possibilities
                    except CatalogSearchError:
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
                        bearing_possibilities = self.catalog.next_bearing_catalog(bearing_classe, d , D)
                        dt.SetCurrentNodeNumberPossibilities(len(bearing_possibilities))
                        list_bearing_possibilities[dt.current_depth] = bearing_possibilities
                    except CatalogSearchError:
                        valid = False
                        
                        dt.SetCurrentNodeNumberPossibilities(0)
                
                elif dt.current_depth == nb_bearings:
                  
                    dt.SetCurrentNodeNumberPossibilities(0) 
                    bc_left = conceptual_bearing_combination_left.bearing_combination(bearings[0: nb_bearings_left])
                    bc_right = conceptual_bearing_combination_right.bearing_combination(bearings[nb_bearings_left:])
                    ba = BearingAssembly([bc_left, bc_right])
                    compt_same_configuration += 1
                    yield ba
                    
                    if compt_same_configuration > max_bearing_assemblies - 1:
                        break
                else:
                    
                    dt.SetCurrentNodeNumberPossibilities(0)      
            else:
                dt.SetCurrentNodeNumberPossibilities(0)
            
            dt.NextNode(valid)
    
    def OptimizeGeneric(self, max_solutions=10, nb_solutions_family=10,
                        progress_callback=lambda x:0,
                        verbose=False):
       
        L10_objective = 0
        
        for speed, time in zip(self.speeds, self.operating_times):
            L10_objective += speed/(2*math.pi)*time
        L10_objective = L10_objective/1e6
        if verbose:
            print('the L10 objective is {}'.format(L10_objective))
        
        pos1_min = self.axial_positions[0]
        pos1_max = self.axial_positions[0] + self.lengths[0]
        pos2_min = self.axial_positions[1]
        pos2_max = self.axial_positions[1] + self.lengths[1]
        pos1_moy = (pos1_min + pos1_max)/2.
        pos2_moy = (pos2_min + pos2_max)/2.
        
        radial_load_left = []
        radial_load_right = []
        for load in self.loads:
            radial_load_left_temp = []
            radial_load_right_temp = []
            for pos1, pos2 in product([pos1_min, pos1_moy, pos1_max], [pos2_min, pos2_moy, pos2_max]):
                load_simul = BearingAssembly.quick_shaft_load((pos1, pos2), [load])
                for rl in load_simul:
                    radial_load_left_temp.append((rl[0]**2 + (0)**2)**0.5)
                    radial_load_right_temp.append((rl[1]**2 + (0)**2)**0.5)
            radial_load_left.append(min(radial_load_left_temp))
            radial_load_right.append(min(radial_load_right_temp))
            
#        nb_solutions_family = max_solutions_family[0]
            
        bearing_assembly_generic = []
        combination_number_bearings = list(product(*self.number_bearings))
        combination_number_bearings = [combination_number_bearings[j] for j in npy.argsort([sum(i) for i in combination_number_bearings])]
        ncnb = float(len(combination_number_bearings))
       
        for icnb, (max_bearings_left, max_bearings_right) in enumerate(combination_number_bearings):
            if verbose:
                print('number of bearings analyzed: {} left and {} right'.format(max_bearings_left, max_bearings_right))
            progress_callback(icnb/ncnb)
            
            li_bearing_assembly_configurations = []
            li_bearing_assembly_L10 = []
            
            bearing_combinations_possibilities = self.ConceptualBearingCombinations(max_bearings = (max_bearings_left, max_bearings_right))
            
            if len(bearing_combinations_possibilities) == 0:
                continue
            for (left, right), bearing_combinations_possibility in bearing_combinations_possibilities.items():
                bearing_assembly_configurations, bearing_assembly_L10 = self.SelectBearingCombinations(bearing_combinations_possibility, 
                                                                    radial_load = (radial_load_left, radial_load_right), 
                                                                    L10_objective = L10_objective)

                li_bearing_assembly_configurations.extend(bearing_assembly_configurations)
                li_bearing_assembly_L10.extend(bearing_assembly_L10)
                
            bearing_assembly_configurations_sort = li_bearing_assembly_configurations
#            bearing_assembly_configurations_sort = [li_bearing_assembly_configurations[i] for i in npy.argsort(li_bearing_assembly_L10)[::-1]]
            if verbose:
                print('number bearing assemblies configurations {}'.format(len(bearing_assembly_configurations_sort)))
            
            for bearing_assembly_configurations in bearing_assembly_configurations_sort:
                bearing_assemblies = self.AnalyzeBearingCombinations(bearing_assembly_configurations, 
                                                                     L10_objective = L10_objective,
                                                                     max_bearing_assemblies=nb_solutions_family)
                cas_bearing_assembly_simulations = []
                
                for i_bearing_assembly, bearing_assembly in enumerate(bearing_assemblies):
                    cas_bearing_assembly_simulations.append(bearing_assembly)
                if cas_bearing_assembly_simulations != []:
                    bearing_assembly_generic.append(cas_bearing_assembly_simulations)
                if verbose:
                   
                    print('size solutions {}'.format(len(bearing_assembly_generic)))
        return bearing_assembly_generic
                
    def optimize(self, max_solutions:int=10, 
                 nb_solutions_family:int=10,
                 verbose:bool=False)->None:
#        progress_callback=lambda x:0,
        
        bearing_assembly_generic = self.OptimizeGeneric(max_solutions,verbose=True)
     
        L10_objective = 0
        for speed, time in zip(self.speeds, self.operating_times):
            L10_objective += speed/(2*math.pi)*time
        L10_objective = L10_objective/1e6
        
        #Cost optimization
        #Initial sort before continuous optim
        list_cost = []
        for cas_bearing_assembly_simulations in bearing_assembly_generic:
            li_cost = [ba.cost for ba in cas_bearing_assembly_simulations]
            list_cost.append(min(li_cost))
        bearing_assembly_simulations_sort = [bearing_assembly_generic[i] for i in npy.argsort(list_cost)]
       
        bearing_assembly_simulations = []
        sort_bearing_assembly_simulations = []
       
        for bearing_assemblies in bearing_assembly_simulations_sort:
            list_cost_temp = [ba.cost for ba in bearing_assemblies]
            bearing_assemblies_temp = [bearing_assemblies[i] for i in npy.argsort(list_cost_temp)]
            cas_bearing_assembly_simulations = []
            for i_bearing_assembly, bearing_assembly in enumerate(bearing_assemblies_temp):
                try:
                    bearing_assembly_simulation = self.ContinuousOptimization(bearing_assembly, L10_objective)
                    L10 = bearing_assembly_simulation.bearing_assembly_simulation_result.L10
                    mass = bearing_assembly_simulation.bearing_assembly.mass
                    cost = bearing_assembly_simulation.bearing_assembly.cost
                    if L10 >= L10_objective:
                        cas_bearing_assembly_simulations.append(bearing_assembly_simulation)
                        sort_bearing_assembly_simulations.append(mass)
                        if verbose:
                            print('solution with L10 {}, nb solutions {}, cost {}, mass {}'.format(L10, len(cas_bearing_assembly_simulations), cost, mass))
                except AxialPositionConvergenceError:
                    pass
                if len(bearing_assembly_simulations) > nb_solutions_family:
                    break
            if cas_bearing_assembly_simulations!= []:
                bearing_assembly_simulations.append(cas_bearing_assembly_simulations)
            if len(bearing_assembly_simulations) > nb_solutions_family:
                break
        #Final sort after optimization
        list_cost = []
        for bas in bearing_assembly_simulations:
            li_cost = [ba.bearing_assembly.cost for ba in bas]
            list_cost.append(min(li_cost))
        bas_sort = [bearing_assembly_simulations[i] for i in npy.argsort(list_cost)]
        list_cost_best = []
        for bas in bas_sort[0]:
            list_cost_best.append(bas.bearing_assembly.cost)
            
        self.bearing_assembly_simulations = [bas_sort[0][i] for i in npy.argsort(list_cost_best)]
        self.bearing_assembly_simulations.extend([bas[0] for bas in bas_sort[1:]])
                
        #Mass optimization
        #Initial sort before continuous optim
        list_mass = []
        for cas_bearing_assembly_simulations in bearing_assembly_generic:
            li_mass = [ba.mass for ba in cas_bearing_assembly_simulations]
            list_mass.append(min(li_mass))
        bearing_assembly_simulations_sort = [bearing_assembly_generic[i] for i in npy.argsort(list_mass)]
            

        bearing_assembly_simulations = []
        sort_bearing_assembly_simulations = []
        for bearing_assemblies in bearing_assembly_simulations_sort:
            list_mass_temp = [ba.mass for ba in bearing_assemblies]
            bearing_assemblies_temp = [bearing_assemblies[i] for i in npy.argsort(list_mass_temp)]
            cas_bearing_assembly_simulations = []
            for i_bearing_assembly, bearing_assembly in enumerate(bearing_assemblies_temp):
                try:
                    bearing_assembly_simulation = self.ContinuousOptimization(bearing_assembly, L10_objective)
                    L10 = bearing_assembly_simulation.bearing_assembly_simulation_result.L10
                    mass = bearing_assembly_simulation.bearing_assembly.mass
                    cost = bearing_assembly_simulation.bearing_assembly.cost
#                    print('L10 {}, nb solutions {}, cost {}, mass {}'.format(L10, len(cas_bearing_assembly_simulations), cost, mass))
                    if L10 >= L10_objective:
                        cas_bearing_assembly_simulations.append(bearing_assembly_simulation)
                        sort_bearing_assembly_simulations.append(mass)
#                        print('solution with L10 {}, nb solutions {}, cost {}, mass {}'.format(L10, len(cas_bearing_assembly_simulations), cost, mass))
                except AxialPositionConvergenceError:
                    pass
                if len(bearing_assembly_simulations) > nb_solutions_family:
                    break
            if cas_bearing_assembly_simulations!= []:
                bearing_assembly_simulations.append(cas_bearing_assembly_simulations)
            if len(bearing_assembly_simulations) > nb_solutions_family:
                break
        #Final sort after optimization
        list_mass = []
        for bas in bearing_assembly_simulations:
            li_mass = [ba.bearing_assembly.mass for ba in bas]
            list_mass.append(min(li_mass))
        bas_sort = [bearing_assembly_simulations[i] for i in npy.argsort(list_mass)]
        list_mass_best = []
        for bas in bas_sort[0]:
            list_mass_best.append(bas.bearing_assembly.mass)
            
        self.bearing_assembly_simulations.extend([bas_sort[0][i] for i in npy.argsort(list_mass_best)])
        self.bearing_assembly_simulations.extend([bas[0] for bas in bas_sort[1:]])
        
        
    def ContinuousOptimization(self, bearing_assembly, L10_objective):
        
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
        pos1_max = max(pos1_min, self.axial_positions[0] + self.lengths[0] - l1/2.)
        pos2_min = self.axial_positions[1] + l2/2.
        pos2_max = max(pos2_min, self.axial_positions[1] + self.lengths[1] - l2/2.)
        pos1_moy = (pos1_min + pos1_max)/2.
        pos2_moy = (pos2_min + pos2_max)/2.
        Bound = [[pos1_min, pos1_max], [pos2_min, pos2_max]]
        def fun(x,bound=Bound):
            obj = 0
            try:
                bearing_assembly.shaft_load([x[0], x[1]], bearing_assembly_simulation_result)
                L10 = bearing_assembly_simulation_result.L10
                obj += 1/(L10)**2
               
                if x[0]<bound[0][0]:
                    x[0]=bound[0][0]
                if x[0]>bound[0][1]:
                    x[0]=bound[0][1]
                if x[1]<bound[1][0]:
                    x[1]=bound[1][0]
                if x[1]>bound[1][1]:
                    x[1]=bound[1][1]
                return obj
            except BearingL10Error:
                return 1e6
            
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
                sol_x = [float(x) for x in res.x]
                status = res.status
        for itera in range(0, 5):
            x0 = (npy.array(Bound)[:,1]-npy.array(Bound)[:,0])*npy.random.random(2)+npy.array(Bound)[:,0]
#            cons = {'type': 'ineq','fun' : fineq}
            res = minimize(fun, x0, method='SLSQP', bounds=Bound)
            if fun(res.x) < sol_fun:
                sol_fun = fun(res.x)
                sol_x = [float(x) for x in res.x]
                status = res.status
           
        if status >= 0:
            bearing_assembly.update(sol_x, self.inner_diameters, self.axial_positions, 
                                self.outer_diameters, self.lengths)
            bearing_assembly.shaft_load(sol_x, bearing_assembly_simulation_result)
            
            bearing_assembly_simulation = BearingAssemblySimulation(bearing_assembly, bearing_assembly_simulation_result)
            return bearing_assembly_simulation
        else:
            raise AxialPositionConvergenceError()