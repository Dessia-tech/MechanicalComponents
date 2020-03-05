#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# cython: language_level=3
"""
Created on Fri Aug 17 02:14:21 2018

@author: steven
"""

try:
    _open_source = True
    import mechanical_components.optimization.bearings_protected as protected_module
except (ModuleNotFoundError, ImportError) as _:
    _open_source = False

from mechanical_components.bearings import BearingCombination, \
        BearingCatalog,\
        BearingCombinationSimulationResult, BearingSimulationResult,\
        BearingAssemblySimulation, BearingCombinationSimulation, \
        bearing_classes, dict_bearing_classes, \
        strength_bearing_classes, RadialBearing, \
        BearingL10Error, CatalogSearchError, Linkage, Mounting
        
import numpy as npy

from dessia_common import DessiaObject, dict_merge, Evolution
from typing import TypeVar, List

npy.seterr(divide='raise', over='ignore', under='ignore', invalid='ignore')

#from scipy.optimize import fsolve
#from copy import deepcopy
from itertools import product

#from dectree import DecisionTree

#import math
#
from mechanical_components.tools import StringifyDictKeys


class BearingCombinationOptimizer(protected_module.BearingCombinationOptimizer if _open_source==True else DessiaObject):
    
    def __init__(self, radial_loads:List[float], axial_loads:List[float], 
                 speeds:List[float], operating_times:List[float],
                 inner_diameter:float, outer_diameter:float,
                 length:float,
                 linkage_types:List[str]=['ball_joint', 'cylindric_joint'],
                 mounting_types:List[str]=['both', 'right', 'left', 'free'],
                 number_bearings:List[int]=[1, 2],
                 bearing_classes:List[RadialBearing]=bearing_classes,
                 bearing_combination_simulations:List[BearingCombinationSimulation]=None,
                 catalog:BearingCatalog=None, name:str=''):
        
        if linkage_types == ['all']:
            linkage_types = ['ball_joint', 'cylindric_joint']
            
        self.radial_loads = radial_loads
        self.axial_loads = axial_loads
        self.speeds = speeds
        self.operating_times = operating_times
        self.inner_diameter = inner_diameter
        self.outer_diameter = outer_diameter
        self.length = length
        self.linkage_types = linkage_types
        self.mounting_types = mounting_types
        self.number_bearings = number_bearings
        self.bearing_classes = bearing_classes
        self.bearing_combination_simulations = bearing_combination_simulations
        self.catalog = catalog
        
        DessiaObject.__init__(self, name=name)
        
#    def __eq__(self, other_eb):
#        equal = (self.radial_loads == other_eb.radial_loads
#                 and self.axial_loads == other_eb.axial_loads
#                 and self.speeds == other_eb.speeds
#                 and self.operating_times == other_eb.operating_times
#                 and self.inner_diameter == other_eb.inner_diameter
#                 and self.outer_diameter == other_eb.outer_diameter
#                 and self.length == other_eb.length
#                 and self.linkage_types == other_eb.linkage_types
#                 and self.mounting_types == other_eb.mounting_types
#                 and self.number_bearings == other_eb.number_bearings
#                 and self.bearing_classes == other_eb.bearing_classes
#                 and self.catalog == other_eb.catalog)
#        return equal
#    
#    def __hash__(self):
#        h = int(sum(self.operating_times) % 230080000)
#        return h
        


#class ConceptualBearingCombinationOptimizer(protected_module.ConceptualBearingCombinationOptimizer if _open_source==True else object):

class ConceptualBearingCombinationOptimizer(protected_module.ConceptualBearingCombinationOptimizer if _open_source==True else DessiaObject):
        
    def __init__(self, linkage, mounting, d, D, length,
                 bearing_classes):
        
        self.bearing_classes = bearing_classes
        self.linkage = linkage
        self.mounting = mounting
        self.d = d
        self.D = D
        self.length = length
        

    
class BearingAssemblyOptimizer(protected_module.BearingAssemblyOptimizer if _open_source==True else DessiaObject):
    
    _dessia_methods = ['optimize']
    
    def __init__(self, loads:List[float], speeds:List[float], operating_times:List[float],
                 inner_diameters:List[float],
                 outer_diameters:List[float],
                 axial_positions:List[float],
                 lengths:List[float],
                 linkage_types=[['ball_joint', 'cylindric_joint'],
                                ['ball_joint', 'cylindric_joint']],
                 mounting_types=list(product(['both', 'right', 'left', 'free'],
                                             ['both', 'right', 'left', 'free'])),
                 number_bearings=[[1, 2], [1, 2]],
                 bearing_classes=bearing_classes,
                 bearing_assembly_simulations=None,
                 catalog=None):
                     
        self.linkage_types = linkage_types
        for i_linkage, linkage_type in enumerate(linkage_types):
            if linkage_type == ['all']:
                self.linkage_types[i_linkage] = ['ball_joint', 'cylindric_joint']
                     
        self.loads = loads
        self.speeds = speeds
        self.operating_times = operating_times
        self.inner_diameters = inner_diameters
        self.axial_positions = axial_positions
        self.outer_diameters = outer_diameters
        self.lengths = lengths
        self.mounting_types = mounting_types
        self.number_bearings = number_bearings
        self.bearing_classes = bearing_classes
        self.bearing_assembly_simulations = bearing_assembly_simulations
        self.catalog = catalog
        
    def __eq__(self, other_eb):
        equal = (self.loads == other_eb.loads
                 and self.speeds == other_eb.speeds
                 and self.operating_times == other_eb.operating_times
                 and self.inner_diameters == other_eb.inner_diameters
                 and self.outer_diameters == other_eb.outer_diameters
                 and self.axial_positions == other_eb.axial_positions
                 and self.lengths == other_eb.lengths
                 and self.linkage_types == other_eb.linkage_types
                 and self.mounting_types == other_eb.mounting_types
                 and self.number_bearings == other_eb.number_bearings
                 and self.bearing_classes == other_eb.bearing_classes
                 and self.catalog == other_eb.catalog)
        
        if (self.bearing_assembly_simulations is not None) and (other_eb.bearing_assembly_simulations is not None):
            for bearing_assembly_simulation, other_bearing_assembly_simulation in zip(self.bearing_assembly_simulations, other_eb.bearing_assembly_simulations):
                equal = equal and bearing_assembly_simulation == other_bearing_assembly_simulation
        elif (self.bearing_assembly_simulations is None) and (other_eb.bearing_assembly_simulations is None):
            pass
        elif (self.bearing_assembly_simulations is None) or (other_eb.bearing_assembly_simulations is None):
            equal = False
        return equal
    
    def __hash__(self):
        br_hash = int(sum(self.operating_times)/300000)
        br_hash += int(sum(self.outer_diameters*145))
#        for loads in self.loads:
#            for load in loads:
#                for item in load:
#                    br_hash += hash(tuple(item))
#        br_hash += hash(tuple(self.speeds)) + hash(tuple(self.operating_times))
#        br_hash += hash(tuple(self.inner_diameters)) + hash(tuple(self.outer_diameters))
#        br_hash += hash(tuple(self.axial_positions)) + hash(tuple(self.lengths))
#        for linkage_type in self.linkage_types:
#            br_hash += hash(tuple(linkage_type))
#        for mounting_type in self.mounting_types:
#            br_hash += hash(tuple(mounting_type))
#        for number_bearing in self.number_bearings:
#            br_hash += hash(tuple(number_bearing))
#        for bearing_classe in self.bearing_classes:
#            br_hash += hash(bearing_classe)
        br_hash += hash(self.catalog)
        return br_hash

        
    def Dict(self, subobjects_id = {}, stringify_keys=True):
        """
        Export dictionary
        """
        d = {}
        d['loads'] = self.loads
        d['speeds'] = self.speeds
        d['operating_times'] = self.operating_times
        d['inner_diameters'] = self.inner_diameters
        d['axial_positions'] = self.axial_positions
        d['outer_diameters'] = self.outer_diameters
        d['lengths'] = self.lengths
        d['linkage_types'] = self.linkage_types
        d['mounting_types'] = self.mounting_types
        d['number_bearings'] = self.number_bearings
                
                
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
        
        # TODO: serialize classes as string -> see another todo in dict to object
#        bc = []
#        for bearing_classe in self.bearing_classes:
#            bc.append(str(bearing_classe))
#        d['bearing_classes'] = bc
        
        d['catalog'] = self.catalog.Dict()
        
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
            
        if 'bearing_classes' in d:
            bearing_classes_ = []
            for bearing_classe in d['bearing_classes']:
                bearing_classes_.append(dict_bearing_classes[bearing_classe])
        else:
            bearing_classes_ = bearing_classes
            
        if not 'catalog' in d:
            catalog = schaeffler_catalog# TODO: change this??
        else:
            catalog = BearingCatalog.DictToObject(d['catalog'])
            
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
                 bearing_classes = bearing_classes_,
                 bearing_assembly_simulations = li_bar,
                 catalog = catalog)
        return obj


        