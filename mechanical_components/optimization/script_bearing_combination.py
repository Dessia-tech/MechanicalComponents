#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 18:12:27 2021

@author: dasilva
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""
#import sys
import dessia_common as dc
import mechanical_components.bearings as bearings
import mechanical_components.optimization.bearings as bearings_opt
# from dessia_api_client import Client

import plot_data

import pkg_resources

with pkg_resources.resource_stream(pkg_resources.Requirement('mechanical_components'),
                           'mechanical_components/catalogs/schaeffler.json') as schaeffler_json:
    schaeffler_catalog = bearings.BearingCatalog.load_from_file(schaeffler_json)
    

bis2 = bearings_opt.BearingCombinationOptimizer(
                    radial_loads = [500], 
                    axial_loads = [0],
                    speeds = [1000],
                    operating_times = [3600000],
                    inner_diameter = 0.15,
                    outer_diameter = 0.4, 
                    length = 0.2,
                    linkage_types = [bearings.Linkage(ball_joint=True), bearings.Linkage(cylindric_joint=True)],
                    mounting_types = [bearings.Mounting(left=True), bearings.Mounting(right=True), bearings.Mounting(left=True, right=True), bearings.Mounting()],
                    number_bearings =[1, 2, 3],
                    bearing_classes = [bearings.RadialBallBearing, 
                                      bearings.AngularBallBearing,
                                      bearings.TaperedRollerBearing,
                                      bearings.NUP, bearings.N, bearings.NU,
                                      # bearings_opt.NF
                                      ]
                    )

bis2.optimize(max_solutions = 10)

# for num_sol, ba_simulation in enumerate(bis2.bearing_combination_simulations):
#     hash_ = hash(ba_simulation)
#     equak = ba_simulation.bearing_combination == ba_simulation.bearing_combination
#     d = ba_simulation.to_dict()
#     obj = bearings.BearingCombinationSimulation.dict_to_object(d)
#     ba_simulation == obj
    
# plot_data.plot_canvas(ba_simulation.plot_data()[0][0])

d = bis2.to_dict()
obj = bearings_opt.BearingCombinationOptimizer.dict_to_object(d)


# d = bis2.to_dict()
# obj = dc.dict_to_object(d)

if not obj == bis2:
    raise KeyError('Non esqual object BearingAssemblyOptimizer with dict_to_object')


    # def to_dict(self, subobjects_id = {}, stringify_keys=True):
    #     """
    #     Export dictionary
    #     """
    #     d = {}
    #     d['radial_loads'] = self.radial_loads
    #     d['axial_loads'] = self.axial_loads
    #     d['speeds'] = self.speeds
    #     d['operating_times'] = self.operating_times
    #     d['inner_diameter'] = self.inner_diameter
    #     d['outer_diameter'] = self.outer_diameter
    #     d['length'] = self.length
    #     d['linkage_types'] = [lt.to_dict() for lt in self.linkage_types]
    #     d['mounting_types'] = [mt.to_dict() for mt in self.mounting_types]
    #     d['number_bearings'] = self.number_bearings
    #     d['bearing_classes'] = [bc.__module__ + '.' + bc.__name__ for bc in self.bearing_classes]
        
    #     if self.bearing_combination_simulations is not None:
    #         bar_dict = []
    #         for bar in self.bearing_combination_simulations:
    #             if bar in subobjects_id:
    #                 bar_dict.append(subobjects_id[bar])
    #             else:                
    #                 bar_dict.append(bar.to_dict())
    #     else:
    #         bar_dict = None
    #     d['bearing_combination_simulations'] = bar_dict
    #     d['catalog'] = self.catalog.to_dict()
    #     d['name'] = self.name
    #     d['object_class'] = 'mechanical_components.optimization.bearings.BearingCombinationOptimizer'
        
    #     if stringify_keys:
    #         return StringifyDictKeys(d)

    #     return d
    
    # @classmethod
    # def dict_to_object(cls, d):
        
    #     if 'bearing_combination_simulations' in d:
    #         if d['bearing_combination_simulations'] is None:
    #             li_bar = None
    #         else:
    #             li_bar = []
    #             for bar in d['bearing_combination_simulations']:
    #                 li_bar.append(BearingCombinationSimulation.dict_to_object(bar))
    #     else:
    #         li_bar = None
            
    #     # if 'bearing_classes' in d:
    #     #     bearing_classes_ = []
    #     #     for bearing_classe in d['bearing_classes']:
    #     #         bearing_classes_.append(dict_bearing_classes[bearing_classe])
    #     # else:
    #     #     bearing_classes_ = bearing_classes
            
    #     if not 'catalog' in d:
    #         catalog = schaeffler_catalog# TODO: change this??
    #     else:
    #         catalog = BearingCatalog.dict_to_object(d['catalog'])
            
    #     bearing_classes = []
    #     for bearing_classe in d['bearing_classes']:
    #         module = bearing_classe.split('.')
    #         mod = ''
    #         for m in module[0:-1]:
    #             mod += m + '.'
    #         bearing_classes.append(getattr(import_module(mod[0:-1]), module[-1]))
        
    #     obj = cls(radial_loads = d['radial_loads'], 
    #               axial_loads = d['axial_loads'],
    #               speeds = d['speeds'], 
    #               operating_times = d['operating_times'],
    #               inner_diameter = d['inner_diameter'],
    #               outer_diameter = d['outer_diameter'],
    #               length = d['length'],
    #               linkage_types = [Linkage.dict_to_object(lt) for lt in d['linkage_types']],
    #               mounting_types = [Mounting.dict_to_object(mt) for mt in d['mounting_types']],
    #               number_bearings = d['number_bearings'],
    #               bearing_classes = bearing_classes,
    #               bearing_combination_simulations = li_bar,
    #               catalog = catalog, name = d['name'])
    #     return obj
            
        