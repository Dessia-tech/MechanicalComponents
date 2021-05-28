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
                    radial_loads = [6817.5526036989395, 4746.947246540115, 3157.1905110758025, 1263.3588630217332, 758.0756501369663, 630.4727850323359, 352.5982680434697, 424.03310100685485, 706.7218350114246, 1766.0419380663182, 3813.2473112127177], 
                    axial_loads = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    speeds = [467.75533333333334, 254.23533333333336, 335.5613333333333, 839.008, 1271.1766666666667, 1677.9113333333335, 1677.9113333333335, 1271.1766666666667, 839.008, 335.5613333333333, 467.75533333333334],
                    operating_times = [1440000, 1242503.9999999998, 1882692.0, 188280.0, 198792.0, 376524.00000000006, 170208.0, 205956.0, 113436.00000000001, 425448.0, 508751.99999999994],
                    inner_diameter = 0.01,
                    outer_diameter = 0.150, 
                    length = 0.08,
                    linkage_types = [bearings.Linkage(ball_joint=True), bearings.Linkage(cylindric_joint=True)],
                    mounting_types = [ bearings.Mounting()],
                    number_bearings =[1, 2],
                    bearing_classes = ['mechanical_components.bearings.RadialBallBearing', 
                                      # 'mechanical_components.bearings.AngularBallBearing',
                                      # 'mechanical_components.bearings.TaperedRollerBearing',
                                      # 'mechanical_components.bearings.NUP', 
                                      # 'mechanical_components.bearings.N', 
                                      # 'mechanical_components.bearings.NU',
                                      # bearings_opt.NF
                                      ]
                    )

bis2.optimize(max_solutions = 10,max_speed=1600)

for num_sol, ba_simulation in enumerate(bis2.bearing_combination_simulations):
    print(ba_simulation.bearing_combination.bearings)
    hash_ = hash(ba_simulation)
    equak = ba_simulation.bearing_combination == ba_simulation.bearing_combination
    d = ba_simulation.to_dict()
    obj = bearings.BearingCombinationSimulation.dict_to_object(d)
    ba_simulation == obj
    
    plot_data.plot_canvas(ba_simulation.plot_data()[0])


# d = bis2.to_dict()

# obj = bearings_opt.BearingCombinationOptimizer.dict_to_object(d)
# print(obj.bearing_classes)

# print(bis2.bearing_classes)

# d = bis2.to_dict()
# obj = dc.dict_to_object(d)

# if not obj == bis2:
#     raise KeyError('Non esqual object BearingCombinationOptimizer with dict_to_object')

