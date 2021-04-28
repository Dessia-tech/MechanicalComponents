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
                    bearing_classes = ['mechanical_components.bearings.RadialBallBearing', 
                                      'mechanical_components.bearings.AngularBallBearing',
                                      'mechanical_components.bearings.TaperedRollerBearing',
                                      'mechanical_components.bearings.NUP', 
                                      'mechanical_components.bearings.N', 
                                      'mechanical_components.bearings.NU',
                                      # bearings_opt.NF
                                      ]
                    )

bis2.optimize(max_solutions = 10)

for num_sol, ba_simulation in enumerate(bis2.bearing_combination_simulations):
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

