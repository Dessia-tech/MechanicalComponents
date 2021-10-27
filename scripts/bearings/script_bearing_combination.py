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
                    radial_loads = [7232.741086672101, 5036.036006189049, 3349.4631952703703, 1340.2973305744622, 804.2424049028032, 668.8685341214476, 739.5367613853548, 889.3636034540525, 1482.2726724234185, 3704.0821062321766, 7997.874131781144], 
                    axial_loads = [6436.67886590166, 4481.751266987825, 2980.8088942705544, 1192.779251784766, 715.7245128323597, 595.2503906949066, 395.02980473508575, 475.0610773598313, 791.7684612317834, 1978.5667292635906, 4272.131990643152],
                    speeds = [34.49695581168109, 18.74985582156486, 24.747648317800305, 61.87683996116258, 93.74927910782431, 123.74596075566335, 123.74596075566335, 93.74927910782431, 61.87683996116258, 24.747648317800305, 34.49695581168109],
                    operating_times = [1440000, 1242503.9999999998, 1882692.0, 188280.0, 198792.0, 376524.00000000006, 170208.0, 205956.0, 113436.00000000001, 425448.0, 508751.99999999994],
                    inner_diameter = 0.005,
                    outer_diameter = 0.150, 
                    length = 0.08,
                    linkage_types = [bearings.Linkage(ball_joint=True), bearings.Linkage(cylindric_joint=True)],
                    mounting_types = [bearings.Mounting(right=True)],
                    number_bearings =[1, 2],
                    bearing_classes = [
                                        'mechanical_components.bearings.RadialBallBearing', 


                                        # 'mechanical_components.bearings.AngularBallBearing',
                                         # 'mechanical_components.bearings.TaperedRollerBearing',
                                         # 'mechanical_components.bearings.NUP', 
                                       #  'mechanical_components.bearings.N', 
                                       #  'mechanical_components.bearings.NU',
                                       # 'mechanical_components.bearings.AngularBallBearing',
                                      # 'mechanical_components.bearings.TaperedRollerBearing',
                                      # 'mechanical_components.bearings.NUP', 
                                      # 'mechanical_components.bearings.N', 
                                      # 'mechanical_components.bearings.NU',

                                      # bearings_opt.NF
                                      ]
                    )


[6436.67886590166, 4481.751266987825, 2980.8088942705544, 1192.779251784766, 715.7245128323597, 595.2503906949066, 395.02980473508575, 475.0610773598313, 791.7684612317834, 1978.5667292635906, 4272.131990643152]
[7232.741086672101, 5036.036006189049, 3349.4631952703703, 1340.2973305744622, 804.2424049028032, 668.8685341214476, 739.5367613853548, 889.3636034540525, 1482.2726724234185, 3704.0821062321766, 7997.874131781144]
[34.49695581168109, 18.74985582156486, 24.747648317800305, 61.87683996116258, 93.74927910782431, 123.74596075566335, 123.74596075566335, 93.74927910782431, 61.87683996116258, 24.747648317800305, 34.49695581168109]
[1440000, 1242503.9999999998, 1882692.0, 188280.0, 198792.0, 376524.00000000006, 170208.0, 205956.0, 113436.00000000001, 425448.0, 508751.99999999994]
bis2.optimize(max_solutions = 10)

for num_sol, ba_simulation in enumerate(bis2.bearing_combination_simulations):
    print(ba_simulation.bearing_combination.bearings)
    hash_ = hash(ba_simulation)
    equak = ba_simulation.bearing_combination == ba_simulation.bearing_combination
    d = ba_simulation.to_dict()
    obj = bearings.BearingCombinationSimulation.dict_to_object(d)
    ba_simulation == obj
    ba_simulation.plot()
# print(bearings.RadialBallBearing.estimate_base_life_time([9431.539915048863, 6567.022659434021, 4367.721095172622, 1747.7561577090514, 1048.7371595804343, 872.2087797545174, 638.8705251685427, 768.3028377657752, 1280.5047296096266, 3199.8799843877973, 6909.198181346897], 
#                                                    [467.75533333333334, 254.23533333333336, 335.5613333333333, 839.008, 1271.1766666666667, 1677.9113333333335, 1677.9113333333335, 1271.1766666666667, 839.008, 335.5613333333333, 467.75533333333334], 
#                                                    [1440000, 1242503.9999999998, 1882692.0, 188280.0, 198792.0, 376524.00000000006, 170208.0, 205956.0, 113436.00000000001, 425448.0, 508751.99999999994], 
#                                                    20000.0))

    
    

# d = bis2.to_dict()

# obj = bearings_opt.BearingCombinationOptimizer.dict_to_object(d)
# print(obj.bearing_classes)

# print(bis2.bearing_classes)

# d = bis2.to_dict()
# obj = dc.dict_to_object(d)

# if not obj == bis2:
#     raise KeyError('Non esqual object BearingCombinationOptimizer with dict_to_object')

