#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""
#import sys
import mechanical_components.bearings as bearings
import mechanical_components.optimization.bearings as bearings_opt
# from dessia_api_client import Client

import plot_data

import pkg_resources

with pkg_resources.resource_stream(pkg_resources.Requirement('mechanical_components'),
                           'mechanical_components/catalogs/schaeffler.json') as schaeffler_json:
    schaeffler_catalog = bearings.BearingCatalog.load_from_file(schaeffler_json)
from mechanical_components.models.catalogs import schaeffler_catalog,ntn_catalog

loads= [[[[0.09385651325269537, 0, 0], 
          [4441.17255936, 8704.26197882,  102.90050552], 
          [0, 0, 0]], 
         [[0.052499999999810316, 0, 0], 
          [-12873.35771282, -18195.29729408, -17481.39271479],
          [0, 0, 0]]], 
        [[[0.09385651325269537, 0, 0], 
          [3092.3137785 , 6060.631261  ,   71.64789181], 
          [0, 0, 0]], 
         [[0.052499999999810316, 0, 0],
          [ -8963.50252075, -12669.07956724, -12171.99981239], 
          [0, 0, 0]]],
        [[[0.09385651325269537, 0, 0],
          [2056.69522152, 4030.92061373,   47.65295092], 
          [0, 0, 0]], 
         [[0.052499999999810316, 0, 0],
          [-5961.61777975, -8426.19387079, -8095.58655544], 
          [0, 0, 0]]], 
        [[[0.09385651325269537, 0, 0], 
          [ 822.9925079 , 1612.98447643,   19.06846536], 
          [0, 0, 0]], 
         [[0.052499999999810316, 0, 0], 
          [-2385.55850006, -3371.76571096, -3239.4722429 ], 
          [0, 0, 0]]],
        [[[0.09385651325269537, 0, 0], 
          [493.83480715, 967.86771472,  11.44198984], 
          [0, 0, 0]], 
         [[0.052499999999810316, 0, 0], 
          [-1431.44902355, -2023.2204469 , -1943.83804833], 
          [0, 0, 0]]], 
        [[[0.09385651325269537, 0, 0], 
          [410.71020571, 804.95166087,   9.5160202 ],
          [0, 0, 0]],
         [[0.052499999999810316, 0, 0],
          [-1190.50077964, -1682.66244888, -1616.64206964], 
          [0, 0, 0]]], 
        [[[0.09385651325269537, 0, 0], 
          [-272.56222743, -293.94093287,  446.09765052], 
          [0, 0, 0]],
         [[0.052499999999810316, 0, 0], 
          [ 790.0596083 , -299.44002282, 1519.3205074 ], 
          [0, 0, 0]]],
        [[[0.09385651325269537, 0, 0],
          [-327.78211633, -353.49205193,  536.47504042], 
          [0, 0, 0]], 
         [[0.052499999999810316, 0, 0], 
          [ 950.12215332, -360.10523293, 1827.12805072], 
          [0, 0, 0]]], 
        [[[0.09385651325269537, 0, 0], 
          [-546.30352721, -589.15341989,  894.12506737], 
          [0, 0, 0]], 
         [[0.052499999999810316, 0, 0], 
          [1583.5369222 , -600.17538821, 3045.21341786], 
          [0, 0, 0]]], 
        [[[0.09385651325269537, 0, 0], 
          [-1365.16928185, -1472.24777266,  2234.34778527], 
          [0, 0, 0]],
         [[0.052499999999810316, 0, 0], 
          [ 3957.1334527 , -1499.79079925,  7609.74734313], 
          [0, 0, 0]]], 
        [[[0.09385651325269537, 0, 0], 
          [-2947.68090223, -3178.88535913,  4824.41583115], 
          [0, 0, 0]], 
         [[0.052499999999810316, 0, 0], 
         [ 8544.26396869, -3238.35641122, 16431.00765034], 
         [0, 0, 0]]]]
operating_time=[400,345.14,522.97,52.3,55.22,04.59,47.28,57.21,31.51,118.18,141.32]

bis2 = bearings_opt.BearingAssemblyOptimizer(
                    loads = loads, 
                    speeds = [142.7452128307315, 77.58516938148286, 102.40348004818196, 256.0406413182242, 387.9258469074142, 512.0493414386791, 512.0493414386791, 387.9258469074142, 256.0406413182242, 102.40348004818196, 142.7452128307315],
                    operating_times = operating_time,
                    inner_diameters = [0.005, 0.005],
                    axial_positions = [-0.01, 0.13], 
                    outer_diameters = [0.2, 0.2], 
                    lengths = [0.08, 0.08],
                    linkage_types = [bearings.SelectionLinkage([bearings.Linkage(ball_joint=True), bearings.Linkage(cylindric_joint=True)]),
                                      bearings.SelectionLinkage([bearings.Linkage(ball_joint=True), bearings.Linkage(cylindric_joint=True)])],
                    mounting_types = [bearings.CombinationMounting([bearings.Mounting(right=True), bearings.Mounting(left=True)])],
                    number_bearings = [[1], [1]],
                    catalog = ntn_catalog,
                    bearing_classes = [
                        bearings.RadialBallBearing, 
                                      bearings.AngularBallBearing,
                                       bearings.TaperedRollerBearing,
                                      # bearings.NUP, bearings.N, bearings.NU,
#                                       bearings_opt.NF
                                      ]
                    )

bis2.optimize(max_solutions = 10,verbose=True)

for num_sol, ba_simulation in enumerate(bis2.bearing_assembly_simulations):
    hash_ = hash(ba_simulation)
    equak = ba_simulation.bearing_assembly == ba_simulation.bearing_assembly
    d = ba_simulation.to_dict()
    obj = bearings.BearingAssemblySimulation.dict_to_object(d)
    print(obj.bearing_assembly.D)
    ba_simulation == obj
    
# ba_simulation.bearing_assembly.bearing_combinations[0].plot()
# plots = ba_simulation.bearing_assembly.plot_data()
# pdg = plot_data.plot_canvas(plots[0])

# d = bis2.to_dict()
# obj = bearings_opt.BearingAssemblyOptimizer.dict_to_object(d)

# if not obj == bis2:
#     raise KeyError('Non esqual object BearingAssemblyOptimizer with dict_to_object')
    
# vol1 = ba_simulation.bearing_assembly.bearing_combinations[0].bearings[0].volmdlr_volume_model()
# vol1.babylonjs()    


# vol1 = ba_simulation.bearing_assembly.bearing_combinations[0].volmdlr_volume_model()
# vol1.babylonjs()

# vol1 = ba_simulation.bearing_assembly.bearing_combinations[0].volmdlr_volume_model()
#vol1.babylonjs()   

# vol1 = ba_simulation.volmdlr_volume_model().babylonjs()
# plot_data.plot_canvas(ba_simulation.plot_data()[0])

#c = Client()
#c.api_url = 'http://localhost:5000'
## c.api_url = 'https://api.platform.dessia.tech'
#r = c.CreateObject(bis2)
