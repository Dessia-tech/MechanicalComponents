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

bis2 = bearings_opt.BearingAssemblyOptimizer(
                    loads = [[[[0.1595, 0, 0], [10, -14000, 0], [0, 0, 0]]]], 
                    speeds = [157.07],
                    operating_times = [3600000],
                    inner_diameters = [0.035, 0.035],
                    axial_positions = [0, 0.3], 
                    outer_diameters = [0.072, 0.072], 
                    lengths = [0.03, 0.03],
                    linkage_types = [bearings.SelectionLinkage([bearings.Linkage(ball_joint=True), bearings.Linkage(cylindric_joint=True)]),
                                     bearings.SelectionLinkage([bearings.Linkage(ball_joint=True), bearings.Linkage(cylindric_joint=True)])],
                    mounting_types = [bearings.CombinationMounting([bearings.Mounting(right=True), bearings.Mounting(left=True)])],
                    number_bearings = [[1], [1]],
                    catalog = schaeffler_catalog,
                    bearing_classes = [bearings.RadialBallBearing, 
                                      bearings.AngularBallBearing,
                                      bearings.TaperedRollerBearing,
                                      bearings.NUP, bearings.N, bearings.NU,
#                                       bearings_opt.NF
                                      ]
                    )

bis2.optimize(max_solutions = 10)

for num_sol, ba_simulation in enumerate(bis2.bearing_assembly_simulations):
    hash_ = hash(ba_simulation)
    equak = ba_simulation.bearing_assembly == ba_simulation.bearing_assembly
    d = ba_simulation.to_dict()
    obj = bearings.BearingAssemblySimulation.dict_to_object(d)
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
