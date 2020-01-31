#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""
#import sys
import mechanical_components.bearings as bearings
import mechanical_components.optimization.bearings as bearings_opt

import pkg_resources

#import numpy as npy

#bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer(
#                    loads = [[[(0.1, 0, 0), (2000, 10000, 0), (0, 0, 0)], 
#                              [(0.3, 0, 0), (1000, 3000, 0), (0, 0, 0)]]], 
#                    speeds = [1000],
#                    operating_times = [10000*3600],
#                    inner_diameters = [0.02, 0.02],
#                    axial_positions = [0, 0.2], 
#                    outer_diameters = [0.1, 0.1], 
#                    lengths = [0.08, 0.08],
#                    linkage_types = [['cylindric_joint'], ['cylindric_joint']],
#                    mounting_types = [['free', 'both'], ['right', 'left']],
#                    number_bearings = [[2], [2]],
#                    bearing_classes = [bearings.RadialBallBearing, 
#                                       bearings.AngularBallBearing,
#                                       bearings.TaperedRollerBearing,
#                                       bearings.NUP, bearings.N, bearings.NU,
##                                       bearings_opt.NF
#                                       ])


bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer(
                    loads = [[[[0.15, 0, 0], [0, 2055, 0], [0, 0, 0]]]], 
                    speeds = [500],
                    operating_times = [100000000],
                    inner_diameters = [0.02, 0.025],
                    axial_positions = [0, 0.3], 
                    outer_diameters = [0.1, 0.1], 
                    lengths = [0.1, 0.1],
                    linkage_types = [['cylindric_joint'], ['ball_joint']],
                    mounting_types = [['left', 'right']],
                    number_bearings = [[1, 2], [1]],

#                    bearing_classes = [bearings.RadialBallBearing, 
#                                       bearings.AngularBallBearing,
#                                       bearings.TaperedRollerBearing,
#                                       bearings.NUP, bearings.N, bearings.NU,
##                                       bearings_opt.NF
#                                       ]
                    )

bis = bearings_opt.BearingAssemblyOptimizer(
                    loads = [[[[0.15, 0, 0], [0, 2000, 0], [0, 0, 0]]]], 
                    speeds = [500],
                    operating_times = [10000000],
                    inner_diameters = [0.03, 0.03],
                    axial_positions = [0, 0.3], 
                    outer_diameters = [0.065, 0.065], 
                    lengths = [0.015, 0.015],
                    linkage_types = [['all'], ['all']],
                    mounting_types = [['left', 'right']],
                    number_bearings = [[1, 2], [1, 2]],

#                    bearing_classes = [bearings.RadialBallBearing, 
#                                       bearings.AngularBallBearing,
#                                       bearings.TaperedRollerBearing,
#                                       bearings.NUP, bearings.N, bearings.NU,
##                                       bearings_opt.NF
#                                       ]
                    )

with pkg_resources.resource_stream(pkg_resources.Requirement('mechanical_components'),
                           'mechanical_components/catalogs/schaeffler.json') as schaeffler_json:
    schaeffler_catalog = bearings.BearingCatalog.LoadFromFile(schaeffler_json)

bis2 = bearings_opt.BearingAssemblyOptimizer(
                    loads = [[[[0.1595, 0, 0], [0, -14000, 0], [0, 0, 0]]]], 
                    speeds = [157.07],
                    operating_times = [3600000],
                    inner_diameters = [0.035, 0.035],
                    axial_positions = [0, 0.3], 
                    outer_diameters = [0.072, 0.072], 
                    lengths = [0.03, 0.03],
                    linkage_types = [['ball_joint', 'cylindric_joint'], ['ball_joint', 'cylindric_joint']],
                    mounting_types = [['free', 'left']],
                    number_bearings = [[1], [1]],
                    catalog = schaeffler_catalog,
#                    bearing_classes = [bearings.RadialBallBearing, 
#                                       bearings.AngularBallBearing,
#                                       bearings.TaperedRollerBearing,
#                                       bearings.NUP, bearings.N, bearings.NU,
##                                       bearings_opt.NF
#                                       ]
                    )

#'axial_positions': ,
# 'bearing_assembly_simulations': [],
# 'inner_diameters': ,
# 'lengths': ,
# 'linkage_types': ,
# 'loads': ,
# 'mounting_types': ,
# 'number_bearings': ,
# 'operating_times': ,
# 'outer_diameters': ,
# 'speeds': 

#print(hash(bearing_assembly_opt))
#print(bearing_assembly_opt == bis)

#d = bearing_assembly_opt.Dict()
#del bearing_assembly_opt
#bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer.DictToObject(d)

#bag = bis2.OptimizeGeneric(10)

#for num_sol, bas in enumerate(bag):
#    for ba in bas:
#    #    print(num_sol, ba_simulation.bearing_assembly.mass, ba_simulation.bearing_assembly_simulation_result.L10)
##        ba.Plot()    
#        print(num_sol, ba.mass, ba.cost)

bis2.Optimize(max_solutions = 10)
#0 0.396 4.6776 97.80460440116983
for num_sol, ba_simulation in enumerate(bis2.bearing_assembly_simulations):
#    print(num_sol, ba_simulation.bearing_assembly.mass, ba_simulation.bearing_assembly_simulation_result.L10)
    ba_simulation.bearing_assembly.Plot()    
    print(num_sol, ba_simulation.bearing_assembly.mass, ba_simulation.bearing_assembly.cost, ba_simulation.bearing_assembly_simulation_result.L10)
#    print(hash(ba_simulation))
    equal = (ba_simulation == ba_simulation)
#    
#print(bearing_assembly_opt == bis2)
    
#for num_sol, ba_simulation in enumerate(r.bearing_assembly_simulations):
#    print(num_sol, ba_simulation.bearing_assembly.mass, ba_simulation.bearing_assembly_simulation_result.L10)
#    ba_simulation.bearing_assembly.Plot()    
#    print(hash(ba_simulation))
#    print(ba_simulation == ba_simulation)
#    
#from dessia_api_client import Client
#c=Client()
#c.api_url='https://api-dev.software.dessia.tech'
#r=c.CreateObject(bearing_assembly_opt)
#r=c.GetAllClassObjects('mechanical_components.optimization.bearings.BearingAssemblyOptimizer')
#print(r)

    
#d = bearing_assembly_opt.Dict()
#del bearing_assembly_opt
#bearing_assembly_opt = bearings_opt.BearingAssemblyOptimizer.DictToObject(d)
#bearing_assembly_opt.Optimize(3)