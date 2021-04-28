#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""
#import sys
import mechanical_components.bearings as bearings
from mechanical_components.optimization.bearings import BearingAssemblyOptimizer
from mechanical_components.models.catalogs import schaeffler_catalog
import plot_data

import pkg_resources
# from .catalogs import schaeffler_catalog

bearing_assembly_optimizer = BearingAssemblyOptimizer(
                    loads = [[[[0.1595, 0, 0], [0, -14000, 0], [0, 0, 0]]]], 
                    speeds = [157.07],
                    operating_times = [3600000],
                    inner_diameters = [0.035, 0.035],
                    axial_positions = [0, 0.3], 
                    outer_diameters = [0.072, 0.072], 
                    lengths = [0.03, 0.03],
                    linkage_types = [bearings.SelectionLinkage([bearings.Linkage(ball_joint=True), bearings.Linkage(cylindric_joint=True)]),
                                     bearings.SelectionLinkage([bearings.Linkage(ball_joint=True), bearings.Linkage(cylindric_joint=True)])],
                    mounting_types = [bearings.CombinationMounting([bearings.Mounting(), bearings.Mounting(left=True)])],
                    number_bearings = [[1], [1]],
                    catalog = schaeffler_catalog,
                    )

# bearing_assembly_optimizer.optimize(max_solutions = 10)
