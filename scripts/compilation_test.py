#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compilation test
"""

import sys


# =============================================================================
#  Tests
# =============================================================================

if sys.argv[1] == 'test_valid_license':
    license_should_be_valid = True
    print('Testing for valid license')
elif sys.argv[1] == 'test_unvalid_license':
    license_should_be_valid = False
    print('Testing for unvalid license')
else:
    raise ValueError('Invalid option: use either test_valid_license ot test_unvalid_license')


import mechanical_components.bearings as bearings
import mechanical_components.optimization.bearings as bearings_opt
# from dessia_api_client import Client


import pkg_resources

with pkg_resources.resource_stream(pkg_resources.Requirement('mechanical_components'),
                           'mechanical_components/catalogs/schaeffler.json') as schaeffler_json:
    schaeffler_catalog = bearings.BearingCatalog.load_from_file(schaeffler_json)

bis2 = bearings_opt.BearingAssemblyOptimizer(
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
                    catalog = schaeffler_catalog)






test_passed = license_should_be_valid
try:
    bis2.optimize(max_solutions = 10)
    print('No license error triggered')
except RuntimeError:
    print('License error triggered')
    if license_should_be_valid:
        print('License was expected to be valid, test failed')
        test_passed = False
    else:
        print('License was expected to be unvalid, test passed')
        test_passed = True
    
assert test_passed