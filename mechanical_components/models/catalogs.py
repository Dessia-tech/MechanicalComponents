#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 14:16:58 2020

@author: dumouchel
"""

import mechanical_components.bearings as bearings
import mechanical_components.optimization.bearings as bearings_opt

from volmdlr import plot_data

import pkg_resources

with pkg_resources.resource_stream(pkg_resources.Requirement('mechanical_components'),
                           'mechanical_components/catalogs/schaeffler.json') as schaeffler_json:
    schaeffler_catalog = bearings.BearingCatalog.load_from_file(schaeffler_json)