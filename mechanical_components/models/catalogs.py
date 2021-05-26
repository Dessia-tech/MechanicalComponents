#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import mechanical_components.bearings as bearings


import pkg_resources

with pkg_resources.resource_stream(pkg_resources.Requirement('mechanical_components'),
                           'mechanical_components/catalogs/schaeffler.json') as schaeffler_json:
    schaeffler_catalog = bearings.BearingCatalog.load_from_file(schaeffler_json)

with pkg_resources.resource_stream(pkg_resources.Requirement('mechanical_components'),
                            'mechanical_components/catalogs/ntn.json') as ntn_json:
    ntn_catalog = bearings.BearingCatalog.load_from_file(ntn_json)
