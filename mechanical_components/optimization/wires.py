#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# cython: language_level=3

"""

"""

#import networkx as nx
#import mechanical_components.wires as wires
#import matplotlib.pyplot as plt
import dessia_common as dc

try:
    _open_source = False
    import mechanical_components.optimization.wires_protected as protected_module

except (ModuleNotFoundError, ImportError) as e:
    _open_source = True

from mechanical_components.optimization.common import RoutingOptimizer


class WiringOptimizer(protected_module.WiringOptimizer if not _open_source
                      else RoutingOptimizer):
    def __init__(self, routes):
        super().__init__(routes)
