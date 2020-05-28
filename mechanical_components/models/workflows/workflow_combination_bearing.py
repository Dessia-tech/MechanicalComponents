#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""
from dessia_common import DessiaObject, dict_merge
from dessia_common import Evolution, CombinationEvolution
import dectree
# from dessia_api_client import Client
from dessia_common import workflow as wf
import dessia_common as dc
from volmdlr import plot_data

import mechanical_components.bearings as bearings
import mechanical_components.optimization.bearings as bearings_opt
import mechanical_components
schaeffler_catalog = mechanical_components.models.schaeffler_catalog

from itertools import product
import networkx as nx
from random import random
from scipy.optimize import minimize, fsolve
from dataclasses import dataclass
from numpy import allclose
import matplotlib.pyplot as plt
from scipy import interpolate
import math
import json
import pkg_resources

blockA = wf.InstanciateModel(bearings_opt.BearingCombinationOptimizer,  name='BearingCombinationOptimizer')
optimizeA = wf.ModelMethod(bearings_opt.BearingCombinationOptimizer, 'optimize', name='BearingCombinationOptimizer-optimize')
attribute_selection1 = wf.ModelAttribute('bearing_combination_simulations')

filters = [
          {'attribute' : 'bearing_combination.B', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'bearing_combination.d', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'bearing_combination.D', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'bearing_combination.number_bearing', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'bearing_combination.mass', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'bearing_combination.cost', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'bearing_combination_simulation_result.max_axial_load', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'bearing_combination_simulation_result.max_radial_load', 'operator' : 'gt', 'bound' : -100},
           ]

filter_analyze= wf.Filter(filters)

input_values = {}
blocks = []

blocks.extend([blockA, optimizeA, attribute_selection1, 
                filter_analyze
                ])

pipes = [wf.Pipe(blockA.outputs[0], optimizeA.inputs[0]),
         wf.Pipe(optimizeA.outputs[1], attribute_selection1.inputs[0]),
          wf.Pipe(attribute_selection1.outputs[0], filter_analyze.inputs[0]),
         ]

workflow_combination_bearing = wf.Workflow(blocks, pipes, filter_analyze.outputs[0])

input_values = {workflow.index(blockA.inputs[0]): [500],
                workflow.index(blockA.inputs[1]): [0],
                workflow.index(blockA.inputs[2]): [1000],
                workflow.index(blockA.inputs[3]): [1e6],
                workflow.index(blockA.inputs[4]): 0.15,
                workflow.index(blockA.inputs[5]): 0.4,
                workflow.index(blockA.inputs[6]): 0.2,
                workflow.index(blockA.inputs[7]): [bearings.Linkage(ball_joint=True), bearings.Linkage(cylindric_joint=True)],
                workflow.index(blockA.inputs[8]): [bearings.Mounting(left=True), bearings.Mounting(right=True), bearings.Mounting(left=True, right=True), bearings.Mounting()],
                workflow.index(blockA.inputs[9]): [1, 2, 3],
                workflow.index(blockA.inputs[10]): [bearings.RadialBallBearing,
                                                    bearings.AngularBallBearing,
                                                    bearings.TaperedRollerBearing,
                                                    bearings.NUP],
                workflow.index(blockA.inputs[12]): schaeffler_catalog,
                workflow.index(optimizeA.inputs[1]): 10,
                }


workflow_combination_bearing = workflow.run(input_values)
