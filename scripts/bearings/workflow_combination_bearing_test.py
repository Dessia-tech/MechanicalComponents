#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2018

@author: Pierrem
"""
from dessia_common import DessiaObject, dict_merge
from dessia_common import Evolution, CombinationEvolution
import dectree
from dessia_api_client import Client
from dessia_common import workflow as wf
import dessia_common as dc
import plot_data

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

# filters = [
#           {'attribute' : 'bearing_combination.B', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'bearing_combination.d', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'bearing_combination.D', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'bearing_combination.number_bearing', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'bearing_combination.mass', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'bearing_combination.cost', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'bearing_combination_simulation_result.max_axial_load', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'bearing_combination_simulation_result.max_radial_load', 'operator' : 'gt', 'bound' : -100},
#            ]

# filter_analyze= wf.Filter(filters)

list_attribute = ['mass', 'cost', 'number_bearing', 'max_axial_load','max_radial_load']


display = wf.MultiPlot(list_attribute,order = 1, name= 'Display')

input_values = {}
blocks = []

blocks.extend([blockA, optimizeA, attribute_selection1, 

                display
                ])

pipes = [wf.Pipe(blockA.outputs[0], optimizeA.inputs[0]),
         wf.Pipe(optimizeA.outputs[1], attribute_selection1.inputs[0]),
          wf.Pipe(attribute_selection1.outputs[0], display.inputs[0])
         ]

workflow = wf.Workflow(blocks, pipes, attribute_selection1.outputs[0])


input_values = {workflow.index(blockA.inputs[0]): [500],
                workflow.index(blockA.inputs[1]): [0],
                workflow.index(blockA.inputs[2]): [1000],
                workflow.index(blockA.inputs[3]): [3600000],
                workflow.index(blockA.inputs[4]): 0.15,
                workflow.index(blockA.inputs[5]): 0.4,
                workflow.index(blockA.inputs[6]): 0.2,
                workflow.index(blockA.inputs[7]): [bearings.Linkage(ball_joint=True), bearings.Linkage(cylindric_joint=True)],
                workflow.index(blockA.inputs[8]): [bearings.Mounting(left=True), bearings.Mounting(right=True), bearings.Mounting(left=True, right=True), bearings.Mounting()],
                workflow.index(blockA.inputs[9]): [1, 2, 3],
                workflow.index(blockA.inputs[10]):  ['mechanical_components.bearings.RadialBallBearing', 
                                      'mechanical_components.bearings.AngularBallBearing',
                                      'mechanical_components.bearings.TaperedRollerBearing',
                                      'mechanical_components.bearings.NUP', 
                                      'mechanical_components.bearings.N', 
                                      'mechanical_components.bearings.NU'],
                workflow.index(blockA.inputs[12]): schaeffler_catalog,
                workflow.index(optimizeA.inputs[1]): 1000,
                }


##
workflow.plot_jointjs()
workflow_run = workflow.run(input_values)

a = workflow_run.to_dict()
obj = wf.WorkflowRun.dict_to_object(a)
    

