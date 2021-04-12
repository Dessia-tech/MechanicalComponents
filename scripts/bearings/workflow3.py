#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 10:48:16 2021

@author: dasilva
"""
from dessia_api_client import Client
from dessia_common import workflow as wf

import plot_data
import mechanical_components.bearings as bearings
import mechanical_components.optimization.bearings as bearings_opt
import mechanical_components
schaeffler_catalog = mechanical_components.models.schaeffler_catalog

block_optimizer = wf.InstanciateModel(bearings_opt.BearingAssemblyOptimizer,  name='BearingAssemblyOptimizer')
block_optimize = wf.ModelMethod(bearings_opt.BearingAssemblyOptimizer, 'optimize', name='BearingAssemblyOptimizer-optimize')
attribute_selection1 = wf.ModelAttribute('bearing_assembly_simulations') 

workflow_block = [block_optimizer, block_optimize, attribute_selection1] 

workflow_pipe = [wf.Pipe(block_optimizer.outputs[0], block_optimize.inputs[0]),
                 wf.Pipe(block_optimize.outputs[1], attribute_selection1.inputs[0])]

workflow = wf.Workflow(workflow_block, workflow_pipe,attribute_selection1.outputs[0])
workflow.plot_jointjs()

input_values = {workflow.index(block_optimizer.inputs[0]): [[[[0.1595, 0, 0], [0, -14000, 0], [0, 0, 0]]]],
                workflow.index(block_optimizer.inputs[1]): [157.07],
                workflow.index(block_optimizer.inputs[2]): [3600000],
                workflow.index(block_optimizer.inputs[3]): [0.035, 0.035],
                workflow.index(block_optimizer.inputs[4]): [0.072, 0.072],
                workflow.index(block_optimizer.inputs[5]): [0, 0.3],
                workflow.index(block_optimizer.inputs[6]): [0.1, 0.1],
                workflow.index(block_optimizer.inputs[7]): [bearings.SelectionLinkage([bearings.Linkage(ball_joint=True), bearings.Linkage(cylindric_joint=True)]),
                                                   bearings.SelectionLinkage([bearings.Linkage(ball_joint=True), bearings.Linkage(cylindric_joint=True)])],
                workflow.index(block_optimizer.inputs[8]): [bearings.CombinationMounting([bearings.Mounting(), bearings.Mounting(left=True)])],
                workflow.index(block_optimizer.inputs[9]): [[1, 2], [1, 2]],
                workflow.index(block_optimizer.inputs[12]): schaeffler_catalog,
                workflow.index(block_optimize.inputs[1]): 10,
                }

##
workflow_run = workflow.run(input_values)

d1 = workflow_run.to_dict()
obj = wf.WorkflowRun.dict_to_object(d1)


import json

object1=json.dumps(d1)

object2=json.loads(object1)

c = Client(api_url = 'https://api.demo.dessia.tech')
r = c.create_object_from_python_object(workflow_run)
