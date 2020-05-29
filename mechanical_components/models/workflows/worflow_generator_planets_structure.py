#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 14:13:57 2020

@author: launay
"""
import mechanical_components.planetary_gears_generator as pg_generator
import dessia_common.workflow as wf
from dessia_common.workflow import Filter

block_planet_structure = wf.InstanciateModel(pg_generator.GeneratorPlanetsStructure, name='GeneratorPlanetsStructure')
generate_planet_structure = wf.ModelMethod(pg_generator.GeneratorPlanetsStructure,  'decision_tree', name='GeneratorPlanetsStructure-decision_tree')

block_planetary_gears_architecture = wf.InstanciateModel(pg_generator.GeneratorPlanetaryGearsArchitecture, name='GeneratorPlanetaryGearsArchitecture')
generate_planetary_gears_architecture = wf.ModelMethod(pg_generator.GeneratorPlanetaryGearsArchitecture,  'decision_tree', name='GeneratorPlanetaryGearsArchitecture-decision_tree')

block_planetary_gears_z_number = wf.InstanciateModel(pg_generator.GeneratorPlanetaryGearsZNumber, name='GeneratorPlanetaryGearsZNumber')
generate_planetary_gears_z_number = wf.ModelMethod(pg_generator.GeneratorPlanetaryGearsZNumber,  'decision_tree', name='GeneratorPlanetaryGearsZNumber')


# attribute_selection1 = wf.ModelAttribute('generator_Planets_structure')

# filters = [
#           {'attribute' : 'generator_planets_structure.number_max_planet', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'generator_planets_structure.number_junction', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'generator_planets_struture.number_max_junction_by_planet', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'generator_planets_struture.min_planet_branch', 'operator' : 'gt', 'bound' : -100},
#           {'attribute' : 'generator_planets_struture.number_max_meshing_plan', 'operator' : 'gt', 'bound' : -100},
#            ]

# filter_analyze= wf.Filter(filters)

input_values = {}
blocks = []

blocks.extend([block_planet_structure , generate_planet_structure, block_planetary_gears_architecture, 
                generate_planetary_gears_architecture])

pipes = [wf.Pipe(block_planet_structure.outputs[0], generate_planet_structure.inputs[0]),
          wf.Pipe(block_planetary_gears_architecture.outputs[0], generate_planetary_gears_architecture.inputs[0]),
          wf.Pipe(generate_planet_structure.outputs[0], block_planetary_gears_architecture.inputs[0]),
         ]
# wf.ForEach()
workflow_generator_planet_structure = wf.Workflow(blocks, pipes, generate_planetary_gears_architecture.outputs[0])

input_values = {workflow_generator_planet_structure.index(block_planet_structure.inputs[0]): 3,
                workflow_generator_planet_structure.index(block_planet_structure.inputs[1]): 0,
                workflow_generator_planet_structure.index(block_planet_structure.inputs[2]): 2,
                workflow_generator_planet_structure.index(block_planet_structure.inputs[3]): 1,
                workflow_generator_planet_structure.index(block_planet_structure.inputs[4]):2,
                workflow_generator_planet_structure.index(block_planetary_gears_architecture.inputs[1]):[[500,550],[600,650],[300,350],[200,250]] }
                # workflow_generator_planet_structure.index(block_planetary_gears_z_number.inputs[1]):[[500,550],[600,650],[300,350],[200,250]] , 
                # workflow_generator_planet_structure.index(block_planetary_gears_z_number.inputs[2]):[7, 80] ,
                # workflow_generator_planet_structure.index(block_planetary_gears_z_number.inputs[3]):[40,100] ,
                # workflow_generator_planet_structure.index(block_planetary_gears_z_number.inputs[4]):3 }

a = workflow_generator_planet_structure.to_dict()
obj = wf.Workflow.dict_to_object(a)
# ##
workflow_generator_planetary_gears_architecture = workflow_generator_planet_structure.run(input_values)
workflow_generator_planetary_gears_architecture.output_value[1].plot_kinematic_graph()