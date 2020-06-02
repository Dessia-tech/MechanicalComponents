#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 14:13:57 2020

@author: launay
"""
import mechanical_components.planetary_gears_generator as pg_generator
import dessia_common.workflow as wf
from dessia_common.workflow import Filter
import mechanical_components.planetary_gears as pg
block_planet_structure = wf.InstanciateModel(pg_generator.GeneratorPlanetsStructure, name='GeneratorPlanetsStructure')
generate_planet_structure = wf.ModelMethod(pg_generator.GeneratorPlanetsStructure,  'decision_tree', name='GeneratorPlanetsStructure-decision_tree')

block_planetary_gears_architecture = wf.InstanciateModel(pg_generator.GeneratorPlanetaryGearsArchitecture, name='GeneratorPlanetaryGearsArchitecture')
generate_planetary_gears_architecture = wf.ModelMethod(pg_generator.GeneratorPlanetaryGearsArchitecture,  'decision_tree', name='GeneratorPlanetaryGearsArchitecture-decision_tree')

block_planetary_gears_z_number = wf.InstanciateModel(pg_generator.GeneratorPlanetaryGearsZNumber, name='GeneratorPlanetaryGearsZNumber')
generate_planetary_gears_z_number = wf.ModelMethod(pg_generator.GeneratorPlanetaryGearsZNumber,  'decision_tree', name='GeneratorPlanetaryGearsZNumber')


input_values = {}
blocks_generator_planetary_gears_architecture = []

blocks_generator_planetary_gears_architecture.extend([block_planet_structure , generate_planet_structure, block_planetary_gears_architecture, 
                generate_planetary_gears_architecture])

pipes_generator_planetary_gears_architecture = [wf.Pipe(block_planet_structure.outputs[0], generate_planet_structure.inputs[0]),
          wf.Pipe(block_planetary_gears_architecture.outputs[0], generate_planetary_gears_architecture.inputs[0]),
          wf.Pipe(generate_planet_structure.outputs[0], block_planetary_gears_architecture.inputs[0]),
         ]

workflow_generator_planetary_gears_architecture = wf.Workflow(blocks_generator_planetary_gears_architecture, pipes_generator_planetary_gears_architecture, 
                                                              generate_planetary_gears_architecture.outputs[0])



# input_values_planetar_gears_architecture = {workflow_generator_planetary_gears_architecture.index(block_planet_structure.inputs[0]): 3,
#                                             workflow_generator_planetary_gears_architecture .index(block_planet_structure.inputs[1]): 0,
#                                             workflow_generator_planetary_gears_architecture .index(block_planet_structure.inputs[2]): 2,
#                                             workflow_generator_planetary_gears_architecture .index(block_planet_structure.inputs[3]): 1,
#                                             workflow_generator_planetary_gears_architecture .index(block_planet_structure.inputs[4]):2,
#                                             workflow_generator_planetary_gears_architecture .index(block_planetary_gears_architecture.inputs[1]):[[500,550],[600,650],[300,350],[200,250]]}


# workflow_generator_planetary_gears_architecture = workflow_generator_planetary_gears_architecture.run(input_values_planetar_gears_architecture )
# workflow_generator_planetary_gears_architecture.output_value[1].plot_kinematic_graph()





blocks_generator_planetary_gears_z_number=[]
blocks_generator_planetary_gears_z_number.extend([block_planetary_gears_z_number,generate_planetary_gears_z_number])


pipes_generator_planetary_gears_z_number=[wf.Pipe(block_planetary_gears_z_number.outputs[0], generate_planetary_gears_z_number.inputs[0])]


workflow_generator_planetary_gears_z_number= wf.Workflow(blocks_generator_planetary_gears_z_number, pipes_generator_planetary_gears_z_number, generate_planet_structure.outputs[0])

# sun = pg.Planetary(36, 'Sun', 'sun')
# sun_2 = pg.Planetary(60, 'Sun', 'sun_2')
# ring = pg.Planetary(84, 'Ring', 'ring')
# planet_carrier = pg.PlanetCarrier('planet_carrier')
# planet_1 = pg.Planet( 12, 'planet_1')
# planet_2 = pg.Planet( 12, 'planet_2')
# planet_3 = pg.Planet( 16, 'planet_3')
# connection = [pg.Connection([sun, planet_1], 'GE'), 
#               pg.Connection([planet_1, planet_2], 'GE'), 
#               pg.Connection([planet_2, ring], 'GE'),
#               pg.Connection([planet_2, planet_3], 'D'), 
#               pg.Connection([planet_3, sun_2], 'GI')]

# planetary_gears_1 = pg.PlanetaryGear([sun, ring, sun_2], [planet_1, planet_2, planet_3], \
#                                      planet_carrier, connection, 'pl_1')


# input_values_planetar_gears_z_number = {workflow_generator_planetary_gears_z_number.index(block_planetary_gears_z_number.inputs[0]): planetary_gears_1,
#                                         workflow_generator_planetary_gears_z_number.index(block_planetary_gears_z_number.inputs[1]): [[500,550],[600,650],[300,350],[200,250]] ,
#                                         workflow_generator_planetary_gears_z_number.index(block_planetary_gears_z_number.inputs[2]): [7, 80],
#                                         workflow_generator_planetary_gears_z_number.index(block_planetary_gears_z_number.inputs[3]): [40,100],
#                                         workflow_generator_planetary_gears_z_number.index(block_planetary_gears_z_number.inputs[4]):3,
#                                         }


# workflow_generator_planetary_gears_z_number = workflow_generator_planetary_gears_z_number.run(input_values_planetar_gears_z_number )




block_for_each_planetary_gears_architecture= wf.ForEach(workflow_generator_planetary_gears_z_number,
                                                        block_planetary_gears_z_number.inputs[0])

block_generator_planetary_gears=[block_planet_structure , generate_planet_structure, block_planetary_gears_architecture, 
                generate_planetary_gears_architecture,block_for_each_planetary_gears_architecture]

pipes_generator_planetary_gears=[wf.Pipe(block_planet_structure.outputs[0], generate_planet_structure.inputs[0]),
                                 wf.Pipe(block_planetary_gears_architecture.outputs[0], generate_planetary_gears_architecture.inputs[0]),
                                 wf.Pipe(generate_planet_structure.outputs[0], block_planetary_gears_architecture.inputs[0]),
                                 wf.Pipe(generate_planetary_gears_architecture.outputs[0],block_for_each_planetary_gears_architecture.inputs[0])]



workflow_generator_planetary_gears=wf.Workflow(block_generator_planetary_gears,
                                              pipes_generator_planetary_gears,block_for_each_planetary_gears_architecture.outputs[0])

input_values = {workflow_generator_planetary_gears_architecture.index(block_planet_structure.inputs[0]): 3,
                workflow_generator_planetary_gears_architecture .index(block_planet_structure.inputs[1]): 0,
                workflow_generator_planetary_gears_architecture .index(block_planet_structure.inputs[2]): 2,
                workflow_generator_planetary_gears_architecture .index(block_planet_structure.inputs[3]): 1,
                workflow_generator_planetary_gears_architecture .index(block_planet_structure.inputs[4]):2,
                workflow_generator_planetary_gears_architecture .index(block_planetary_gears_architecture.inputs[1]):[[500,550],[600,650],[300,350],[200,250]],

                workflow_generator_planetary_gears.index(block_for_each_planetary_gears_architecture.inputs[1]):[[500,550],[600,650],[300,350],[200,250]] , 
                workflow_generator_planetary_gears.index(block_for_each_planetary_gears_architecture.inputs[2]):[7, 80] ,
                workflow_generator_planetary_gears.index(block_for_each_planetary_gears_architecture.inputs[3]):[40,100] ,
                workflow_generator_planetary_gears.index(block_for_each_planetary_gears_architecture.inputs[4]):3}

# a = workflow_generator_planetary_gears.to_dict()
# obj = wf.Workflow.dict_to_object(a)
# ##

workflow_generator_planetary_gears = workflow_generator_planetary_gears.run(input_values)
workflow_generator_planetary_gears.output_value[1].plot_kinematic_graph()