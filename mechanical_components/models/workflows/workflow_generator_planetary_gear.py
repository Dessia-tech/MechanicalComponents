#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 14:13:57 2020

@author: launay
"""
import mechanical_components.planetary_gears_generator as pg_generator
import mechanical_components.planetary_gears as pg
import dessia_common.workflow as wf
from dessia_common.workflow import Filter
import mechanical_components.planetary_gears as pg
from dessia_api_client import Client
from dessia_common.vectored_objects import Catalog, Objective, ParetoSettings, ObjectiveSettings, from_csv

import os
block_planet_structure = wf.InstanciateModel(pg_generator.GeneratorPlanetsStructure, name='GeneratorPlanetsStructure')
generate_planet_structure = wf.ModelMethod(pg_generator.GeneratorPlanetsStructure,  'decision_tree', name='GeneratorPlanetsStructure-decision_tree')

block_planetary_gears_architecture = wf.InstanciateModel(pg_generator.GeneratorPlanetaryGearsArchitecture, name='GeneratorPlanetaryGearsArchitecture')
generate_planetary_gears_architecture = wf.ModelMethod(pg_generator.GeneratorPlanetaryGearsArchitecture,  'decision_tree', name='GeneratorPlanetaryGearsArchitecture-decision_tree')

block_planetary_gears_z_number = wf.InstanciateModel(pg_generator.GeneratorPlanetaryGearsZNumber, name='GeneratorPlanetaryGearsZNumber')
generate_planetary_gears_z_number = wf.ModelMethod(pg_generator.GeneratorPlanetaryGearsZNumber,  'decision_tree', name='GeneratorPlanetaryGearsZNumber')

block_planetary_gears_geometry = wf.InstanciateModel(pg_generator.GeneratorPlanetaryGearsGeometry, name='GeneratorPlanetaryGearsGeometry')
generate_planetary_gears_geometry = wf.ModelMethod(pg_generator.GeneratorPlanetaryGearsGeometry,  'verification', name='GeneratorPlanetaryGearsGeometry')
optimize_min_planetary_gears_geometry=  wf.ModelMethod(pg_generator.GeneratorPlanetaryGearsGeometry,  'optimize_min', name='OptimizePlanetaryGearsGeometryMin')
optimize_max_planetary_gears_geometry=  wf.ModelMethod(pg_generator.GeneratorPlanetaryGearsGeometry,  'optimize_max', name='OptimizePlanetaryGearsGeometryMax')


block_planetary_gear_result=wf.InstanciateModel(pg.PlanetaryGearResult, name='PlanetaryGearResult')

# workflow_planetary_gear_result= wf.Workflow([block_planetary_gear_result],[],block_planetary_gear_result.outputs[0])

# block_workflow_planetary_gear_result=wf.WorkflowBlock(workflow_planetary_gear_result)

# block_for_each_planetary_gear_result= wf.ForEach(block_workflow_planetary_gear_result,
#                                                         block_workflow_planetary_gear_result.inputs[0])



block_solution_sort = wf.InstanciateModel(pg_generator.SolutionSort, name='SolutionSort')
solution_sort=wf.ModelMethod(pg_generator.SolutionSort,  'solution_sort', name='solution_sort')
list_planetary_gears=wf.ModelMethod(pg_generator.SolutionSort,  'list_planetary_gears', name='list_planetary_gears')




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


workflow_generator_planetary_gears_z_number= wf.Workflow(blocks_generator_planetary_gears_z_number, pipes_generator_planetary_gears_z_number, generate_planetary_gears_z_number.outputs[0])

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

blocks_generator_planetary_gears_geometry=[]
blocks_generator_planetary_gears_geometry.extend([block_planetary_gears_geometry,generate_planetary_gears_geometry,optimize_min_planetary_gears_geometry,
                                                  optimize_max_planetary_gears_geometry,block_planetary_gear_result])


pipes_generator_planetary_gears_geometry=[wf.Pipe(block_planetary_gears_geometry.outputs[0], generate_planetary_gears_geometry.inputs[0]),
                                          wf.Pipe(generate_planetary_gears_geometry.outputs[1], optimize_min_planetary_gears_geometry.inputs[0]),
                                          wf.Pipe(optimize_min_planetary_gears_geometry.outputs[1], optimize_max_planetary_gears_geometry.inputs[0]),
                                          wf.Pipe(generate_planetary_gears_geometry.outputs[0], block_planetary_gear_result.inputs[0]),
                                          wf.Pipe(optimize_max_planetary_gears_geometry.outputs[0], block_planetary_gear_result.inputs[1])]

workflow_generator_planetary_gears_geometry = wf.Workflow(blocks_generator_planetary_gears_geometry, pipes_generator_planetary_gears_geometry, block_planetary_gear_result.outputs[0])

block_workflow_generator_planetary_gears_geometry=wf.WorkflowBlock(workflow_generator_planetary_gears_geometry)

 

block_for_each_planetary_gears_geometry= wf.ForEach(block_workflow_generator_planetary_gears_geometry,0)


block_workflow_generator_planetary_gears_z_number=wf.WorkflowBlock(workflow_generator_planetary_gears_z_number)





block_for_each_planetary_gears_z_number= wf.ForEach(block_workflow_generator_planetary_gears_z_number,
                                                        0)

filters = [
          {'attribute' : 'sum_Z_planetary', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'min_Z_planetary', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'max_Z_planetary', 'operator' : 'gt', 'bound' : -100},
          {'attribute' : 'D_train', 'operator' : 'gt', 'bound' : -100},
           {'attribute' : 'speed_max_planet', 'operator' : 'gt', 'bound' : -100}]
          
filter_analyze= wf.Filter(filters)

list_attribute=['sum_Z_planetary','min_Z_planetary','max_Z_planetary','D_train','speed_max_planet', 'recycle_power']

block_parallel_plot=wf.ParallelPlot(list_attribute, 4 ,name = 'Parallel_Plot')
minimized_attributes = {'sum_Z_planetary':False,'D_train': True,'min_Z_planetary':True,'max_Z_planetary':True,'speed_max_planet': True}

pareto_settings = ParetoSettings(minimized_attributes=minimized_attributes,
                                  enabled=True)


list_attribute2=['number_meshing','number_double', 'number_group_solution_planet_structure']
display_architecture_planets_structure=wf.ParallelPlot(list_attribute2,1,name='Display planet structure')

list_attribute3=['number_ring','number_sun','number_group_solution_planet_structure','number_group_solution_architecture']
display_architecture_planetary_gear=wf.ParallelPlot(list_attribute3,2,name='Display Architecture planetary gear')

list_attribute4=['sum_Z_planetary','min_Z_planetary','max_Z_planetary','speed_max_planet','number_group_solution_architecture']
display_planetary_gear_Z=wf.ParallelPlot(list_attribute4,3,'Display planetary gear Z')

block_generator_planetary_gears=[block_planet_structure , generate_planet_structure, block_planetary_gears_architecture, 
                generate_planetary_gears_architecture,block_for_each_planetary_gears_z_number,block_solution_sort,solution_sort,
                block_for_each_planetary_gears_geometry,filter_analyze,block_parallel_plot,display_architecture_planets_structure,
                display_architecture_planetary_gear,list_planetary_gears,display_planetary_gear_Z]



pipes_generator_planetary_gears=[wf.Pipe(block_planet_structure.outputs[0], generate_planet_structure.inputs[0]),
                                  wf.Pipe(block_planetary_gears_architecture.outputs[0], generate_planetary_gears_architecture.inputs[0]),
                                  wf.Pipe(generate_planet_structure.outputs[0], block_planetary_gears_architecture.inputs[0]),
                                  wf.Pipe(generate_planet_structure.outputs[0], display_architecture_planets_structure.inputs[0]),
                                  wf.Pipe(generate_planetary_gears_architecture.outputs[0],block_for_each_planetary_gears_z_number.inputs[0]),
                                  wf.Pipe(generate_planetary_gears_architecture.outputs[0],display_architecture_planetary_gear.inputs[0]),
                                  wf.Pipe(block_for_each_planetary_gears_z_number.outputs[0],block_solution_sort.inputs[0]),
                                  wf.Pipe(block_solution_sort.outputs[0],solution_sort.inputs[0]),
                                  wf.Pipe(block_solution_sort.outputs[0],list_planetary_gears.inputs[0]),
                                  wf.Pipe(list_planetary_gears.outputs[0],display_planetary_gear_Z.inputs[0]),
                                  wf.Pipe(solution_sort.outputs[0],block_for_each_planetary_gears_geometry.inputs[0]),
                                  wf.Pipe(block_for_each_planetary_gears_geometry.outputs[0],filter_analyze.inputs[0]),
                                  wf.Pipe(filter_analyze.outputs[0],block_parallel_plot.inputs[0])]
    
                                  # wf.Pipe(filter_analyze.outputs[0],block_catalog.inputs[0]),
                                  # wf.Pipe(block_catalog.outputs[0],catalog.inputs[0])]
                                  


workflow_generator_planetary_gears=wf.Workflow(block_generator_planetary_gears,
                                              pipes_generator_planetary_gears,filter_analyze.outputs[0])



input_values = {workflow_generator_planetary_gears_architecture.index(block_planet_structure.inputs[0]): 3,
                workflow_generator_planetary_gears_architecture .index(block_planet_structure.inputs[1]): 0,
                workflow_generator_planetary_gears_architecture .index(block_planet_structure.inputs[2]): 2,
                workflow_generator_planetary_gears_architecture .index(block_planet_structure.inputs[3]): 1,
                workflow_generator_planetary_gears_architecture .index(block_planet_structure.inputs[4]):2,
                workflow_generator_planetary_gears_architecture .index(block_planetary_gears_architecture.inputs[1]):[[500,501],[610,611],[310,315],[380,385]],

                workflow_generator_planetary_gears.index(block_for_each_planetary_gears_z_number.inputs[1]):[[500,505],[610,615],[310,315],[380,385]] , 
                workflow_generator_planetary_gears.index(block_for_each_planetary_gears_z_number.inputs[2]):[[10,60],[-70,-10],[-80,-10],[100,200]],
                workflow_generator_planetary_gears.index(block_for_each_planetary_gears_z_number.inputs[3]):[7, 80] ,
                workflow_generator_planetary_gears.index(block_for_each_planetary_gears_z_number.inputs[4]):[40,100] ,
                workflow_generator_planetary_gears.index(block_for_each_planetary_gears_z_number.inputs[5]):3,
                workflow_generator_planetary_gears.index(block_for_each_planetary_gears_z_number.inputs[6]):1000,
                
                workflow_generator_planetary_gears.index(block_for_each_planetary_gears_geometry.inputs[1]):3,
                workflow_generator_planetary_gears.index(block_for_each_planetary_gears_geometry.inputs[2]):1,
                workflow_generator_planetary_gears.index(block_for_each_planetary_gears_geometry.inputs[3]):10,
                workflow_generator_planetary_gears.index(block_for_each_planetary_gears_geometry.inputs[4]):250}

                
                # workflow_generator_planetary_gears.index(block_parallel_plot.inputs[1]):pareto_settings

# a = workflow_generator_planetary_gears.to_dict()
# obj = wf.Workflow.dict_to_object(a)
# ##
workflow_generator_run = workflow_generator_planetary_gears.run(input_values)

# variables=['name','sum_Z_planetary','d_min','sum_speed_planetary','speed_planet_carrer','min_Z_planetary','max_Z_planetary']
# array=[]

# choice_args = ['sum_Z_planetary','d_min','sum_speed_planetary','speed_planet_carrer','min_Z_planetary','max_Z_planetary']




# for i,planetary_gear in enumerate(workflow_generator_run.output_value):

    
#     array.append((planetary_gear.name+str(i),planetary_gear.sum_Z_planetary,planetary_gear.d_min,planetary_gear.sum_speed_planetary,planetary_gear.speed_planet_carrer,planetary_gear.min_Z_planetary,planetary_gear.max_Z_planetary)) 

# print(array)
# print(variables)
# catalog = Catalog(array=array,
#                   variables=variables,
#                   choice_variables=choice_args,
#                   objectives=[],
#                   pareto_settings=pareto_settings,
#                   name='Planetary_gears')

planetary_gear_1=workflow_generator_run.output_value[0]
# planetary_gear_1.planetary_gear.mesh_generation()


c = Client(api_url = 'https://api.platform-dev.dessia.tech')
r = c.create_object_from_python_object(workflow_generator_run)

# r2= c.create_object_from_python_object(catalog)

# <<<<<<< HEAD
# workflow_generator_planetary_gears = workflow_generator_planetary_gears.run(input_values)
# =======
# >>>>>>> planatery_gears
# workflow_generator_planetary_gears.output_value[1].plot_kinematic_graph()