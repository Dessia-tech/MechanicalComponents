#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 12:43:48 2020

@author: launay
"""



import mechanical_components.planetary_gears_generator as pg_generator
import mechanical_components.planetary_gears as pg



# volumic_mass=7800
# data_coeff_YB_Iso={'data':[[0.0,1.0029325508401201],
#                            [4.701492563229561,0.9310850480431024],
#                            [23.955224059269884,0.7609970656504502],
#                            [40.0,0.7492668574805859]
#                           ], 'x':'Linear','y':'Linear'}
# data_wholer_curve={'data':[[4.307791955971963,1.6419147590563592],
#                        [6.518240063668731,1.431665495290182],
#                        [7.989456220850952,1.4353220033111185]
#                       ], 'x':'Log','y':'Log'}
# data_gear_material={'data':[[1.313871566195314,0.7858874572688317],
#                       [1.4294457009773085,0.8802021097895326],
#                       [1.4551288380965028,0.9097910273994609]
#                      ], 'x':'Log','y':'Lèog'}
# material1=meshes.Material(volumic_mass, data_coeff_YB_Iso,
#                            data_wholer_curve, data_gear_material)
# rack=meshes.Rack(0.34,)

# meshes_1=meshes.Mesh(20,0.06,0.01,rack) 
# meshes_1.Contour()

# center_distances=[0.09713462117912072]
# # connections = [(0,1)]
# torques = {0: -16.380372067375156, 1: 'output'}
# cycles={0:1e8}
# # mesh_assembly=meshes.MeshCombination(center_distances,connections,{0:meshes_1,1:copy.copy(meshes_1)},torques,cycles)
# # volumemodel=mesh_assembly.VolumeModel({0:(0,0,0),1:(0,0.2,0.3)})   
# # volumemodel.babylonjs()                               
# # pos = vm.Point3D((0, 0, 1))
# # axis = vm.Vector3D((0,0,1))
# # radius=0.2
# # length=0.5
# # cylinder = p3d.Cylinder(pos, axis, radius, length)

# # volumemodel = vm.Contour2D(meshes_1.Contour(1) )
# # volumemodel.MPLPlot() 

sun=pg.Planetary(30,'Sun','sun')
sun_2=pg.Planetary(30,'Sun','sun_2')
sun_3=pg.Planetary(60,'Sun','sun_2')
ring= pg.Planetary(100,'Ring','ring')
planet_carrier= pg.PlanetCarrier('planet_carrier')
planet_1=pg.Planet(12,'planet_1')

sun.speed_input=[12,100]
sun_2.speed_input=[10,100]
ring.speed_input=[10,100]
planet_carrier.speed_input=[10,100]

sun.torque_input=[10,100]
sun_2.torque_input=[12,100]
ring.torque_input=[13,100]
planet_carrier.torque_input=[10,100]
list1= [5,1,3]
list2= [1,5,3]

planet_2=pg.Planet(30,'planet_2')
planet_3=pg.Planet(13,'planet_3')
planet_4=pg.Planet(5,'planet_4')
planet_5=pg.Planet(5,'planet_5')
planet_6=pg.Planet(5,'planet_5')
planet_7=pg.Planet(5,'planet_5')
connections=[pg.Connection([sun,planet_1],'GE'),pg.Connection([planet_1,planet_2],'GE'),pg.Connection([planet_2,ring],'GE'),
              pg.Connection([planet_2,planet_3],'D'),
              pg.Connection([planet_3,sun_2],'GI')]

connections3=[pg.Connection([sun,planet_1],'GE'),pg.Connection([planet_1,planet_2],'GE'),pg.Connection([planet_2,ring],'GE'),
              pg.Connection([planet_2,planet_3],'D'),pg.Connection([planet_5,sun_2],'GI'),
              pg.Connection([planet_3,planet_4],'GI'),
              pg.Connection([planet_4,planet_5],'GI')]

connections_2=[pg.Connection([sun,planet_1],'GE'),pg.Connection([planet_1,planet_2],'GE'), pg.Connection([planet_2,ring],'GE')]

# generator_planetary_gear=pg_generator.GeneratorPlanetaryGears(4,[[(0,0),(165,175),(-1000,1000),(9000,10000)],[(25,35),(0,0),(-1000,1000),(9000,10000)]],0,10)

# print(generator_planetary_gear.speed_conversion())

# plt.figure()
# sun.volume_model(4.0,(0,0),0,0.1)
# sun_2.volume_model(4.0,(0,0),0,0.1)




planetary_gears_1= pg.PlanetaryGear([sun,ring,sun_2], [planet_1,planet_2,planet_3], planet_carrier,connections,3,'pl_1')
# planetary_gears_2= pg.PlanetaryGear([sun,ring,sun_2,], [planet_1,planet_2,planet_3,planet_4,planet_5], planet_carrier,connections3,3,'pl_1')





# vmp.plot(planetary_gears_1.plot_data())
# planetary_gears_2= pg.PlanetaryGear([sun,ring], [planet_1,planet_2], planet_carrier,connections_2,3,'pl_1')

# input_torque_and_composant={sun:100,ring:100}
# input_speed={sun:150,ring:100}
# result_torque=planetary_gears_1.torque_solve(input_torque_and_composant)
# result_speed=planetary_gears_1.speed_solve(input_speed)
# print(result_torque)
# print(result_speed)
# component=[sun,ring,sun_2,planet_1,planet_2,planet_3,planet_carrier]
# puissance=[]
# for i,element in enumerate(component):
#     puissance.append(result_speed[element]*result_torque[element])

generatorgeometry=pg_generator.GeneratorPlanetaryGearsGeometry(planetary_gears_1,3,1,10,200)
generatorgeometry.verification()
# planetary1=generatorgeometry.planetary_gear.planetaries[0]
# planet=generatorgeometry.planetary_gear.planets[0]
# module=0.45
# Z_1=-100
# Z_2=25
# Z_3=50
# center_distances=[[0,abs((Z_2+Z_1)*module/2)],[0,abs((Z_2+Z_3)*module/2)]]

# list_rack = {0:{'name':'Catalogue_A','module':[module,module],
#               'transverse_pressure_angle_rack':[20/180.*npy.pi,20/180.*npy.pi],
#               'coeff_gear_addendum':[0.85,0.85],'coeff_gear_dedendum':[1,1],
#               'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.4,0.4]}}

# rack=meshes_opt.RackOpti(module=[module,module],transverse_pressure_angle=[20/180.*npy.pi,20/180.*npy.pi],
#              coeff_gear_addendum=[1,1],coeff_gear_dedendum=[0.01,2],coeff_root_radius=[0.01,2],
#              coeff_circular_tooth_thickness=[0.01,2])
# list_rack = {0:rack}

# rack_choices = {0:0, 1:0 , 2:0}
# db2=m.cos(20/180.*npy.pi)*module*Z_2

# db1=m.cos(20/180.*npy.pi)*module*abs(Z_1)
# db3=m.cos(20/180.*npy.pi)*module*Z_3
# print(center_distances)
# torques = {0: 'output', 1: 0, 2:100}
# cycles = {0: 1272321481513.054 }
# material={0:hardened_alloy_steel}
# transverse_pressure_angle={0: [20/180.*npy.pi-0.1, 20/180.*npy.pi],1: [20/180.*npy.pi-0.1, 20/180.*npy.pi]}
# db=[[db1-0.2,db1],[db2-0.2,db2],[db3-0.2,db3]]
# coefficient_profile_shift=[[0.01,0.01],[0.01,0.01],[0.01,0.01]]
# d=mg.ContinuousMeshesAssemblyOptimizer(Z={0:Z_1,1:Z_2,2:Z_3},center_distances=center_distances,connections=[[(0,1)],[(1,2)]],rigid_links=[],
#                                         transverse_pressure_angle=transverse_pressure_angle,rack_list=list_rack,rack_choice=rack_choices,material=material,
#                                         external_torques=torques,cycles=cycles,safety_factor=1,db=db,coefficient_profile_shift=coefficient_profile_shift)
# d.Optimize(verbose=True)
# solution=d.solutions[0]
# # plt.figure()
# Z5=solution.mesh_combinations[0]



# debut=time.time()
# generatorgeometry.verification_2()
# # generatorgeometry.optimize_max()
# generatorgeometry.optimize_min_recirculation_2()
# generatorgeometry.planetary_gear.babylonjs()
# planetary_gears_1.babylonjs()
# print(planetary_gears_1.recirculation_power())

# planetary_gear_result=pg.PlanetaryGearResult(planetary_gears_1,generatorgeometry.position_min_max)
# planetary_gear_result.planetary_gear.babylonjs()
# planetary_gear_result.babylonjs()

generatorgeometry.optimize_min()
# planetary_gear_result=pg.PlanetaryGearResult(planetary_gears_2,generatorgeometry.position_min_max)
# planetary_gear_result.update_geometry()
# print(planetary_gear_result.planetary_gear.recirculation_power())
# fin=time.time()
# print(debut-fin)
# planetary_gear_result.update_torque({sun:[10,100]})
# planetary_gear_result.update_torque_max()
# print(planetary_gears_1.mech)
# c = Client(api_url = 'http://localhost:5000')
# r = c.create_object_from_python_object(planetary_gears_1)

# input_torque_and_composant={planetary_gears_1.planetaries[2]:2,planetary_gears_1.planetaries[1]:-5}
# link=planetary_gears_1.torque_resolution_PFS(input_torque_and_composant)
# print(planetary_gear_result.planetary_gear.recirculation_power())


# print(link)
# result_torque=planetary_gears_1.torque_solve(input_torque_and_composant)
# print(result_torque)
# input_speed_and_composant={planetary_gears_1.planetaries[0]:300,planetary_gears_1.planetaries[2]:200}
# result_speed=planetary_gears_1.speed_solve(input_speed_and_composant)
# print(result_speed)

# vmp.plot_d3(planetary_gears_1.plot_data())
# c = Client(api_url = 'http://localhost:5000')
# r = c.create_object_from_python_object(planetary_gears_1)

# li_box=planetary_gears_1.volmdlr_primitives()
# planetary_gears_1.babylonjs()
# print([sun,ring,sun_2])
# torque_solution=planetary_gears_1.torque_solve({sun:0,planet_carrier:500})
# speed_solution=planetary_gears_1.speed_èsolve({sun:200,planet_carrier:500})
# print(torque_solution)
# print(speed_solution)
# debut=time.time()
# Generator_planet_structure=pg_generator.GeneratorPlanetsStructure(3,0,2,1,2)
# list_planet_structure=Generator_planet_structure.decision_tree()
# # print(len(list_planet_structure))
# # for planet_structure in list_planet_structure:
# #     planet_structure.plot_kinematic_graph()
# # c = Client(api_url = 'http://localhost:5000')
# # r = c.create_object_from_python_object(Generator_planet_structure)
# Generator_planetarie_gears=pg_generator.GeneratorPlanetaryGearsArchitecture(list_planet_structure,[[500,550],[600,650],[300,350],[200,250]])
# list_planetary_gears=Generator_planetarie_gears.decision_tree()

# print(len(list_planetary_gears))
# for planetary_gears in list_planetary_gears:
#     planetary_gears.plot_kinematic_graph()
#     print(planetary_gears)

# for i in range(len(list_planetary_gears)):
# list_solution=[]
# for planetary_gear in list_planetary_gears:
# Generator_planetarie_gear_z=pg_generator.GeneratorPlanetaryGearsZNumber(list_planetary_gears[0],[[500,505],[610,615],[310,315],[380,385]],[7,80],[40,100],3)
# list_solution=Generator_planetarie_gear_z.decision_tree()
# a=list_solution[0].speed_min_max_planets()


# for planetary_gear in list_solution:
#     generatorgeometry=pg_generator.GeneratorPlanetaryGearsGeometry(planetary_gear,3,10,100)
#     print(generatorgeometry.verification())

# Z_planetary=[]
# for planetary_gears in list_solution:
#     Z=[]
#     for planetary in planetary_gears.planetaries:
#        Z.append(planetary.Z)
#     if not Z  in Z_planetary:
#         Z_planetary.append(Z)
    

# print(len(Z_planetary))
            
# print(len(list_solution))
# number=0
# solution=list_solution[7]
# for i,solution in enumerate(list_solution):
    
#     list_range_1=solution.speed_range(solution.planetaries[0],solution.planetaries[1],0.01)
#     print(list_range_1)
#     list_range_2=solution.speed_range(solution.planetaries[0],solution.planetaries[2],0.01)
#     print(list_range_2)
#     list_range_3=solution.speed_range(solution.planetaries[1],solution.planetaries[2],0.01)
#     print(list_range_3)
#     list_range_4=solution.speed_range(solution.planetaries[1],solution.planet_carrier,0.01)
#     print(list_range_4)
#     list_range_5=solution.speed_range(solution.planetaries[2],solution.planet_carrier,0.01)
#     print(list_range_5)
#     if not list_range_4 or not list_range_5 or not list_range_1 or not list_range_2 or not list_range_3:
        
#         print(i)
#         number+=1
        
# print(number)
# list_range_2=list_solution[7].speed_range(list_solution[7].planetaries[0],list_solution[7].planet_carrier,[])
# list_range_1=list_solution[7].speed_range(list_solution[7].planetaries[0],list_solution[7].planetaries[1],[])
# print(list_range_1)
# print(list_range_2)
# print(list_solution[7].speed_solve({list_solution[7].planetaries[1]:list_range_1[list_solution[7].planetaries[1]][0],list_solution[7].planetaries[0]:list_range_2[list_solution[7].planetaries[0]][1]}))

# fin=time.time()
# print(debut-fin)
# for planetary_gear in list_planetary_gear:
#     planetary_gear.plot_kinematic_graph()


# list_pos=[]
# list_solution=Generator.decision_tree_architecture(3,0,2,1)


# # list_intervalle=pg.intervalle_fonction_test([-300,300],[-700,700],[-1000,1100],[],1,1)
# # print(len(list_solution))
# for i,solution in enumerate(list_solution): 
#     print(i)
# # # # list_solution[24].plot()


#     solution_2=Generator.decision_tree_z_number(list_solution[i],[7,80],[40,100],3)
    # if solution_2:
    #      break
# for solution in  solution_2:
#     list_range_1=solution.speed_range(solution.planetaries[0],solution.planetaries[2],[])
#     print(list_range_1)
#     # print(solution.speed_solve({solution.planetaries[1]:306.15,solution.planetaries[2]:251.58}))
#     list_range_2=solution.speed_range(solution.planetaries[0],solution.planet_carrier,[])
#     print(list_range_2)
# print(solution[0])
# list_range_1=solution[0].speed_range(solution[0].planetaries[1],solution[0].planetaries[0],[])
# print(list_range_1)
# list_range_2=solution[0].speed_range(solution[0].planetaries[2],solution[0].planetaries[0],[])
# print(list_range_2)
# print(solution[0].speed_solve({solution[0].planetaries[1]:303,solution[0].planetaries[0]:753}))
# for solution in list_solution_2:
#     print(solution.speed_solve({solution.planetaries[1]:400,solution.planet_carrier:-400}))
#     if solution.speed_solve({solution.planetaries[1]:400,solution.planet_carrier:-400})[0]>680:
        # break
#     # if i==20:
#     #     break
#     solution.plot_kinematic_graph(0.1,1,2,0.2,0.5,2,2,10)

  
# list_solution[14].plot()

    
    # plt.savefig('bifurcation 0_'+str(i))
# list_pos=[]
# list_previous=[]
# # Generator.list_possibilities_planets_by_branch_step_1([0,0,0],list_pos,6,0,3,0,1)
# # print(list_pos)
# global_architecture,number_branch=Generator.list_possibilities_planets_by_branch_step_2([0,2,0],5,2,1)
# print(global_architecture)
# # print(global_architecture,number_branch)
# list_branch=[]
# list_planet_type=[]
# for i in range (number_branch):
#     list_branch.append(0)
#     list_planet_type.append('Simple')
# list_connexion=[]
# for i in range(len(global_architecture)):
#     Generator.list_possibilities_architecture_planet(1,number_branch,global_architecture[i],list_branch,list_pos,[],list_connexion)
# print(list_pos,list_connexion)
#Generator.number_possibilities_planets_type(0,[],['Simple','Simple','Simple','Simple','Simple','Simple','Simple','Simple','Simple','Simple'],10)
# volume=vm.VolumeModel(planetary_gears_1.volume_plot([0,0],0,[0.02,0.01],0.1,0.5,0.05))
# volume.babylonjs()




# planetary_gears_2= pg.PlanetaryGears('pl_2', [sun,ring], [planet_1], planet_carrier)
# planetary_gears_3= pg.PlanetaryGears('pl_3', sun, ring, [planet_1], planet_carrier)
# planetary_gears_4= pg.PlanetaryGears('pl_4', sun, ring, [planet_1], planet_carrier)
# assembly_planetary_gear=pg.AssemblyPlanetaryGears('assembly_planetary_gear', 
#                                                   [planetary_gears_1,planetary_gears_2,planetary_gears_3,planetary_gears_4], 
#                                                   [[[sun,planetary_gears_2],[sun,planetary_gears_3]],
#                                                     [[planet_carrier,planetary_gears_3],[ring,planetary_gears_4]],
#                                                     [[planet_carrier,planetary_gears_1],[planet_carrier,planetary_gears_4]],
#                                                     [[sun,planetary_gears_1],[ring,planetary_gears_2]],
#                                                     [[ring,planetary_gears_1],[sun,planetary_gears_4]]])

# print(planetary_gears_1)

# planet_carrier=pg.PlanetCarrier('PlanetCarrier')
# planetary_1=pg.Planetary('Planetary_1',7,'Ring')
# planetary_2=pg.Planetary('Planetary_2',7,'Ring')
# list_element={'Planet_Carrier': planet_carrier, 'Planetary_1' : planetary_1,'Planetary_2':planetary_2}
# planets=[]
# for i in range(2):
#             planets.append(pg.Planet('Planet'+str(i),'Double',7)) 
                                                   
# pg.test_ratio_max_ratio_min([2,69,60,42,9],planetary_1,planetary_2,planet_carrier,planets,{'Planetary_1':540,'Planet_Carrier':0},'Planetary_2',200,[7 , 80],3)




# planetary_gears_1.solve(500, ring,planet_3)
# assembly_planetary_gear.plot()
# print(assembly_planetary_gear.system_equations()[0])print(node)
                        
#solution,system_matrix=assembly_planetary_gear.solve({planet_carrier:[500,planetary_gears_2] ,ring:[0,planetary_gears_3],sun:[0,planetary_gears_2]},)

#print(assembly_planetary_gear.solve(500,planet_carrier,planetary_gears_2,[ring,sun],[planetary_gears_3,planetary_gears_2]))

# Generator=pg.GeneratorPlanetaryGears([[200,250],[100,150], [300,350],[200,300],[200,350],[11,12],[13,15]],0.1,2,45)

# Generator.decission_tree_()



# solutions=[]
# solutions = pg.decision_tree_planetary_gears({'Planetary_1':500,'Planet_Carrier':0},'Planetary_2',230,2,0.1,3,[7 , 80],[0.1, 1],0.1)

# x = vm.Vector3D((1,0,0))
# y = vm.Vector3D((0,1,0))
# z = vm.Vector3D(npy.cross(x.vector, y.vector))



# Gears3D={0:meshes_1.Contour(3)}
# export=[]
# center_2=(0,0)
# center = vm.Point2D(center_2)
# model_trans =Gears3D[0][0].Translation(center)
# model_trans_rot = model_trans.Rotation(center, 0.1)
# Gears3D_Rotate=[model_trans_rot]
# list_gear=[Gears3D[0]]
# list_center=[center_2]
# list_rot=[-1]
# export=[]
# for (i,center,k) in zip(list_gear,list_center,list_rot):
#             model_export=[]
            
#             for m in i:
                
#                 center = vm.Point2D(center)
#                 model_trans = m.Translation(center)
#                 model_trans_rot = model_trans.Rotation(center, k)
#                 model_export.append(model_trans_rot)
#             export.append(model_export)
# Gears3D_Rotate=export

# vect_x = -0.5*0.01*x 
# extrusion_vector1 = 0.01*x
# C1=vm.Contour2D(Gears3D_Rotate[0])
# i=vm.Vector3D(vect_x)
# t1=p3d.ExtrudedProfile(vm.Vector3D(vect_x), y, z, C1, [], vm.Vector3D(extrusion_vector1))
# modle=vm.VolumeModel([t1],'d')
# modle.babylonjs()

