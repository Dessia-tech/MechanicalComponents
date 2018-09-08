import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy

# 7 gears Test case with fixed modulus to 2
# definition of input data
speeds = {0: [1973.4035805488818, 1982.1057009682193],
               1: [512.9532091805967, 536.0308586733844],
               2: [1917.3652941203281, 2003.6271273936354],
               3: [589.7265856084168, 616.2582520366075],
               4: [705.2891418162362, 737.0199416525495],
               5: [976.5850066096972, 1020.5213463753686],
               6: [930.487760604817, 972.3501956422525]}

torques = [{0: 'output', 1: 0, 2: -410.0, 3: 0, 4: 0, 5: 0, 6: 0},
           {0: 'output', 1: 0, 2: 0.2117260970058401, 3: 0,
            4: 81.66578027368121, 5: 0, 6: 5.63051765539608}]

center_distances = [(0.13142967869478409, 0.1714296786947841),
                   (0.13098751234708322, 0.17098751234708323),
                   (0.11655203593568475, 0.15655203593568476),
                   (0.1651650012440174, 0.2051650012440174),
                   (0.08147365232706552, 0.12147365232706553),
                   (0.11133464347238507, 0.15133464347238507)]


connections = [(0, 1), (1, 2), (2, 3), (3, 4), (0, 5), (5, 6)]




list_rack = {0:{'name':'Catalogue_A','module':[2*1e-3,2*1e-3],
              'transverse_pressure_angle_rack':[20/180.*npy.pi,20/180.*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}

rack_choices = {0:[0], 1:[0], 2:[0],3:[0], 4:[0], 5:[0], 6:[0]}

GA=meshes_opt.MeshAssemblyOptimizer(connections = connections, 
                                  gear_speed = speeds,
                                  center_distance = center_distances,
                                  rack_list = list_rack,
                                  torque = torques,
                                  rack_choice=rack_choices,
                                  verbose = True)

#Optimization for gear set with center-distance closed to the minimum boundary
GA.Optimize(nb_sol=1, verbose=True)
print('Number of solutions:',len(GA.solutions))
solution=GA.solutions[-1]
#solution.SVGExport('name.txt',{6 : [0,0], 4 : [0.5,0]})
solution.FreeCADExport('meshes_agb2')



