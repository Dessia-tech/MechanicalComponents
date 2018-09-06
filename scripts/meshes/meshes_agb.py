import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy

# 7 gears Test case with fixed modulus to 2
# definition of input data
center_distances=[(0.09, 0.14),
                  (0.09, 0.137),
                  (0.09, 0.131),
                  (0.13, 0.163),
                  (0.085, 0.115),
                  (0.12, 0.14)]
speeds = {0: [1950, 2000],
          1: [670, 700],
          2: [1900, 1950],
          3: [750, 780],
          4: [700, 740],
          5: [800, 850],
          6: [950, 980]}
connections = [(0, 1), (1, 2), (2, 3), (3, 4), (0, 5), (5, 6)]
list_rack = {0:{'name':'Catalogue_A','module':[2*1e-3,2*1e-3],
              'transverse_pressure_angle_rack':[20/180.*npy.pi,20/180.*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
rack_choices = {0:[0], 1:[0], 2:[0],3:[0], 4:[0], 5:[0], 6:[0]}
torques = [{0: 16.4, 4: 30, 5: 0, 6: 'output'}]

GA=meshes_opt.MeshAssemblyOptimizer(connections = connections, 
                                  gear_speed = speeds,
                                  center_distance = center_distances,
                                  rack_list = list_rack,
                                  torque = torques,
                                  rack_choice=rack_choices,
                                  verbose = True)

#Optimization for gear set with center-distance closed to the minimum boundary
GA.Optimize(nb_sol=35, verbose=True)
print('Number of solutions:',len(GA.solutions))
#solution=GA.solutions[-1]
#solution.SVGExport('name.txt',{6 : [0,0], 4 : [0.5,0]})
#solution.FreeCADExport('meshes3')




