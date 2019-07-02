import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy

# 7 gears Test case with fixed modulus to 2
# definition of input data

connections = [[(0, 1), (2, 3), (0, 3)]]

speeds = {0: [425.3901758131878, 427.26601944275086],
          1: [416.7824497926327, 435.53339845406],
          2: [727.2667514172009, 759.9863189656231],
          3: [530.5572830702263, 554.4269358872916]}

torques=[{0: 'output', 1: 0, 2: 0, 3: 0},
          {0: 'output', 1: 193.40347514964495, 2: 229.0154116255149, 3: 0}]

center_distances = [(0.23005304306554464, 0.2700530430655446),
                   (0.1563104344511389, 0.1963104344511389),
                   (0.2053813701229723, 0.24538137012297231)]

list_rack = {0:{'name':'Catalogue_A','module':[2*1e-3,2*1e-3],
              'transverse_pressure_angle_rack':[20/180.*npy.pi,20/180.*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}

rack_choices = {0:[0], 1:[0], 2:[0],3:[0], 4:[0], 5:[0], 6:[0]}

GA=meshes_opt.MeshAssemblyOptimizer(connections = connections, 
                                  gear_speeds = speeds,
                                  center_distances = center_distances,
                                  rack_list = list_rack,
                                  torques = torques,
                                  cycles = {0:1e8},
                                  rack_choice=rack_choices,
                                  verbose = True)

#Optimization for gear set with center-distance closed to the minimum boundary
GA.Optimize(nb_sol=3, verbose=True)
print('Number of solutions:',len(GA.solutions))
solution=GA.solutions[-1]
#solution.SVGExport('name.txt',{6 : [0,0], 4 : [0.5,0]})
solution.FreeCADExport('meshes_agb2')



