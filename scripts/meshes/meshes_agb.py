import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy

# 7 gears Test case with fixed modulus to 2
# definition of input data
connections = [[(0, 1)], [(1, 2)], [(2, 3)], [(3, 4)], [(0, 5)], [(5, 6)]]

rigid_links = []

gear_speeds = {0: (1878.1453579221634, 1974.460504482274),
               1: (449.8807309958231, 472.95153771355757),
               2: (1875.9428646237072, 1972.145062809538),
               3: (480.3124107238717, 504.94381640201897),
               4: (700.6720400543617, 736.6039395443289),
               5: (725.5314102562529, 762.738149243753),
               6: (865.8158304033849, 910.216642218943)}

center_distances = [(0.11134984458664793, 0.1293457790652981),
                    (0.11137510812362536, 0.12937512559815068),
                    (0.10568352491185784, 0.1227636905541783),
                    (0.14181884668124436, 0.164739064326698),
                    (0.07722003844178454, 0.0897000446545982),
                    (0.10237916425095688, 0.11892529180666707)]

cycles = {0: 1272321481513.054}
torques = {0: 'output', 1: 0, 2: -410.0, 3: 0, 4: 0, 5: 0, 6: 0}


list_rack = {0:{'name':'Catalogue_A','module':[2*1e-3,2*1e-3],
              'transverse_pressure_angle_rack':[20/180.*npy.pi,20/180.*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
rack_choices = {0:[0], 1:[0], 2:[0],3:[0], 4:[0], 5:[0], 6:[0]}


GA=meshes_opt.MeshAssemblyOptimizer(connections = connections, 
                                  gear_speeds = gear_speeds,
                                  center_distances = center_distances,
                                  cycles = cycles,
                                  rack_list = list_rack,
                                  torques = torques,
                                  rack_choice=rack_choices,
                                  verbose = True)

#Optimization for gear set with center-distance closed to the minimum boundary
GA.Optimize(nb_sol=35, verbose=True)
print('Number of solutions:',len(GA.solutions))
#solution=GA.solutions[-1]
#solution.SVGExport('name.txt',{6 : [0,0], 4 : [0.5,0]})
#solution.FreeCADExport('meshes3')


