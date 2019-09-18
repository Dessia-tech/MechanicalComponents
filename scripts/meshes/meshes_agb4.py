import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy

# 7 gears Test case with fixed modulus to 2
# definition of input data

connections = [[(0, 1)], [(1, 3)], [(0, 2)]]

speeds = {0: (1501.3273041622153, 1562.6059696382242),
          1: (1708.7167736682445, 1778.4603154506217),
          2: (1467.5384133029263, 1527.4379403765154),
          3: (1722.7733871400499, 1793.090668247807)}

torques = {0: 'output', 1: -450, 2: 0, 3: 0}

center_distances = [(0.2583900996409985, 0.2871001107122206),
                    (0.2407104031685444, 0.2674560035206049),
                    (0.27825055112729036, 0.3091672790303226)]

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
                                    cycles = {0: 1526097819603.9133},
                                    rack_choice=rack_choices,
                                    verbose = True,
                                    )

#Optimization for gear set with center-distance closed to the minimum boundary
GA.Optimize(nb_sol=3, verbose=True)
print('Number of solutions:',len(GA.solutions))
solution=GA.solutions[-1]
#solution.SVGExport('name.txt',{6 : [0,0], 4 : [0.5,0]})
solution.FreeCADExport('meshes_agb4')

