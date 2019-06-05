import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy

# 7 gears Test case with fixed modulus to 2
# definition of input data

connections =  [[(0, 1)], [(0, 3)], [(2, 3)]] 

gear_speeds =  {0: (415.3740980732496, 432.32814289256595),
               1: (1774.8921913060299, 1847.336770543011),
               2: (754.2994180129452, 785.0871493604124),
               3: (461.6759910719359, 480.51990907487215)}
center_distances_bounds = [(0.09504015622236828, 0.1056001735804092),
                           (0.14630841589872876, 0.16256490655414307),
                           (0.11170312986121075, 0.12411458873467864)]

torques =  {0: 'output',
           1: 172.08337910792062,
           2: 298.4155182973038,
           3: 0} 

cycles =  {0: 527308013125.3259} 

rigid_links =  []

GA=meshes_opt.MeshAssemblyOptimizer(connections = connections, 
                                  gear_speeds = gear_speeds,
                                  center_distances = center_distances_bounds,
                                  torques = torques,
                                  cycles = cycles,
                                  verbose = True)

#Optimization for gear set with center-distance closed to the minimum boundary
GA.Optimize(nb_sol=1, verbose=True)
print('Number of solutions:',len(GA.solutions))
solution=GA.solutions[-1]
#solution.SVGExport('name.txt',{6 : [0,0], 4 : [0.5,0]})
#solution.FreeCADExport('meshes_agb2')



