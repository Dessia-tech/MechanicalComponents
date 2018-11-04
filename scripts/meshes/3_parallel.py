import sys
#del sys.modules['mechanical_components.optimization']
import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy
from interval import interval

#3 gears meshes test


connections = [[(0, 1), (2, 3), (4, 5)]]
rigid_links = [(0, 2), (2, 4)]
speeds = {0: (393.24146705049753, 394.8175851549084),
          1: (466.28776359111964, 468.1566524231482),
          2: (393.24146705049753, 394.8175851549084),
          3: (244.35334237247366, 245.33271448619098),
          4: (393.24146705049753, 394.8175851549084),
          5: (932.5242979862064, 936.2618703228245)}

center_distances = [(0.166883767540157, 0.1854264083779522)]
cycles = {0: 15918810165.932156}
torques = {0: 'output', 1: -25.300350840558927, 2: 'output',
           3: -48.27944605534433, 4: 'output', 5: -12.650870370875234}


GA = meshes_opt.MeshAssemblyOptimizer(Z={},
                                  connections = connections, 
                                  strong_link = rigid_links,
                                  gear_speed = speeds,
                                  center_distance = center_distances,
                                  cycle = cycles,verbose=True)



#Optimization for gear set with center-distance closed to the minimum boundary
GA.Optimize(nb_sol=10, verbose=True)
print('Number of solutions:',len(GA.solutions))
solution=GA.solutions[-1]
solution.SVGExport('meshes2.txt',{1 : [0,0]})
#solution.FreeCADExport('meshes2')
 

