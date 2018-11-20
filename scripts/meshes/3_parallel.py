import sys
#del sys.modules['mechanical_components.optimization']
import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy
from interval import interval

#3 gears meshes test


connections = [[(0, 1), (2, 3), (4, 5)]]
rigid_links = [(0, 2), (2, 4)]
speeds = {0: (375.39001972643246, 398.6100209466242),
          1: (183.38456636758707, 194.7279416068193),
          2: (375.39001972643246, 398.6100209466242),
          3: (366.66143803370795, 389.3415269842466),
          4: (375.39001972643246, 398.6100209466242),
          5: (687.3592382738575, 729.8763045588385)}

  
{0: [71, 115], 1: [26, 95], 2: [31, 111], 3: [66, 115], 4: [63, 115], 5: [33, 115]}


center_distances = [(0.12779958830381627, 0.16136311654522256)]
cycles = {0: 15606510690.100288}
torques = {0: 'output', 1: -31.690621345499412, 2: -30, 3: -61.410547329865366, 4: -30, 5: -16.38403480073937}

GA = meshes_opt.MeshAssemblyOptimizer(Z={},
                                  connections = connections, 
                                  strong_link = rigid_links,
                                  gear_speed = speeds,
                                  center_distance = center_distances,
                                  cycle = cycles,verbose=True,
                                  torque = torques)



#Optimization for gear set with center-distance closed to the minimum boundary
GA.Optimize(nb_sol=10, verbose=True)
print('Number of solutions:',len(GA.solutions))
solution=GA.solutions[-1]
solution.SVGExport('meshes2.txt',{1 : [0,0]})
#solution.FreeCADExport('meshes2')
 

