import sys
#del sys.modules['mechanical_components.optimization']
import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy
from interval import interval

#3 gears meshes test


connections = [[(0, 1), (2, 3), (4, 5)]]
rigid_links = [(0, 2), (2, 4)]
speeds = {0: (183.38381539454952, 194.72714418184125),
          1: (469.23779247859505, 498.26281057005457),
          2: (375.3900158272758, 398.61001680628254),
          3: (183.38381539454952, 194.72714418184125),
          4: (375.3900158272758, 398.61001680628254),
          5: (687.3581881253122, 729.8751894526511)}


  
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
 

