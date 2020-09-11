import sys
#del sys.modules['mechanical_components.optimization']
import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy
# from interval import interval

#3 gears meshes test


connections = [[(0, 1)]]
rigid_links = []
speeds = {0: (248.4937141433299, 263.864459348072), 1: (469.23940372838325, 498.2645214847781)}

center_distances = [(0.19381369521535707, 0.2349256911701298)]
cycles =  {0: 17383903277.550545}
torques = {0: 'output', 1: -10}

   

GA = meshes_opt.MeshAssemblyOptimizer(Z = {},
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
 

