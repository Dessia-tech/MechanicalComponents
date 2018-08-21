import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy

#Optimization of one gear mesh with a fixed center-distance
list_cd=[[0.117,0.117]]
list_gear_set=[(1,0)]
list_speed={1:[1000*npy.pi/30,1500*npy.pi/30],0:[4100*npy.pi/30,
               4300*npy.pi/30]}

GA = meshes_opt.MeshAssemblyOptimizer(connections = list_gear_set,
                                gear_speed = list_speed,
                                center_distance = list_cd)

#Optimization for a short list of architecture generate with the decision tree
for plex in GA.plex_calcul:
    print(plex['Z'])
GA.Optimize(list_sol=[1,2,3,4], verbose=True)

#Optimization for gear set with center-distance closed to the minimum boundary
GA.SearchOptimumCD(nb_sol=1, verbose=True)

#Export SVG and FreeCAD
print('Nombre de solutions convergés:',len(GA.solutions))
solution=GA.solutions[-1]
solution.SVGExport('name.txt',{0 : [0,0]})
solution.FreeCADExport('Gears1',centers = {0 : (0,0.117*npy.sin(0.1),0.117*npy.cos(0.1)),1 : (0,0,0)})

