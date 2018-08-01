import mechanical_components.gears_assembly as gears
import numpy as npy

#cas test avec deux engrenages
list_cd=[[0.117,0.117]]
list_gear_set=[(1,0)]
list_speed={1:[1000*npy.pi/30,1500*npy.pi/30],0:[4100*npy.pi/30,
               4300*npy.pi/30]}

GA=gears.GearAssemblyOptimizer(gear_set=list_gear_set,gear_speed=list_speed,
                               center_distance=list_cd)

#Recherche triée des nb_sol architecture ayant un entraxe mini (nb_sol=-1 pour analyser l'ensemble des solutions)
GA.SearchOptimumCD(nb_sol=1,post_traitement=True)
print('Nombre de solutions convergés:',len(GA.solutions))
solution=GA.solutions[-1]
print(solution.Dict())
#solution.SVGExport('name.txt',{5:[0,0]})
solution.FreeCADExport('Gears1',centers={0:(0,0.117*npy.sin(0.1),0.117*npy.cos(0.1)),1:(0,0,0)})

##Recherche non triée des nb_sol architecture vérifiant le CDC (nb_sol=-1 pour analyser l'ensemble des solutions)
#GA.Optimize(nb_sol=-1)
#print('Nombre de solutions convergés:',len(GA.solutions))
#solution=GA.solutions[-1]
#solution.SVGExport('name.txt',{5:[0,0]})
##solution.FreeCADExport('Gears1')
