import mechanical_components.meshes as meshes
import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy

"""
7 gears Test case with fixed modulus to 2.56

"""
erreur=0.02
list_cd=[[148.6*1e-3,148.6*1e-3+0.01],[130.1*1e-3,130.1*1e-3+0.01],[137*1e-3,137*1e-3+0.01],
         [135*1e-3,135*1e-3+0.01],[142*1e-3,142*1e-3+0.01],[145.4*1e-3,145.4*1e-3+0.01]]

list_gear_set=[(2,4),(4,6),(6,7),(7,0),(0,3),(3,5)]

list_speed={2:[9000*npy.pi/30*(1-erreur),9000*npy.pi/30],4:[20000*npy.pi/30*(1-erreur),
               20000*npy.pi/30],6:[11000*npy.pi/30,11000*npy.pi/30*(1+erreur)],7:[17000*npy.pi/30,
               17000*npy.pi/30*(1+erreur)],0:[1000*npy.pi/30,30000*npy.pi/30],
               3:[15000*npy.pi/30,15000*npy.pi/30*(1+erreur)],5:[10000*npy.pi/30,11000*npy.pi/30]}

list_rack={0:{'name':'Catalogue_A','module':[2.56*1e-3,2.56*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}

list_choice={0:[0],2:[0],3:[0],4:[0],5:[0],6:[0],7:[0]}

list_torque={2:106,5:-85,3:'output'}

GA=meshes_opt.MeshAssemblyOptimizer(connections = list_gear_set, 
                                  gear_speed = list_speed,
                                  center_distance = list_cd,
                                  rack_list = list_rack,
                                  torque = list_torque,
                               rack_choice=list_choice)

#Recherche triée des nb_sol architecture ayant un entraxe mini (nb_sol=-1 pour analyser l'ensemble des solutions)
GA.SearchOptimumCD(nb_sol=10, verbose = True)
print('Number of solutions:',len(GA.solutions))
solution=GA.solutions[0]
#solution.SVGExport('name.txt',{5:[0,0]})
solution.FreeCADExport('meshes3')


#Recherche non triée des nb_sol architecture vérifiant le CDC (nb_sol=-1 pour analyser l'ensemble des solutions)
#GA.Optimize(nb_sol=-1,post_traitement=True)
#print('Nombre de solutions convergés:',len(GA.solutions))
#solution=GA.solutions[-1]
#solution.SVGExport('name.txt',{5:[0,0]})
#solution.FreeCADExport('Gears1')

