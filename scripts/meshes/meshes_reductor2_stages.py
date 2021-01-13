#import sys
#del sys.modules['mechanical_components.optimization']
import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy

"""
7 gears Test case with fixed modulus to 2.56

"""
k_error=0.02

center_distances=[(0.09713462117912072, 0.12713462117912072)]

#speeds_list = [1968.9802959880892, 682.3978991447033, 1909.1202983099663,
#          750.2180371710731, 723.8745865768974, 820.4085146766263, 954.5599648123681]
#
#speeds = {}
## Formatting
#for ispeed, speed in enumerate(speeds_list):
#    speeds[ispeed] = ((1-k_error) * speed, (1+k_error) * speed)

speeds = {0: [243.620149773094, 253.56382935566927], 1: [92.0103009101568, 95.76582339628565]}

connections = [[(0, 1)]]
#list_speed={2:[9000*npy.pi/30*(1-erreur),9000*npy.pi/30],4:[20000*npy.pi/30*(1-erreur),
#               20000*npy.pi/30],6:[11000*npy.pi/30,11000*npy.pi/30*(1+erreur)],7:[17000*npy.pi/30,
#               17000*npy.pi/30*(1+erreur)],0:[1000*npy.pi/30,30000*npy.pi/30],
#               3:[15000*npy.pi/30,15000*npy.pi/30*(1+erreur)],5:[10000*npy.pi/30,11000*npy.pi/30]}

#list_rack = {0:{'name':'Catalogue_A','module':[2.43*1e-3,2.43*1e-3],
#              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
#              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
#              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}

#rack_choices = {0:[0], 1:[0], 2:[0], 3:[0], 4:[0], 5:[0], 6:[0]}

torques = {0: -16.380372067375156, 1: 'output'}

GA=meshes_opt.MeshAssemblyOptimizer(connections = connections, 
                                  gear_speeds = speeds,
                                  center_distances = center_distances,
#                                  rack_list = list_rack,
                                  torques = torques,
#                                  rack_choice=rack_choices)
                                  cycles={0:1e8}
                                  )
# GA.FreeCADExport('meshes')

#Recherche triée des nb_sol architecture ayant un entraxe mini (nb_sol=-1 pour analyser l'ensemble des solutions)
GA.Optimize(nb_sol=3, verbose = True)
print('Number of solutions:',len(GA.solutions))
#for solution in GA.solutions:
#    for a in solution.meshes.values():
#        Z1 = a[0].Z
#        Z2 = a[1].Z
#        print(243*Z1/Z2, 253*Z1/Z2)
#solution=GA.solutions[0]
#solution.SVGExport('name.txt',{5:[0,0]})
#solution.FreeCADExport('meshes_agb')


#Recherche non triée des nb_sol architecture vérifiant le CDC (nb_sol=-1 pour analyser l'ensemble des solutions)
#GA.Optimize(nb_sol=-1,post_traitement=True)
#print('Nombre de solutions convergés:',len(GA.solutions))
#solution=GA.solutions[-1]
#solution.SVGExport('name.txt',{5:[0,0]})
#solution.FreeCADExport('Gears1')

