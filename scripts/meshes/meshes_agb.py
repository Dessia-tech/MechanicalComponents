import mechanical_components.meshes as meshes
import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy

"""
7 gears Test case with fixed modulus to 2.56

"""
k_error=0.02

center_distances=[(0.09713462117912072, 0.12713462117912072),
                  (0.09791849011138422, 0.12791849011138423),
                  (0.09139747947854315, 0.12139747947854315),
                  (0.13361504530822224, 0.16361504530822224),
                  (0.08500000415463765, 0.11500000415463765),
                  (0.11156775444767722, 0.14156775444767722)]

#speeds_list = [1968.9802959880892, 682.3978991447033, 1909.1202983099663,
#          750.2180371710731, 723.8745865768974, 820.4085146766263, 954.5599648123681]
#
#speeds = {}
## Formatting
#for ispeed, speed in enumerate(speeds_list):
#    speeds[ispeed] = ((1-k_error) * speed, (1+k_error) * speed)

speeds = {0: [1943.3278688438672, 2022.6473736946373],
          1: [674.6053970306391, 702.1403111951549],
          2: [1870.9422413238008, 1947.3072307655884],
          3: [744.1576229906221, 774.5314035208515],
          4: [709.3952039966548, 738.3501102822327],
          5: [809.7193567370781, 842.769126399816],
          6: [935.4684378940851, 973.6508231142518]}

connections = [(0, 1), (1, 2), (2, 3), (3, 4), (0, 5), (5, 6)]
#list_speed={2:[9000*npy.pi/30*(1-erreur),9000*npy.pi/30],4:[20000*npy.pi/30*(1-erreur),
#               20000*npy.pi/30],6:[11000*npy.pi/30,11000*npy.pi/30*(1+erreur)],7:[17000*npy.pi/30,
#               17000*npy.pi/30*(1+erreur)],0:[1000*npy.pi/30,30000*npy.pi/30],
#               3:[15000*npy.pi/30,15000*npy.pi/30*(1+erreur)],5:[10000*npy.pi/30,11000*npy.pi/30]}

list_rack = {0:{'name':'Catalogue_A','module':[2.43*1e-3,2.43*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}

rack_choices = {0:[0], 1:[0], 2:[0], 3:[0], 4:[0], 5:[0], 6:[0]}

torques = {0: -16.380372067375156, 1: 0, 2: 3.7052066976022022,
           3: 0, 4: 27.221926757893733, 5: 0, 6: 'output'}

GA=meshes_opt.MeshAssemblyOptimizer(connections = connections, 
                                  gear_speed = speeds,
                                  center_distance = center_distances,
                                  rack_list = list_rack,
                                  torque = torques,
                                  rack_choice=rack_choices)

#Recherche triée des nb_sol architecture ayant un entraxe mini (nb_sol=-1 pour analyser l'ensemble des solutions)
GA.Optimize(nb_sol=10, verbose = True)
print('Number of solutions:',len(GA.solutions))
solution=GA.solutions[0]
#solution.SVGExport('name.txt',{5:[0,0]})
solution.FreeCADExport('meshes_agb')


#Recherche non triée des nb_sol architecture vérifiant le CDC (nb_sol=-1 pour analyser l'ensemble des solutions)
#GA.Optimize(nb_sol=-1,post_traitement=True)
#print('Nombre de solutions convergés:',len(GA.solutions))
#solution=GA.solutions[-1]
#solution.SVGExport('name.txt',{5:[0,0]})
#solution.FreeCADExport('Gears1')

