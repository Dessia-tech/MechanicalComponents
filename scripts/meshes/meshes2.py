import mechanical_components.meshes as meshes
import numpy as npy


"""
3 gear test case with 0.5 to 2 for modulus

"""
list_cd=[[0.05,0.055],[0.05,0.055]]
connections=[(5,1),(1,2)]
list_speed={5:[1000*npy.pi/30,1100*npy.pi/30],1:[2000*npy.pi/30,
               2100*npy.pi/30],2:[1000*npy.pi/30,1100*npy.pi/30]}
list_rack={0:{'name':'Racks_A','module':[0.5*1e-3,2.3*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
list_rack_choice={5:[0],1:[0]}
list_helix_angle={5:[0,0]}
list_material={5:meshes.hardened_alloy_steel}
list_torque={1:186,5:'output',2:20}
list_cycle={1:1e12}

GA = meshes.MeshAssemblyOptimizer(Z={},
                                  connections = connections, 
                                  gear_speed = list_speed,
                                  center_distance = list_cd,
                                  rack_list = list_rack,
                                  rack_choice = list_rack_choice,
                                  helix_angle = list_helix_angle,
                                  material = list_material,torque = list_torque,
                                  cycle = list_cycle)

#Recherche triée des nb_sol architecture ayant un entraxe mini (nb_sol=-1 pour analyser l'ensemble des solutions)
GA.SearchOptimumCD(nb_sol=10, verbose=True)
print('Number of solutions:',len(GA.solutions))
solution=GA.solutions[-1]
solution.FreeCADExport('Gears2')
 

