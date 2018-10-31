import sys
#del sys.modules['mechanical_components.optimization']
import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy
from interval import interval

###############################################################
#input for a specified gear mesh (modification of the modulus)
list_cd=[[0.117,0.16]]
list_gear_set=[[(0,1)]]
list_speed={0:[1000*npy.pi/30,1500*npy.pi/30],1:[4100*npy.pi/30,
               4300*npy.pi/30]}
list_rack={0:{'name':'Catalogue_A','module':[1*1e-3,2.5*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
list_rack_choice={0:[0],1:[0]}
list_helix_angle={0:[0,0]}
list_material={0:meshes_opt.hardened_alloy_steel}
list_torque={1:186,0:'output'}
list_cycle={1:1e12}

print('#####################################')
print('############ Decision Tree ##########')
print('#####################################')

GA=meshes_opt.MeshAssemblyOptimizer(Z={},
                               connections=list_gear_set,
                               gear_speed=list_speed,
                               center_distance=list_cd,
                               rack_list=list_rack,
                               rack_choice=list_rack_choice,
                               helix_angle=list_helix_angle,
                               material=list_material,
                               torque=list_torque,
                               cycle=list_cycle
                               )

#Optimization for gear set with center-distance closed to the minimum boundary
GA.OptimizeCD(nb_sol=5,verbose=True)
print('Nombre de solutions convergés:',len(GA.solutions))
solution=GA.solutions[-1]
#solution.SVGExport('name.txt',{5:[0,0]})
#solution.FreeCADExport('Gears1')

##Recherche non triée des nb_sol architecture vérifiant le CDC (nb_sol=-1 pour analyser l'ensemble des solutions)
#GA.Optimize(nb_sol=-1)
#print('Nombre de solutions convergés:',len(GA.solutions))
#solution=GA.solutions[-1]
#solution.SVGExport('meshes1.txt',{5:[0,0]})
##solution.FreeCADExport('meshes1')

###############################################################
#input for a specified gear mesh (modification of several parameter of the rack)
list_rack={0:{'name':'Catalogue_A','module':[2*1e-3,2.5*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[0.9,1],'coeff_gear_dedendum':[1.2,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
list_rack_choice={0:[0],1:[0]}

GA=meshes_opt.MeshAssemblyOptimizer(Z={},
                               connections=list_gear_set,
                               gear_speed=list_speed,
                               center_distance=list_cd,
                               rack_list=list_rack,
                               rack_choice=list_rack_choice,
                               helix_angle=list_helix_angle,
                               material=list_material,
                               torque=list_torque,
                               cycle=list_cycle)

#Optimization for gear set with center-distance closed to the minimum boundary
GA.OptimizeCD(nb_sol=1,verbose=True)
print('Nombre de solutions convergés:',len(GA.solutions))
solution=GA.solutions[-1]
#solution.SVGExport('meshes1.txt',{5:[0,0]})
#solution.FreeCADExport('meshes1')


###############################################################
#input for a specified gear mesh (optimization of two rack)
list_rack={0:{'name':'Catalogue_A','module':[1.8*1e-3,2*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
           1:{'name':'Catalogue_A','module':[2*1e-3,2*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
           2:{'name':'Catalogue_A','module':[1.9*1e-3,1.9*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
           3:{'name':'Catalogue_A','module':[1.8*1e-3,1.8*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
list_rack_choice={0:[0],1:[1,2,3]}

GA=meshes_opt.MeshAssemblyOptimizer(Z={},
                               connections=list_gear_set,
                               gear_speed=list_speed,
                               center_distance=list_cd,
                               rack_list=list_rack,
                               rack_choice=list_rack_choice,
                               helix_angle=list_helix_angle,
                               material=list_material,
                               torque=list_torque,
                               cycle=list_cycle)

#Optimization for gear set with center-distance closed to the minimum boundary
GA.OptimizeCD(nb_sol=10,verbose=True)
print('Nombre de solutions convergés:',len(GA.solutions))
solution=GA.solutions[-1]
solution.SVGExport('meshes1.txt',{0:[0,0]})
#solution.FreeCADExport('meshes1')