import sys
del sys.modules['mechanical_components.optimization']
import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy

#3 gears meshes test
list_cd=[[0.08,0.12],[0.07,0.12],[0.05,0.07]]
connections=[[(5,1),(6,7)],[(1,2)],[(6,8)]]
list_speed={5:[1000*npy.pi/30,1200*npy.pi/30],1:[2000*npy.pi/30,
               2100*npy.pi/30],2:[1000*npy.pi/30,1200*npy.pi/30]}
list_tpa={0:[18/180.*npy.pi,22/180.*npy.pi],1:[18/180.*npy.pi,22/180.*npy.pi]}
list_rack={0:{'name':'Racks_A','module':[2*1e-3,2*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
list_rack_choice={5 : [0],1 : [0]}
list_helix_angle={5 : [0,0]}
list_material={5 : meshes_opt.hardened_alloy_steel}
list_torque={1 : 186,5 : 'output',2 : 20}
list_cycle={1 : 1e12}

GA = meshes_opt.MeshAssemblyOptimizer(Z={},
                                  connections = connections, 
                                  gear_speed = list_speed,
                                  center_distance = list_cd,
                                  rack_list = list_rack,
                                  rack_choice = list_rack_choice,
                                  helix_angle = list_helix_angle,
                                  material = list_material,torque = list_torque,
                                  cycle = list_cycle,verbose=True,
                                  transverse_pressure_angle=list_tpa)

#Optimization for gear set with center-distance closed to the minimum boundary
GA.SearchOptimumCD(nb_sol=1, verbose=True)
print('Number of solutions:',len(GA.solutions))
solution=GA.solutions[-1]
solution.SVGExport('meshes2.txt',{5 : [0,0], 2 : [0.15,0]})
#solution.FreeCADExport('meshes2')
 

