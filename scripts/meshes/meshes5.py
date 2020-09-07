import sys
#del sys.modules['mechanical_components.optimization']
import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy
# from interval import interval

#3 gears meshes test
list_cd=[[0.1,0.15]]
connections=[[(1,2),(3,4),(5,6)]]
list_strong_link=[[2, 4], [4, 6]]
list_speed={1:[900*npy.pi/30,900*npy.pi/30],
            2:[2000*npy.pi/30,2100*npy.pi/30],
            3:[900*npy.pi/30,920*npy.pi/30],
            5:[900*npy.pi/30,910*npy.pi/30],
            6:[2000*npy.pi/30,2100*npy.pi/30]}
list_tpa={0:[15/180.*npy.pi,25/180.*npy.pi],1:[15/180.*npy.pi,25/180.*npy.pi]}
list_rack={0:{'name':'Racks_A','module':[1.8*1e-3,2*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
           1:{'name':'Racks_A','module':[1.8*1e-3,2*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
list_rack_choice={1 : [0], 2 : [0], 3 : [1], 4: [1]}
list_helix_angle={1 : [0,0]}
list_material={1 : meshes_opt.hardened_alloy_steel}
list_torque={3 : 186,1 : 'output'}
list_cycle={1 : 1e12}

GA = meshes_opt.MeshAssemblyOptimizer(Z={},
                                  connections = connections, 
                                  strong_link = list_strong_link,
                                  gear_speed = list_speed,
                                  center_distance = list_cd,
                                  rack_list = list_rack,
                                  rack_choice = list_rack_choice,
                                  helix_angle = list_helix_angle,
                                  material = list_material,torque = list_torque,
                                  cycle = list_cycle,verbose=True,
                                  transverse_pressure_angle=list_tpa)

# def pgcd(a,b) :
#     while a%b != 0 :
#         a, b = b, a%b
#     return b
# compt_check_false=0
# for plex in GA.plex_calcul:
#     Z=plex['Z']
#     # validation of gear speed
#     w={}    
#     w[1]=900*npy.pi/30
#     w[2]=Z[1]/Z[2]*w[1]
#     w[4]=w[2]
#     w[3]=Z[4]/Z[3]*w[4]
#     w[6]=w[2]
#     w[5]=Z[6]/Z[5]*w[6]
#     check_valid=True
#     for num_mesh,interval_speed in list_speed.items():
#         if round(w[num_mesh],4) not in interval([round(interval_speed[0],4),round(interval_speed[1],4)]):
#             check_valid=False
#     # validation of pgcd
#     if pgcd(Z[1],Z[2])!=1:
#         check_valid=False
#     if pgcd(Z[2],Z[4])!=1:
#         check_valid=False
#     if pgcd(Z[3],Z[4])!=1:
#         check_valid=False
#     if pgcd(Z[4],Z[6])!=1:
#         check_valid=False
#     if pgcd(Z[5],Z[6])!=1:
#         check_valid=False
#     # validation of internal ratio
#     if Z[1]/Z[2] not in interval([1/9.,9]):
#         check_valid=False
#     if Z[3]/Z[4] not in interval([1/9.,9]):
#         check_valid=False
#     if Z[5]/Z[6] not in interval([1/9.,9]):
#         check_valid=False
#     if check_valid==False:
#         compt_check_false+=1
# #        print('plex with teetch number: {}'.format(Z))
# #        for num_mesh,speed_dat in w.items():
# #            print('speed gear {}:{}'.format(num_mesh,speed_dat))
# #        print(list_speed)
# print('Number of False solution is:{}'.format(compt_check_false))

#Optimization for gear set with center-distance closed to the minimum boundary
GA.OptimizeCD(nb_sol=10, verbose=True)
print('Number of solutions:',len(GA.solutions))
solution=GA.solutions[-1]
solution.SVGExport('meshes2.txt',{1 : [0,0]})
#solution.FreeCADExport('meshes2')
 

