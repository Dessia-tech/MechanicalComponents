import mechanical_components.optimization.meshes as meshes_opt
import mechanical_components.meshes as me
import numpy as npy

# 7 gears Test case with fixed modulus to 2
# definition of input data



rack=meshes_opt.RackOpti(module=[1.5*1e-3,2.5*1e-3],transverse_pressure_angle_0=[20/180.*npy.pi,20/180.*npy.pi],
             coeff_gear_addendum=[1,1],coeff_gear_dedendum=[1.25,1.25],coeff_root_radius=[0.38,0.38],
             coeff_circular_tooth_thickness=[0.5,0.5],helix_angle=[21,60],)


meshopti1 = meshes_opt.MeshOpti(rack=rack, torque_input= 'output', speed_input=(1878.1453579221634, 1974.460504482274))
meshopti2 = meshes_opt.MeshOpti(rack=rack, torque_input= 40, speed_input=(449.8807309958231, 472.95153771355757))
meshopti3 = meshes_opt.MeshOpti(rack=rack, torque_input= -10, speed_input=(1875.9428646237072, 1972.145062809538))
meshopti4 = meshes_opt.MeshOpti(rack=rack, torque_input= 0, speed_input=(480.3124107238717, 504.94381640201897))
meshopti5 = meshes_opt.MeshOpti(rack=rack, torque_input= 0, speed_input=(700.6720400543617, 736.6039395443289))
meshopti6 = meshes_opt.MeshOpti(rack=rack, torque_input= 0, speed_input=(725.5314102562529, 762.738149243753))


center_distance1 = meshes_opt.CenterDistanceOpti((0.11134984458664793, 0.1293457790652981),[meshopti1,meshopti2])

center_distance2 = meshes_opt.CenterDistanceOpti((0.11137510812362536, 0.12937512559815068),[meshopti2,meshopti3])
center_distance3 = meshes_opt.CenterDistanceOpti((0.10568352491185784, 0.1227636905541783),[meshopti3,meshopti4])
center_distance4 = meshes_opt.CenterDistanceOpti((0.14181884668124436, 0.164739064326698),[meshopti4,meshopti5])
center_distance5 = meshes_opt.CenterDistanceOpti((0.07722003844178454, 0.0897000446545982),[meshopti1,meshopti6])






center_distances = [center_distance1,center_distance2,center_distance3,
                  center_distance4,center_distance5,]




cycles = {0: 1272321481513.054,2: 1272321481513.054}




torques = {0: 'output', 1: 40, 2: 10.0, 3: 0, 4: 0, 5: 0, 6: 0}






GA=meshes_opt.MeshAssemblyOptimizer(center_distances,cycles)

#Optimization for gear set with center-distance closed to the minimum boundary
GA.Optimize(nb_sol = 5, verbose=True)
print('Number of solutions:',len(GA.solutions))
solution=GA.solutions[0]
# solution.pos_axis({0:(0,0,0)})


    
for mesh_assembly in GA.solutions:
    m=mesh_assembly.mesh_combinations
    

#solution=GA.solutions[-1]
#solution.SVGExport('name.txt',{6 : [0,0], 4 : [0.5,0]})
#solution.FreeCADExport('meshes3')


