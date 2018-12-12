import sys
#del sys.modules['mechanical_components.optimization']
import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy

# 7 gears Test case with fixed modulus to 2
# definition of input data
center_distances=[(0.09, 0.14),
                  (0.09, 0.137),
                  (0.09, 0.131),
                  (0.13, 0.163),
                  (0.085, 0.115),
                  (0.12, 0.14)]
speeds = {0: [1950, 1950],
          1: [670, 700],
          2: [1900, 1950],
          3: [750, 780],
          4: [700, 740],
          5: [800, 850],
          6: [950, 980]}
connections = [[(0, 1)], [(1, 2)], [(2, 3)], [(3, 4)], [(0, 5)], [(5, 6)]]
list_rack = {0:{'name':'Catalogue_A','module':[2*1e-3,2*1e-3],
              'transverse_pressure_angle_rack':[20/180.*npy.pi,20/180.*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
rack_choices = {0:[0], 1:[0], 2:[0],3:[0], 4:[0], 5:[0], 6:[0]}
torques = {0: 16.4, 4: 30, 5: 0, 6: 'output'}

GA=meshes_opt.MeshAssemblyOptimizer(connections = connections, 
                                  gear_speeds = speeds,
                                  center_distances = center_distances,
                                  rack_list = list_rack,
                                  torques = torques,
                                  rack_choice=rack_choices,
                                  verbose = True,
                                  cycles={0:1e8})

def pgcd(a,b) :
    while a%b != 0 :
        a, b = b, a%b
    return b
compt_check_false=0
for plex in GA.plex_calcul:
    Z=plex['Z']
    # validation of gear speed
    w={}    
    w[0]=1950.
    w[1]=Z[0]/Z[1]*w[0]
    w[2]=Z[1]/Z[2]*w[1]
    w[3]=Z[2]/Z[3]*w[2]
    w[4]=Z[3]/Z[4]*w[3]
    w[5]=Z[0]/Z[5]*w[0]
    w[6]=Z[5]/Z[6]*w[5]
    check_valid=True
    for num_mesh,interval_speed in speeds.items():
        if round(w[num_mesh],4) not in interval([round(interval_speed[0],4),round(interval_speed[1],4)]):
            check_valid=False
    # validation of pgcd
    if pgcd(Z[0],Z[1])!=1:
        check_valid=False
    if pgcd(Z[1],Z[2])!=1:
        check_valid=False
    if pgcd(Z[2],Z[3])!=1:
        check_valid=False
    if pgcd(Z[3],Z[4])!=1:
        check_valid=False
    if pgcd(Z[0],Z[5])!=1:
        check_valid=False
    if pgcd(Z[5],Z[6])!=1:
        check_valid=False
    # validation of internal ratio
    if Z[0]/Z[1] not in interval([1/9.,9]):
        check_valid=False
    if Z[1]/Z[2] not in interval([1/9.,9]):
        check_valid=False
    if Z[2]/Z[3] not in interval([1/9.,9]):
        check_valid=False
    if Z[3]/Z[4] not in interval([1/9.,9]):
        check_valid=False
    if Z[0]/Z[5] not in interval([1/9.,9]):
        check_valid=False
    if Z[5]/Z[6] not in interval([1/9.,9]):
        check_valid=False
    if check_valid==False:
        compt_check_false+=1
        print('plex with teetch number: {}'.format(Z))
        for num_mesh,speed_dat in w.items():
            print('speed gear {}:{}'.format(num_mesh,speed_dat))
        print(list_speed)
print('Number of False solution is:{}'.format(compt_check_false))

#Optimization for gear set with center-distance closed to the minimum boundary
GA.OptimizeCD(nb_sol = 1, verbose = True)
print('Number of solutions:', len(GA.solutions))
solution = GA.solutions[-1]
#solution.SVGExport('meshes3.txt', {6 : [0,0], 4 : [0.5,0]})
#solution.FreeCADExport('meshes3')

solution.FreeCADExport(fcstd_filepath = 'meshes3', python_path = '/Applications/FreeCAD.app/Contents/MacOS/FreeCADCmd',
            path_lib_freecad = '/Applications/FreeCAD.app/Contents/lib')


