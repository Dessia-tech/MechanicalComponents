import sys
del sys.modules['mechanical_components.optimization']
import mechanical_components.optimization.meshes as meshes_opt
import numpy as npy

#Optimization of one gear mesh with a fixed center-distance
list_cd=[[0.117,0.117]]
list_gear_set=[[(1,0)]]
list_speed={1:[1000*npy.pi/30,1000*npy.pi/30],0:[4100*npy.pi/30,
               4300*npy.pi/30]}

print('#####################################')
print('############ Decision Tree ##########')
print('#####################################')
GA = meshes_opt.MeshAssemblyOptimizer(Z={}, connections = list_gear_set,
                                gear_speed = list_speed,
                                center_distance = list_cd, verbose = True)

#Optimization for a short list of architecture generate with the decision tree
print('#####################################')
print('############ Simple Optimizer #######')
print('#####################################')
for plex in GA.plex_calcul:
    print(plex['Z'])
GA.Optimize(list_sol=[1,2,3,4], verbose=True)



print('#####################################')
print('############ Analyse Test ###########')
print('#####################################')
def pgcd(a,b) :
    while a%b != 0 :
        a, b = b, a%b
    return b
compt_check_false=0
for plex in GA.plex_calcul:
    Z=plex['Z']
    # validation of gear speed
    w={}    
    w[1]=1000*npy.pi/30
    w[0]=Z[1]/Z[0]*w[1]
    check_valid=True
    for num_mesh,interval_speed in list_speed.items():
        if round(w[num_mesh],4) not in interval([round(interval_speed[0],4),round(interval_speed[1],4)]):
            check_valid=False
    # validation of pgcd
    if pgcd(Z[1],Z[0])!=1:
        check_valid=False
    # validation of internal ratio
    if Z[1]/Z[0] not in interval([1/9.,9]):
        check_valid=False
    if check_valid==False:
        compt_check_false+=1
#        print('plex with teetch number: {}'.format(Z))
#        for num_mesh,speed_dat in w.items():
#            print('speed gear {}:{}'.format(num_mesh,speed_dat))
#        print(list_speed)
print('Number of False solution is:{}'.format(compt_check_false))

print('#####################################')
print('######## Intelligent optimizer ######')
print('#####################################')
#Optimization for gear set with center-distance closed to the minimum boundary
GA.SearchOptimumCD(nb_sol=1, verbose=True)

#Export SVG and FreeCAD
print('Nombre de solutions convergés:',len(GA.solutions))
solution=GA.solutions[-1]
solution.SVGExport('name.txt',{0 : [0,0]})
solution.FreeCADExport('Gears1',centers = {0 : (0,0.117*npy.sin(0.1),0.117*npy.cos(0.1)),1 : (0,0,0)})

