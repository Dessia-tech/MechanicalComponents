import mechanical_components.gears_assembly as gears
import numpy as npy
import networkx as nx

##########################
# script pour AGB Safran
##########################
erreur=0.02
data_cd=[[[148.6*1e-3,148.6*1e-3+0.01],[130.1*1e-3,130.1*1e-3+0.01],[137*1e-3,137*1e-3+0.01],
         [135*1e-3,135*1e-3+0.01],[142*1e-3,142*1e-3+0.01],[145.4*1e-3,145.4*1e-3+0.01]],[[117*1e-3,117*1e-3+0.01]]]
data_gear_set=[[(2,4),(4,6),(6,7),(7,0),(0,3),(3,5)],[(5,1)]]
data_speed={2:[9000*npy.pi/30*(1-erreur),9000*npy.pi/30],4:[20000*npy.pi/30*(1-erreur),
               20000*npy.pi/30],6:[11000*npy.pi/30,11000*npy.pi/30*(1+erreur)],7:[17000*npy.pi/30,
               17000*npy.pi/30*(1+erreur)],0:[1000*npy.pi/30,30000*npy.pi/30],
               3:[15000*npy.pi/30,15000*npy.pi/30*(1+erreur)],5:[1000*npy.pi/30,30000*npy.pi/30],
               1:[4100*npy.pi/30*(1-erreur),4100*npy.pi/30]}
data_rack={0:{'name':'Catalogue_A','module':[2.54*1e-3,2.54*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
data_rack_choice={2:0,4:0,6:0,7:0,0:0,3:0,5:0}
data_helix_angle=[[0,0],[0,0]]
data_material={2:gears.hardened_alloy_steel,4:gears.hardened_alloy_steel}
data_torque={2:106,7:-85,3:'output',1:186}
data_cycle={2:1e12}

list_gear=[]
gear_graph=[]
for list_gear_set in data_gear_set:
    lg=[]
    for le in list_gear_set:
        lg=list(set(list(set(le))+lg))
    list_gear.append(lg)
    gear_graph_elem=nx.Graph()
    gear_graph_elem.add_nodes_from(lg)
    gear_graph_elem.add_edges_from(list_gear_set)
    gear_graph.append(gear_graph_elem)
nb_plan_eng=len(list_gear)


    


#construction du premier plan d'engrenement
erreur=0.02
list_cd=[[148.6*1e-3,148.6*1e-3+0.01],[130.1*1e-3,130.1*1e-3+0.01],[137*1e-3,137*1e-3+0.01],
         [135*1e-3,135*1e-3+0.01],[142*1e-3,142*1e-3+0.01],[145.4*1e-3,145.4*1e-3+0.01]]
list_gear_set=[(2,4),(4,6),(6,7),(7,0),(0,3),(3,5)]
list_speed={2:[9000*npy.pi/30*(1-erreur),9000*npy.pi/30],4:[20000*npy.pi/30*(1-erreur),
               20000*npy.pi/30],6:[11000*npy.pi/30,11000*npy.pi/30*(1+erreur)],7:[17000*npy.pi/30,
               17000*npy.pi/30*(1+erreur)],0:[1000*npy.pi/30,30000*npy.pi/30],
               3:[15000*npy.pi/30,15000*npy.pi/30*(1+erreur)],5:[1000*npy.pi/30,30000*npy.pi/30]}
list_rack={0:{'name':'Catalogue_A','module':[2.54*1e-3,2.54*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
list_rack_choice={2:0,4:0,6:0,7:0,0:0,3:0,5:0}
list_helix_angle={2:[0,0]}
list_material={2:gears.hardened_alloy_steel,4:gears.hardened_alloy_steel}
list_torque={2:106,7:-85,3:300}
list_cycle={2:1e12}
#list_Z={2:[40,100],4:[40,100],6:[40,100],7:[40,100],0:[40,100],3:[40,100],5:[40,100]}

GA1=gears.GearAssemblyOptimizer(gear_set=list_gear_set,gear_speed=list_speed,
                                            center_distance=list_cd,rack_list=list_rack,
                                            rack_choice=list_rack_choice,helix_angle=list_helix_angle,
                                            material=list_material,torque=list_torque,cycle=list_cycle)
GA1.SearchCenterLine(nb_sol=1)
eng1=GA1.solutions_search[-1]
#sol.SVGExport('name.txt',{2:[0,0]})

#construction du second plan d'engrènement
erreur=0.05
list_cd2=[[117*1e-3,117*1e-3+0.01]]
list_gear_set2=[(5,1)]
list_speed2={5:[1000*npy.pi/30,30000*npy.pi/30],1:[4100*npy.pi/30*(1-erreur),
               4100*npy.pi/30]}
list_rack2={0:{'name':'Catalogue_A','module':[2.54*1e-3,2.54*1e-3],
              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
list_rack_choice2={5:0,1:0}
list_helix_angle2={5:[0,0]}
list_material2={5:gears.hardened_alloy_steel}
list_torque2={1:186}
list_cycle2={1:1e12}
#list_Z={2:[40,100],4:[40,100],6:[40,100],7:[40,100],0:[40,100],3:[40,100],5:[40,100]}

GA2=gears.GearAssemblyOptimizer(Z={},gear_set=list_gear_set2,gear_speed=list_speed2,center_distance=list_cd2,rack_list=list_rack2,rack_choice=list_rack_choice2,helix_angle=list_helix_angle2,material=list_material2,torque=list_torque2,cycle=list_cycle2)
GA2.SearchCenterLine(nb_sol=1)
eng2=GA2.solutions_search[-1]

eng2.SVGExport('name.txt',{5:[0,0]})
#v=sol.VolumeModel(name='3D')
#sol.FreeCADExport('Gears')