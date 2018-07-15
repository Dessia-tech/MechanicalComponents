import mechanical_components.gears_assembly as gears
import numpy as npy
import networkx as nx
from scipy.optimize import minimize,fsolve
import itertools
import volmdlr as vm

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
data_rack_choice={2:0,4:0,6:0,7:0,0:0,3:0,5:0,1:0}
data_helix_angle=[[0,0],[0,0]]
data_material={2:gears.hardened_alloy_steel,4:gears.hardened_alloy_steel}
data_torque={2:106,7:-85,3:'output',1:186}
data_cycle={1:1e12,5:1e12}

composants_externe={'FP':{'diam':0.25,'long':0.5},'PMG':{'diam':0.23,'long':0.5},'Manual':{'diam':0.1,'long':0.5},
                    'VFG':{'diam':0.28,'long':0.5},'PMA':{'diam':0.15,'long':0.5},'Starter':{'diam':0.25,'long':0.5},
                    'TG':{'diam':0.25,'long':0.5},'HP':{'diam':0.25,'long':0.5},'GDL':{'diam':0.2,'long':0.5}
                    }

position_composant={2:{1:'bearing',3:'gear',4:'bearing',5:'FP'},
                    4:{0:'PMG',1:'bearing',3:'gear',4:'bearing'},
                    6:{1:'bearing',3:'gear',4:'bearing',5:'Manual'},
                    7:{0:'VFG',1:'bearing',3:'gear',4:'bearing',5:'PMA'},
                    0:{1:'bearing',3:'gear',4:'bearing'},
                    3:{0:'Starter',1:'bearing',3:'gear',4:'bearing',5:'TG'},
                    5:{1:'bearing',2:'gears',3:'gear',4:'bearing'},
                    1:{0:'HP',1:'bearing',2:'gear',4:'bearing',5:'GDL'},
                    }

list_gear=[]
gear_graph=[]
gear_set_totale=[]
list_gear_totale=[]
for num_plan_engr,list_gear_set in enumerate(data_gear_set):
    lg=[]
    for le in list_gear_set:
        list_gear=list(set(list(set(le))+list_gear))
        lg=list(set(list(set(le))+lg))
        gear_set_totale.append(tuple(list(le)+[{'num_plan_engr':num_plan_engr}]))
    gear_graph_elem=nx.Graph()
    gear_graph_elem.add_nodes_from(lg)
    gear_graph_elem.add_edges_from(list_gear_set)
    gear_graph.append(gear_graph_elem)
    list_gear_totale.append(lg)
nb_plan_eng=len(list_gear)

gear_graph_totale=nx.Graph()
gear_graph_totale.add_nodes_from(list_gear)
gear_graph_totale.add_edges_from(gear_set_totale)
ggd=gear_graph_totale.degree(list_gear)

for ne,tq in data_torque.items():
    if tq=='output':
        ne_output=ne
list_node_init=[]
for ne,nb_connexion in ggd:
    if (ne!=ne_output) and (nb_connexion==1) and (ne in data_torque.keys()):
        list_node_init.append(ne)
        
calcul_plan_engr={}
for i,nd_init in enumerate(list_node_init):
    path_elem=nx.shortest_path(gear_graph_totale,source=nd_init,target=ne_output)
    node_init=nd_init
    calcul_plan_engr[i]=[]
    for nd_prop in path_elem:
        edge_elem=gear_graph_totale.edges(nd_prop)
        li=[]
        for nd1,nd2 in edge_elem:
            li.append(gear_graph_totale[nd1][nd2]['num_plan_engr'])
        if ((len(li)>1) and (len(set(li))>1)) or (nd_prop==ne_output):
            num_plan=gear_graph_totale[nd_prop_m][nd_prop]['num_plan_engr']
            calcul_plan_engr[i].append({'num_plan':num_plan,'node_init':node_init,'node_fin':nd_prop})
            node_init=nd_prop
        nd_prop_m=nd_prop

ordre_calcul=[]
liste_node_init=[]
liste_node_fin=[]
for path,path_engr in calcul_plan_engr.items():
    for plan_engr in path_engr:
        ajout=plan_engr['num_plan']
        node_init=plan_engr['node_init']
        node_fin=plan_engr['node_fin']
        if ajout not in ordre_calcul:
            valid=True
            for p,pe in calcul_plan_engr.items():
                li=[]
                for pr in pe:
                    li.append(pr['num_plan'])
                if ajout in li:
                    i=0
                    while li[i]!=ajout:
                        if li[i] not in ordre_calcul:
                            valid=False
                        i+=1
            if valid:
                ordre_calcul.append(ajout)
                liste_node_init.append(node_init)
                liste_node_fin.append(node_fin)
        
erreur=0.02
sol_eng={}
for i,(num_plan,node_output) in enumerate(zip(ordre_calcul,liste_node_fin)):
    list_cd=data_cd[num_plan]
    list_gear_set=data_gear_set[num_plan]
    list_gear_temp=list_gear_totale[num_plan]
    list_rack=data_rack
    list_speed={}
    list_rack_choice={}
    list_helix_angle={}
    list_torque={}
    list_cycle={}
    for node in list_gear_temp:
        if node in data_speed.keys():
            list_speed[node]=data_speed[node]
        list_rack_choice[node]=data_rack_choice[node]
        if node in data_torque.keys():
            list_torque[node]=data_torque[node]
        if node in data_cycle.keys():
            list_cycle[node]=data_cycle[node]
    list_torque[node_output]='output'
    if i==0:
        list_helix_angle[num_plan]=[0,0]
        
    GA1=gears.GearAssemblyOptimizer(gear_set=list_gear_set,gear_speed=list_speed,Z={},
                                            center_distance=list_cd,rack_list=list_rack,
                                            rack_choice=list_rack_choice,helix_angle=list_helix_angle,
                                            torque=list_torque,cycle=list_cycle,
                                            safety_factor=4)
    GA1.SearchCenterLine(nb_sol=1)
    sol_eng[num_plan]=GA1.solutions_search[-1]
    for node,tq in sol_eng[num_plan].torque2.items():
        data_torque[node]=-tq
        
nb_eng_total=0
list_gear_complete=[]
for i in list_gear_totale:
    nb_eng_total+=len(i)
    for j in i:
        if j not in list_gear_complete:
            list_gear_complete.append(j)
    
diam_max={}
for node,composants in position_composant.items():
    dia_max=0.1
    if 0 in composants.keys():
        dia_max=max(dia_max,composants_externe[composants[0]]['diam'])
    if 5 in composants.keys():
        dia_max=max(dia_max,composants_externe[composants[5]]['diam'])
    diam_max[node]=dia_max
print(diam_max)
    
ordre_rangement=nx.dfs_edges(gear_graph_totale,2)

Rint=2
Rext=2.6
def fun(x):
    obj=0
    ine=ineg(x)
    for itera in ine:
        obj+=itera**2
    return obj
def eg(x):
    ine=[]
    for i,ne in enumerate(list_gear_complete):
        if ne==ne_output:
            ine.append(x[2*int(i)]-0)
            ine.append(x[2*int(i)+1]+2.3)
    return ine
def ineg(x):
    ine=[]
    for num_plan in ordre_calcul:
        list_gear_set=data_gear_set[num_plan]
        list_gear=sol_eng[num_plan].list_gear
        list_cd=sol_eng[num_plan].center_distance
        for num,it in enumerate(list_gear_set):
            eng1=(list_gear_complete).index(it[0])
            eng2=(list_gear_complete).index(it[1])
            ine.append(((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5-0.99999*list_cd[num])
            ine.append(1.00001*list_cd[num]-((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5)
#    for it1,it2 in itertools.combinations(list_gear_complete,2):
#        eng1=(list_gear_complete).index(it1)
#        eng2=(list_gear_complete).index(it2)
#        ine.append(((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5-(diam_max[it1]/2+diam_max[it2]/2))
    for (it1,it2) in ordre_rangement:
        eng1=(list_gear_complete).index(it1)
        eng2=(list_gear_complete).index(it2)
        ine.append(x[2*eng2]-x[2*eng1])
    for eng,node in enumerate(list_gear_complete):
        ine.append(((x[2*eng])**2+(x[2*eng+1])**2)**0.5-Rint)
        ine.append(Rext-((x[2*eng])**2+(x[2*eng+1])**2)**0.5)
    return ine
cons = ({'type': 'eq','fun' : eg},{'type': 'ineq','fun' : ineg})
drap=1
while drap==1:
    x0=tuple(npy.random.random(2*nb_eng_total)*3)
    Bound=[[-1,1],[-2.6,-1.5]]*(nb_eng_total)
    res = minimize(fun,x0, method='SLSQP', bounds=Bound,constraints=cons)
    if (min(ineg(res.x))>0) and (max(eg(res.x))<1e-7):
        drap=0
x_opt=res.x

data_SVG={}
for i,ne in enumerate(list_gear_complete):
    data_SVG[ne]=[0,x_opt[2*i],x_opt[2*i+1]]

list_gear=sol_eng[0].list_gear
centers=[]
for ne in list_gear:
    centers.append(data_SVG[ne])
model,primitives1=sol_eng[0].VolumeModel(centers)

list_gear=sol_eng[1].list_gear
centers=[]
for ne in list_gear:
    t1=data_SVG[ne]
    t1[0]=0.05
    centers.append(t1)
model,primitives2=sol_eng[1].VolumeModel(centers)

primitives=primitives1
primitives.extend(primitives2)
model=vm.VolumeModel(primitives)
model.FreeCADExport('python' ,'Gears1', '/usr/lib/freecad/lib', ['fcstd'])
#sol_eng[0].FreeCADExport('Gears1',centers)

##construction du premier plan d'engrenement
#erreur=0.02
#list_cd=[[148.6*1e-3,148.6*1e-3+0.001],[130.1*1e-3,130.1*1e-3+0.001],[137*1e-3,137*1e-3+0.001],
#         [135*1e-3,135*1e-3+0.001],[142*1e-3,142*1e-3+0.001],[145.4*1e-3,145.4*1e-3+0.001]]
#list_gear_set=[(2,4),(4,6),(6,7),(7,0),(0,3),(3,5)]
#list_speed={2:[9000*npy.pi/30*(1-erreur),9000*npy.pi/30],4:[20000*npy.pi/30*(1-erreur),
#               20000*npy.pi/30],6:[11000*npy.pi/30,11000*npy.pi/30*(1+erreur)],7:[17000*npy.pi/30,
#               17000*npy.pi/30*(1+erreur)],0:[1000*npy.pi/30,30000*npy.pi/30],
#               3:[15000*npy.pi/30,15000*npy.pi/30*(1+erreur)],5:[1000*npy.pi/30,30000*npy.pi/30]}
#list_rack={0:{'name':'Catalogue_A','module':[2.54*1e-3,2.54*1e-3],
#              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
#              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
#              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
#list_rack_choice={2:0,4:0,6:0,7:0,0:0,3:0,5:0}
#list_helix_angle={2:[0,0]}
#list_material={2:gears.hardened_alloy_steel,4:gears.hardened_alloy_steel}
#list_torque={2:106,7:-85,3:300}
#list_cycle={2:1e12}
#list_Z={2:[80,90],4:[35,45],6:[60,70],7:[40,50],0:[60,70],3:[45,55],5:[60,70]}
#
#GA1=gears.GearAssemblyOptimizer(gear_set=list_gear_set,gear_speed=list_speed,Z={},
#                                            center_distance=list_cd,rack_list=list_rack,
#                                            rack_choice=list_rack_choice,helix_angle=list_helix_angle,
#                                            material=list_material,torque=list_torque,cycle=list_cycle,
#                                            safety_factor=4)
#GA1.SearchCenterLine(nb_sol=1)
#eng1=GA1.solutions_search[-1]
#eng1.FreeCADExport('Gears1')
##sol.SVGExport('name.txt',{2:[0,0]})
#
##construction du second plan d'engrènement
#erreur=0.05
#list_cd2=[[117*1e-3,117*1e-3+0.01]]
#list_gear_set2=[(5,1)]
#list_speed2={5:[1000*npy.pi/30,30000*npy.pi/30],1:[4100*npy.pi/30*(1-erreur),
#               4100*npy.pi/30]}
#list_rack2={0:{'name':'Catalogue_A','module':[2.54*1e-3,2.54*1e-3],
#              'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
#              'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
#              'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
#list_rack_choice2={5:0,1:0}
#list_helix_angle2={5:[0,0]}
#list_material2={5:gears.hardened_alloy_steel}
#list_torque2={1:186}
#list_cycle2={1:1e12}
##list_Z={2:[40,100],4:[40,100],6:[40,100],7:[40,100],0:[40,100],3:[40,100],5:[40,100]}
#
#GA2=gears.GearAssemblyOptimizer(Z={},gear_set=list_gear_set2,gear_speed=list_speed2,center_distance=list_cd2,rack_list=list_rack2,rack_choice=list_rack_choice2,helix_angle=list_helix_angle2,material=list_material2,torque=list_torque2,cycle=list_cycle2)
#GA2.SearchCenterLine(nb_sol=1)
#eng2=GA2.solutions_search[-1]
#eng2.FreeCADExport('Gears2')
##eng2.SVGExport('name.txt',{5:[0,0]})
##v=sol.VolumeModel(name='3D')
#sol.FreeCADExport('Gears')