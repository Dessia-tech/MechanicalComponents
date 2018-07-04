import mechanical_components.gears_assembly as gears
import numpy as npy

input_dict=(
          {'type':'gear_speed','data':{'0':[1000,1200],'2':[4000,4500],'5':[3000,3200],'3':[2000,2200]}},
          {'type':'gear_set','data':[(0,1),(1,2),(3,5),(3,4),(1,3)]},
          {'type':'center_distance','data':[[70,80],[60,100],[70,100],[60,100],[75,100]]},
          {'type':'Z','data':{'0':[10,90],'1':[10,90]}},
          {'type':'transverse_pressure_angle','data':[[15,22]]},
          {'type':'helix_angle','data':{'0':[20,30]}},
          {'type':'gear_width','data':{'3':[40,50],'2':[10,20]}},
          {'type':'frequency','data':[[1000,1200],[500,550],[3000,3500]]},
          {'type':'coefficient_profile_shift','data':{'0':[-1,1],'1':[-1,1],'2':[-1,1]}},
          )

input_dict=(
          {'type':'gear_speed','data':{'0':[800,820],'1':[970,1000],'3':[1050,1080]}},
          {'type':'gear_set','data':[(0,1),(1,2),(1,3)]},
          {'type':'center_distance','data':[[80,90],[80,90],[77,80]]},
          {'type':'Z','data':{'0':[10,80],'1':[10,80],'2':[10,80],'3':[10,80]}},
          {'type':'transverse_pressure_angle','data':[[10,21],[10,21],[10,21]]},
          {'type':'helix_angle','data':{'0':[20,30]}},
          {'type':'gear_width','data':{'0':[40,50],'1':[10,20]}},
          {'type':'frequency','data':[[500,550],[700,1200]]},
          {'type':'coefficient_profile_shift','data':{'0':[-1,1],'1':[-1,1]}},
#          {'type':'rack_list','data':{'0':{'module':1,'transverse_pressure_angle':20},
#                                      '1':{'module':1.25,'transverse_pressure_angle':20},
#                                      '2':{'module':1.5,'transverse_pressure_angle':20},
#                                      '3':{'module':2,'transverse_pressure_angle':20},
#                                      '4':{'module':2.5,'transverse_pressure_angle':20},
#                                      'x':None,'y':None
#                                      }},
#          {'type':'rack','data':{'0':['0','1','2','3'],'1':['0','1','2','3']}}
          )

input_dict=(
          {'type':'gear_speed','data':{'0':[15500,15900],'1':[15000,15500],'2':[16000,16500]}},
          {'type':'gear_set','data':[(0,1),(1,2)]},
          {'type':'center_distance','data':[[130,130],[140,145]]},
          {'type':'Z','data':{'0':[40,150],'1':[40,150],'2':[40,150]}},
          {'type':'transverse_pressure_angle','data':[[15,30],[15,30]]},
          {'type':'helix_angle','data':{'0':[20,30]}},
          {'type':'gear_width','data':{'0':[40,50],'1':[10,20]}},
          {'type':'frequency','data':[[500,550],[700,1200]]},
          {'type':'coefficient_profile_shift','data':{'0':[-1,1],'1':[-1,1]}},
          {'type':'rack_list','data':{'0':{'type':'catalog','name':'Catalogue_A','module':[0.8,0.8],'transverse_pressure_angle_rack':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
                                      '1':{'type':'catalog','name':'Catalogue_A','module':[0.9,0.9],'transverse_pressure_angle_rack':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
                                      '2':{'type':'catalog','name':'Catalogue_A','module':[1,1],'transverse_pressure_angle_rack':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
                                      '3':{'type':'catalog','name':'Catalogue_A','module':[1.25,1.25],'transverse_pressure_angle_rack':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
                                      '4':{'type':'catalog','name':'Catalogue_A','module':[1.5,1.5],'transverse_pressure_angle_rack':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
                                      '5':{'type':'catalog','name':'Catalogue_A','module':[1.75,1.75],'transverse_pressure_angle_rack':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
                                      '6':{'type':'catalog','name':'Catalogue_A','module':[2,2],'transverse_pressure_angle_rack':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
                                      '7':{'type':'catalog','name':'Catalogue_A','module':[2.5,2.5],'transverse_pressure_angle_rack':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
                                      '10':{'type':'optim','name':'Rack_A','module':[1,2.6],'transverse_pressure_angle_rack':[19,21],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
                                      '11':{'type':'optim','name':'Rack_B','module':[1,2.6],'transverse_pressure_angle_rack':[19,21],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
                                      '12':{'type':'optim','name':'Rack_C','module':[1,2.6],'transverse_pressure_angle_rack':[19,21],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
                                      }},
          {'type':'rack_choice','data':{'0':['10'],'1':['11'],'2':['12']}}
          )

#input_dict=(
#          {'type':'gear_speed','data':{'0':[800,820],'1':[980,1000],'2':[1150,1160]}},
#          {'type':'gear_set','data':[(0,1),(1,2)]},
#          {'type':'center_distance','data':[[80,100],[70,100]]},
#          {'type':'Z','data':{'0':[10,80],'1':[10,80],'2':[10,80]}},
#          {'type':'transverse_pressure_angle','data':[[15,30],[15,30]]},
#          {'type':'helix_angle','data':{'0':[20,30]}},
#          {'type':'gear_width','data':{'0':[40,50],'1':[10,20]}},
#          {'type':'frequency','data':[[500,550],[700,1200]]},
#          {'type':'coefficient_profile_shift','data':{'0':[-1,1],'1':[-1,1]}},
#          {'type':'rack_list','data':{'0':{'type':'catalog','name':'Catalogue_A','module':[0.8,0.8],'transverse_pressure_angle':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
#                                      '1':{'type':'catalog','name':'Catalogue_A','module':[0.9,0.9],'transverse_pressure_angle':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
#                                      '2':{'type':'catalog','name':'Catalogue_A','module':[1,1],'transverse_pressure_angle':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
#                                      '3':{'type':'catalog','name':'Catalogue_A','module':[1.25,1.25],'transverse_pressure_angle':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
#                                      '4':{'type':'catalog','name':'Catalogue_A','module':[1.5,1.5],'transverse_pressure_angle':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
#                                      '5':{'type':'catalog','name':'Catalogue_A','module':[1.75,1.75],'transverse_pressure_angle':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
#                                      '6':{'type':'catalog','name':'Catalogue_A','module':[2,2],'transverse_pressure_angle':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
#                                      '10':{'type':'optim','name':'Rack_A','module':[1,2],'transverse_pressure_angle':[18,22],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
#                                      '11':{'type':'optim','name':'Rack_B','module':[1,2],'transverse_pressure_angle':[18,22],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
#                                      }},
#          {'type':'rack_choice','data':{'0':['10'],'1':['11'],'2':['11']}}
#          )

#input_dict=(
#          {'type':'gear_speed','data':{'0':[7000,9000],'1':[16000,20000],'2':[11000,15000],'3':[17000,19000],'4':[15000,19000],'5':[15000,17000],'6':[4000,17000]}},
#          {'type':'gear_set','data':[(0,1),(1,2),(2,3),(3,4),(4,5),(5,6)]},
#          {'type':'center_distance','data':[[148.6,158.6],[130.1,140.1],[130.1,140.1],[130.1,140.1],[130.1,140.1],[130.1,140.1]]},
#          {'type':'Z','data':{'0':[30,100],'1':[30,100],'2':[30,100],'3':[30,100],'4':[30,100],'5':[30,100],'6':[30,100]}},
#          {'type':'transverse_pressure_angle','data':[[15,30],[15,30],[15,30],[15,30],[15,30],[15,30]]},
#          {'type':'helix_angle','data':{'0':[20,30]}},
#          {'type':'gear_width','data':{'0':[40,50],'1':[10,20]}},
#          {'type':'frequency','data':[[0,0]]},
#          {'type':'coefficient_profile_shift','data':{'0':[-1.2,1.2],'1':[-1.2,1.2]}},
#          {'type':'rack_list','data':{'0':{'type':'catalog','name':'Catalogue_A','module':[2.54,2.54],'transverse_pressure_angle_rack':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
#                                      }},
#          {'type':'rack_choice','data':{'0':['0'],'1':['0'],'2':['0'],'3':['0'],'4':['0'],'5':['0'],'6':['0']}}
#          )
erreur=0.02
list_cd=[[148.6,148.6+10],[130.1,130.1+10],[137,137+10],[135,135+10],[142,142+10],[145.4,145.4+10]]
list_gear_set=[(0,1),(1,2),(2,3),(3,4),(4,5),(5,6)]
#list_speed={'0':[7000,9000],'1':[16000,20000],'2':[11000,15000],'3':[17000,19000],'4':[15000,19000],'5':[15000,17000],'6':[4000,17000]}
list_speed={'0':[9000*(1-erreur),9000],'1':[20000*(1-erreur),20000],'2':[11000,11000*(1+erreur)],'3':[17000,17000*(1+erreur)],'4':[1000,30000],'5':[15000,15000*(1+erreur)],'6':[1000,30000]}
list_Z={}
module_safran=2.54

nb_gear=len(list_speed.keys())
nb_gset=len(list_cd)
for i,gs in enumerate(list_gear_set):
    cd_min=list_cd[i][0]
    cd_max=list_cd[i][1]
    demul_min=list_speed[str(gs[0])][0]/list_speed[str(gs[1])][1]
    demul_max=list_speed[str(gs[0])][1]/list_speed[str(gs[1])][0]
    DF1_max=2*cd_max/(1+demul_min)
    Z1_max=int(DF1_max/module_safran)+1
    DF2_max=2*cd_max*demul_max/(1+demul_max)
    Z2_max=int(DF2_max/module_safran)+1
    DF1_min=2*cd_min/(1+demul_max)
    Z1_min=int(DF1_min/module_safran)-1
    DF2_min=2*cd_min*demul_min/(1+demul_min)
    Z2_min=int(DF2_min/module_safran)-1
    if str(gs[0]) not in list_Z.keys():
        list_Z[str(gs[0])]=[Z1_min,Z1_max]
    else:
        list_Z[str(gs[0])]=[min(Z1_min,list_Z[str(gs[0])][0]),max(Z1_max,list_Z[str(gs[0])][1])]
    if str(gs[1]) not in list_Z.keys():
        list_Z[str(gs[1])]=[Z2_min,Z2_max]
    else:
        list_Z[str(gs[1])]=[min(Z2_min,list_Z[str(gs[1])][0]),max(Z2_max,list_Z[str(gs[1])][1])]


input_dict=(
#          {'type':'gear_speed','data':{'0':[9000*(1-erreur),9000],'1':[20000*(1-erreur),20000],'2':[11000,11000*(1+erreur)],'3':[17000,17000*(1+erreur)],'4':[0,30000],'5':[15000,15000*(1+erreur)],'6':[0,30000]}},
          {'type':'gear_speed','data':list_speed},
          {'type':'gear_set','data':list_gear_set},
          {'type':'center_distance','data':list_cd},
          {'type':'Z','data':list_Z},
          {'type':'transverse_pressure_angle','data':[[15,30]*nb_gset]},
          {'type':'helix_angle','data':{'0':[20,30]}},
          {'type':'gear_width','data':{'0':[40,50],'1':[10,20]}},
          {'type':'frequency','data':[[0,0]]},
          {'type':'coefficient_profile_shift','data':{'0':[-1.2,1.2],'1':[-1.2,1.2]}},
          {'type':'rack_list','data':{'0':{'type':'catalog','name':'Catalogue_A','module':[module_safran,module_safran],'transverse_pressure_angle_rack':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
                                      }},
          {'type':'rack_choice','data':{'0':['0'],'1':['0'],'2':['0'],'3':['0'],'4':['0'],'5':['0'],'6':['0']}}
          )

GA_wizard=gears.GearAssemblyOptimizerWizard(input_dict)
GA_wizard.SearchCenterLine()

compt=0
solutions=[]
sort=[]
for i,el in enumerate(GA_wizard.solutions_search):
    Z_data=GA_wizard.solutions_search[i]['Z_data'].copy()
    cd_data=GA_wizard.solutions_search[i]['center_distance'].copy()
    valid=True
    cd_input=[]
    fonctionnel=0
    for v in input_dict:
        if v['type']=='center_distance':
            for j,v2 in enumerate(v['data']):
#                print(cd_data,v2,cd_data[i][1]<v2[0],valid)
                if (cd_data[j][0]*0.99)<(v2[0]):
                    valid=False
                fonctionnel+=cd_data[j][0]-v2[0]
                cd_input.append([cd_data[j][0]*(0.99),cd_data[j][0]*1.1])
#            print(v['data'],cd_data,Z_data,valid)
    if valid:
        sort.append([fonctionnel,cd_input,Z_data])
        compt+=1
print('Nombre de solution avec entraxe supérieur au CDC:',compt)

sort_np=npy.array(sort)
for ind in npy.argsort(sort_np[:,0]):
    cd_input=sort[ind][1]
    Z_data=sort[ind][2]
    input_dict2=(
#              {'type':'gear_speed','data':{'0':[9000*(1-erreur),9000],'1':[20000*(1-erreur),20000],'2':[11000,11000*(1+erreur)],'3':[17000,17000*(1+erreur)],'4':[0,30000],'5':[15000,15000*(1+erreur)],'6':[0,30000]}},
              {'type':'gear_speed','data':list_speed},
              {'type':'gear_set','data':list_gear_set},
              {'type':'center_distance','data':cd_input},
              {'type':'Z','data':Z_data},
              {'type':'transverse_pressure_angle','data':[[15,30],[15,30],[15,30],[15,30],[15,30],[15,30]]},
              {'type':'helix_angle','data':{'0':[20,30]}},
              {'type':'gear_width','data':{'0':[40,50],'1':[10,20]}},
              {'type':'frequency','data':[[0,0]]},
              {'type':'coefficient_profile_shift','data':{'0':[-1.2,1.2],'1':[-1.2,1.2]}},
              {'type':'rack_list','data':{'0':{'type':'catalog','name':'Catalogue_A','module':[module_safran,module_safran],'transverse_pressure_angle_rack':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
                                          }},
              {'type':'rack_choice','data':{'0':['0'],'1':['0'],'2':['0'],'3':['0'],'4':['0'],'5':['0'],'6':['0']}}
              )
    GA_wizard2=gears.GearAssemblyOptimizerWizard(input_dict2)
    GA_wizard2.Optimize()
#        print(GA_wizard2.solutions[-1].center_distance)
    if len(GA_wizard2.solutions)>0:
        valid2=True
        for v in input_dict:
            if v['type']=='center_distance':
                for j,v2 in enumerate(v['data']):
                    if (GA_wizard2.solutions[-1].center_distance[j]*1e3)<(v2[0]):
                        valid2=False
        if valid2:
            solutions.append(GA_wizard2.solutions[-1])
            compt+=1
            print(solutions[-1].center_distance,Z_data)
            break

    

#print(GA_wizard.solutions[0].Dict())
#GA_wizard2.solutions[-1].SVGGearSet('name.txt',{'0':[0,0],'2':[0.2,0],'6':[0.5,0]})
#print(GA_wizard.solutions[0].linear_backlash[0]*GA_wizard.solutions[0].DF[0][1]/2)

#gear_set=[(0, 1)]
#nb_gear=npy.array(gear_set).max()+1
#gear_graph=nx.Graph()
#gear_graph.add_nodes_from(range(nb_gear))
#gear_graph.add_edges_from(gear_set)
#sol={'Z': {0: 10, 1: 11}, 
#     'gear_set': [(0, 1)], 
#     'gear_graph': gear_graph, 
#     'center_distance': [0.05], 
#     'transverse_pressure_angle': 15/180*npy.pi,
#     'coefficient_profile_shift': {'0': -0.1, '1': -0.01}}
#ga=gears.GearAssembly(**sol)
#ga.SVGGearSet('name.txt',{'0':[0,0],'1':[0.05,0]})
#print(ga.linear_backlash)

# ### Sorties Obj et CSV
# results=gears.GearAssemblyOptimizationResults(GA_wizard.solutions,input_dict)
# results.CSVExport('data.csv','w')
#
# ### Sorties graphiques SVG
# print('Nb solutions:',len(results.solutions))
# for i,ga in enumerate(results.solutions):
#     ga.Gear1.rack.SVGExport(5,'rack1-s{}.html'.format(str(i)))
#     ga.Gear2.rack.SVGExport(5,'rack2-s{}.html'.format(str(i)))
#     ga.Gear1.GearGenerationSVGExport('Cremaillere_Z1-s{}.html'.format(str(i)))
#     ga.Gear2.GearGenerationSVGExport('Cremaillere_Z2-s{}.html'.format(str(i)))
#     ga.MeshingSVGExport('Creation_Dent1-s{}.html'.format(str(i)),'Z1')
#     ga.MeshingSVGExport('Creation_Dent2-s{}.html'.format(str(i)),'Z2')
#     ga.SVGExport()
#     r=ga.FreeCADExport('GearAssembly_{}'.format(i),(0,0),(ga.center_distance,0),
#                      'python','/usr/lib/freecad/lib',['fcstd','stl'])
#     print(r)
#
#
# d=results.Dict()
# #print(d)
# import json
# json.dumps(d)
