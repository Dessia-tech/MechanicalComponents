import mechanical_components.gears_assembly as gears
import numpy as npy

##########################
# script pour AGB Safran
##########################

#definition data
erreur=0.02
list_cd=[[148.6,148.6+10],[130.1,130.1+10],[137,137+10],[135,135+10],[142,142+10],[145.4,145.4+10]]
list_gear_set=[(0,1),(1,2),(2,3),(3,4),(4,5),(5,6)]
#list_speed={'0':[7000,9000],'1':[16000,20000],'2':[11000,15000],'3':[17000,19000],'4':[15000,19000],'5':[15000,17000],'6':[4000,17000]}
list_speed={'0':[9000*(1-erreur),9000],'1':[20000*(1-erreur),20000],'2':[11000,11000*(1+erreur)],'3':[17000,17000*(1+erreur)],'4':[1000,30000],'5':[15000,15000*(1+erreur)],'6':[1000,30000]}
list_Z={}
module_safran=2.54

#nombre de dents adaptatif
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


#recherche exhaustive
list_rack={}
for k,v in list_speed.items():
    list_rack[k]=['0']
input_dict=(
#          {'type':'gear_speed','data':{'0':[9000*(1-erreur),9000],'1':[20000*(1-erreur),20000],'2':[11000,11000*(1+erreur)],'3':[17000,17000*(1+erreur)],'4':[0,30000],'5':[15000,15000*(1+erreur)],'6':[0,30000]}},
          {'type':'gear_speed','data':list_speed},
          {'type':'gear_set','data':list_gear_set},
          {'type':'center_distance','data':list_cd},
          {'type':'Z','data':list_Z},
          {'type':'transverse_pressure_angle','data':[[15,30]*nb_gset]},
          {'type':'helix_angle','data':{'0':[20,30]}},
          {'type':'gear_width','data':{'0':[40,50]}},
          {'type':'frequency','data':[[0,0]]},
          {'type':'coefficient_profile_shift','data':{'0':[-1.2,1.2]}},
          {'type':'rack_list','data':{'0':{'type':'catalog','name':'Catalogue_A','module':[module_safran,module_safran],'transverse_pressure_angle_rack':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
                                      }},
          {'type':'rack_choice','data':list_rack}
          )
GA_wizard=gears.GearAssemblyOptimizerWizard(input_dict)
GA_wizard.SearchCenterLine()

#recherche des solutions ayant des entraxes supérieures au CDC Safran
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

#optimisation de la denture optimale (ayant les entraxes les plus proches du CDC)
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
              {'type':'gear_width','data':{'0':[40,50]}},
              {'type':'frequency','data':[[0,0]]},
              {'type':'coefficient_profile_shift','data':{'0':[-1.2,1.2]}},
              {'type':'rack_list','data':{'0':{'type':'catalog','name':'Catalogue_A','module':[module_safran,module_safran],'transverse_pressure_angle_rack':[20,20],'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]},
                                          }},
              {'type':'rack_choice','data':list_rack}
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
            print('Liste des entraxes obtenues:',solutions[-1].center_distance)
            print('Liste des dents retenues:',Z_data)
            break
        else:
            print('Solution avec entraxe non conforme')

    
solutions[-1].SVGGearSet('name.txt',{'0':[0,0],'2':[0.2,0],'6':[0.5,0]})
