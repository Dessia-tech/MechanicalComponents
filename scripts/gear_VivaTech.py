import mechanical_components.meshes as meshes
import numpy as npy

input_dict=(
          {'type':'ratio','nom':0.3,'erreur':0.04},
          {'type':'Z1','min':20,'max':40},
          {'type':'Z2','min':30,'max':70},
          {'type':'center_distance','min':60,'max':60},
          {'type':'transverse_pressure_angle','min':15,'max':30},
          {'type':'helix_angle','min':20,'max':20},
          {'type':'gear_width','min':20,'max':20},
          {'type':'maximum_torque','nom':1},
          {'type':'coefficient_profile_shift1','min':-1,'max':1},
          {'type':'coefficient_profile_shift2','min':-1,'max':1},
#          {'type':'material1','nom':'hardened_alloy_steel'},
#          {'type':'material2','nom':'hardened_alloy_steel'},
#          {'type':'nb_cycle1','nom':10000000},
          {'type':'rim1','nom':'shaft_gear'},
          {'type':'rim2','nom':'shaft_gear'},
#          {'type':'alpha_rim1','nom':1.1},
#          {'type':'alpha_rim2','nom':1.1},
#          {'type':'boring_diameter1','nom':10},
#          {'type':'boring_diameter2','nom':30},
#          {'type':'thickness_rim1','nom':6},
#          {'type':'thickness_rim2','nom':6},
          )

#constraint_dict=({},)

GA_wizard=meshes.MeshesAssemblyOptimizerWizard(input_dict)
if GA_wizard.ok:
    GA_wizard.Optimize()

### Sorties Obj et CSV
results=meshes.GearAssemblyOptimizationResults(GA_wizard.solutions,input_dict)
#results.CSVExport('data.csv','w')

### Sorties graphiques SVG
#print('Nb solutions:',len(results.solutions))
min_backlash=[]
for i,ga in enumerate(results.solutions):
    min_backlash.append([i,ga.linear_backlash])
sort=npy.argsort(min_backlash,0)
    
solu=results.solutions[sort[0][1]]
r=solu.FreeCADExport('GearAssembly',(0,0),(solu.center_distance,0),
                 'python','/usr/lib/freecad/lib',['fcstd','stl'])


d=results.Dict()
#print(d)
import json
json.dumps(d)
