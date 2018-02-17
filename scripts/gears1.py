import mechanical_components.gears as gears

input_dict=(
#          {'type':'ratio','nom':0.3},
          {'type':'Z1','min':20,'max':20},
          {'type':'Z2','min':50,'max':50},
          {'type':'center_distance','min':30,'max':100},
          {'type':'transverse_pressure_angle','min':15,'max':30},
          {'type':'helix_angle','min':20,'max':30},
          {'type':'gear_width','min':20,'max':30},
          {'type':'maximum_torque','nom':50},
          {'type':'coefficient_profile_shift1','min':-1.2,'max':1.2},
          {'type':'coefficient_profile_shift2','min':-1.2,'max':1.2},
          {'type':'material1','nom':'hardened_alloy_steel'},
          {'type':'material2','nom':'hardened_alloy_steel'},
          {'type':'nb_cycle1','nom':10000000},
          {'type':'rim1','nom':'rim_gear'},
          {'type':'rim2','nom':'rim_gear'},
          {'type':'alpha_rim1','nom':1.1},
          {'type':'alpha_rim2','nom':1.1},
          {'type':'boring_diameter1','nom':10},
          {'type':'boring_diameter2','nom':30},
          {'type':'thickness_rim1','nom':6},
          {'type':'thickness_rim2','nom':6},
          )

#constraint_dict=({},)

GA_wizard=gears.GearAssemblyOptimizerWizard(input_dict)
if GA_wizard.ok:
    GA_wizard.Optimize()

### Sorties Obj et CSV
results=gears.GearAssemblyOptimizationResults(GA_wizard.solutions,input_dict)
results.CSVExport('data.csv','w')

### Sorties graphiques SVG
print('Nb solutions:',len(results.solutions))
for i,ga in enumerate(results.solutions):
    ga.SVGExport('Assembly_{}.html'.format(i),(0,0),(ga.center_distance,0))
    ga.Gear1.rack.SVGExport(5,'rack1-s{}.html'.format(str(i)))
    ga.Gear2.rack.SVGExport(5,'rack2-s{}.html'.format(str(i)))
    ga.Gear1.GearGenerationSVGExport('Cremaillere_Z1-s{}.html'.format(str(i)))
    ga.Gear2.GearGenerationSVGExport('Cremaillere_Z2-s{}.html'.format(str(i)))
    ga.MeshingSVGExport('Creation_Dent1-s{}.html'.format(str(i)),'Z1')
    ga.MeshingSVGExport('Creation_Dent2-s{}.html'.format(str(i)),'Z2')
    r=ga.FreeCADExport('GearAssembly_{}'.format(i),(0,0),(ga.center_distance,0),
                     'python','/usr/lib/freecad/lib',['fcstd','stl'])
    print(r)
    

d=results.Dict()
#print(d)
import json
json.dumps(d)
