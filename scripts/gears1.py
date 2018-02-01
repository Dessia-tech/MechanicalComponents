import mechanical_components.gears as gears

input_dict=({'type':'Z1','nom':30},
          {'type':'Z2','nom':51},
          {'type':'center_distance','min':60,'max':100},
          {'type':'transverse_pressure_angle','min':15,'max':25},
          {'type':'helix_angle','min':20,'max':30},
          {'type':'gear_width','min':20,'max':40},
          {'type':'maximum_torque','nom':50},
          {'type':'coefficient_profile_shift1','min':-1.2,'max':1.2},
          {'type':'coefficient_profile_shift2','min':-1.2,'max':1.2},
#          {'type':'material1','nom':'acier_allie_cementation'},
#          {'type':'material2','nom':'acier_allie_cementation'},
#          {'type':'nb_cycle1','nom':10000000},
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
    ga.FreeCADExport('Assembly_{}'.format(i),(0,0),(ga.center_distance,0),'python','/usr/lib/freecad/lib',['fcstd','stl'])

d=results.Dict()
print(d)

