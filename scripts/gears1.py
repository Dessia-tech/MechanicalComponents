import mechanical_components.gears as gears

input_dict=({'type':'ratio'},
          {'type':'Z1','min':30,'max':30},
          {'type':'Z2','min':52,'max':52},
          {'type':'center_distance','min':30,'max':70},
          {'type':'helix_angle','nom':0.3},
          {'type':'gear_width','nom':20})



GA_wizard=gears.GearAssemblyOptimizerWizard(input_dict)
if GA_wizard.error:
    GA_wizard.Optimize()


### Sorties Obj et CSV
results=gears.GearAssemblyOptimizationResults(GA_wizard.solutions,input_dict)
results.CSVExport('data.csv','w')



### Sorties graphiques SVG
print('Nb solutions:',len(results.solutions))
for i,ga in enumerate(results.solutions):
#    AG1.SVGExport('Assembly_{}.html'.format(i),(0,0),(AG1.center_distance,0))
#    AG1.Gear1.rack.SVGExport(5,'rack1-s{}.html'.format(str(i)))
#    AG1.Gear2.rack.SVGExport(5,'rack2-s{}.html'.format(str(i)))
#    AG1.Gear1.GearGenerationSVGExport('Cremaillere_Z1-s{}.html'.format(str(i)))
#    AG1.Gear2.GearGenerationSVGExport('Cremaillere_Z2-s{}.html'.format(str(i)))
#    AG1.MeshingSVGExport('Creation_Dent1-s{}.html'.format(str(i)),'Z1')
#    AG1.MeshingSVGExport('Creation_Dent2-s{}.html'.format(str(i)),'Z2')
    ga.FreeCADExport('Assembly_{}'.format(i),(0,0),(ga.center_distance,0))
    
print(results.Dict())

