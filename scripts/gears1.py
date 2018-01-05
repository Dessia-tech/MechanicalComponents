import mechanical_components.gears as gears

input_dict=({'type':'ratio'},
          {'type':'Z1','min':30,'max':30},
          {'type':'Z2','min':45,'max':45},
          {'type':'center_distance','min':60,'max':70},
          {'type':'helix_angle','nom':0.3},
          {'type':'gear_width','nom':20})

F1=gears.GearAssemblyOptimizerWizard(input_dict)

### Sorties Obj et CSV
E1=gears.GearAssemblyOptimizationResults(F1.solutions,input_dict,'Famille_A')
E1.CSVExport('data.csv','Famille_A')

### Sorties graphiques SVG
print('Nb solutions:',len(F1.solutions))
for i,AG1 in enumerate(F1.solutions):
    AG1.SVGExport('Assembly_{}.html'.format(i),(0,0),(AG1.center_distance,0))
    AG1.Gear1.rack.SVGExport(5,'rack1-s{}.html'.format(str(i)))
    AG1.Gear2.rack.SVGExport(5,'rack2-s{}.html'.format(str(i)))
    AG1.Gear1.GearGenerationSVGExport('Cremaillere_Z1-s{}.html'.format(str(i)))
    AG1.Gear2.GearGenerationSVGExport('Cremaillere_Z2-s{}.html'.format(str(i)))
    AG1.MeshingSVGExport('Creation_Dent1-s{}.html'.format(str(i)),'Z1')
    AG1.MeshingSVGExport('Creation_Dent2-s{}.html'.format(str(i)),'Z2')
    
### Sorties data
results=[]
for GA in F1.solutions:
    results.append(GA.Dict())
