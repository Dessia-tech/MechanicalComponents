import mechanical_components.gears_assembly as gears

input_dict=(
          {'type':'gear_speed','data':{'0':[1000,1200],'2':[4000,4500],'5':[3000,3200],'3':[2000,2200]}},
          {'type':'gear_set','data':((0,1),(1,2),(3,5),(3,4),(1,3))},
          {'type':'center_distance','data':([30,50],[30,50],[30,50],[30,50],[30,50])},
          {'type':'Z','data':{'0':[10,150],'1':[10,150]}},
          {'type':'transverse_pressure_angle','data':([15,30],[15,30],[15,30],[15,30],[15,30])},
          {'type':'helix_angle','data':{'0':[20,30]}},
          {'type':'gear_width','data':{'3':[40,50],'2':[10,20]}},
          {'type':'frequency','data':([1000,1200],[500,550],[3000,3500])},
          {'type':'coefficient_profile_shift','data':{'0':[-1,1],'1':[-1,1],'2':[-1,1]}},
          )

#constraint_dict=({},)

GA_wizard=gears.GearAssemblyOptimizerWizard(input_dict)
GA_wizard.Optimize()

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
