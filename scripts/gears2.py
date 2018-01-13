import mechanical_components.gears as gears

input_dict=({'type':'ratio','min':0.3,'max':0.4},
          {'type':'Z1','min':20,'max':60},
          {'type':'Z2','min':45,'max':70},
          {'type':'center_distance','min':50,'max':70},
          {'type':'helix_angle','nom':0.3},
          {'type':'gear_width','nom':20})


F1=gears.GearAssemblyOptimizerWizard(input_dict)
if F1.error:
    F1.Optimize()
    
if not F1.solutions==[]:
    E1=gears.GearAssemblyOptimizationResults(F1.solutions,input_dict)
    E1.CSVExport('data.csv','w')
    
input_dict=({'type':'ratio','min':0.4,'max':0.5},
          {'type':'Z1','min':20,'max':60},
          {'type':'Z2','min':45,'max':70},
          {'type':'center_distance','min':50,'max':70},
          {'type':'helix_angle','nom':0.3},
          {'type':'gear_width','nom':20})

F1=gears.GearAssemblyOptimizerWizard(input_dict)
if F1.error:
    F1.Optimize()
    
if not F1.solutions==[]:
    E1.Add(F1.solutions,input_dict)
    E1.CSVExport('data.csv','a')
    
