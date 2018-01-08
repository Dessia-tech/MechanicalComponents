import mechanical_components.gears as gears

input_dict=({'type':'ratio','nom':0.29},
          {'type':'Z1','min':20,'max':40},
          {'type':'Z2','min':45,'max':70},
          {'type':'center_distance','min':60,'max':70},
          {'type':'helix_angle','nom':0.3},
          {'type':'gear_width','nom':20})


F1=gears.GearAssemblyOptimizerWizard(input_dict)
if F1.error:
    F1.Optimize()
    
if not F1.solutions==[]:
    E1=gears.GearAssemblyOptimizationResults(F1.solutions,input_dict,'Famille_A')
    E1.CSVExport('data.csv','w','Famille_A')
    
input_dict=({'type':'ratio','nom':0.23},
          {'type':'Z1','min':20,'max':40},
          {'type':'Z2','min':45,'max':70},
          {'type':'center_distance','min':60,'max':70},
          {'type':'helix_angle','nom':0.3},
          {'type':'gear_width','nom':20})

F1=gears.GearAssemblyOptimizerWizard(input_dict)
if F1.error:
    F1.Optimize()
    
if not F1.solutions==[]:
    E1.Add(F1.solutions,input_dict,'Famille_A')
    E1.CSVExport('data.csv','a','Famille_A')
    
input_dict=({'type':'ratio','nom':0.31},
          {'type':'Z1','min':20,'max':40},
          {'type':'Z2','min':45,'max':70},
          {'type':'center_distance','min':60,'max':70},
          {'type':'helix_angle','nom':0.3},
          {'type':'gear_width','nom':20})

F1=gears.GearAssemblyOptimizerWizard(input_dict)
if F1.error:
    F1.Optimize()
    
if not F1.solutions==[]:
    E1.Add(F1.solutions,input_dict,'Famille_A')
    E1.CSVExport('data.csv','a','Famille_A')

input_dict=({'type':'ratio','nom':0.32},
          {'type':'Z1','min':20,'max':40},
          {'type':'Z2','min':45,'max':70},
          {'type':'center_distance','min':60,'max':70},
          {'type':'helix_angle','nom':0.3},
          {'type':'gear_width','nom':20})

F1=gears.GearAssemblyOptimizerWizard(input_dict)
if F1.error:
    F1.Optimize()
    
if not F1.solutions==[]:
    E1.Add(F1.solutions,input_dict,'Famille_A')
    E1.CSVExport('data.csv','a','Famille_A')
