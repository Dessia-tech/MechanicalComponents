import numpy as npy
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math
from scipy.linalg import norm,solve,LinAlgError
from scipy.optimize import *
from scipy.interpolate import splprep, splev
from sympy import *
import itertools

import mechanical_components.LibSvg as LibSvg
import mechanical_components.Gears as Gears

dico_ref=({'type':'ratio'},
          {'type':'Z1','min':30,'max':30},
          {'type':'Z2','min':45,'max':45},
          {'type':'center_distance','min':60,'max':70},
          {'type':'helix_angle','nom':0.3},
          {'type':'gear_width','nom':20})
F1=Gears.FrontGearAssembly(dico_ref)

### Sorties Obj et CSV
E1=Gears.ExportGearAssembly(F1.solutions,dico_ref,'Famille_A')
E1.CSVExport('data.csv','Famille_A')

### Sorties graphiques SVG
print('Nb solutions:',len(F1.solutions))
for i,AG1 in enumerate(F1.solutions):
    AG1.SVGExportGearAssembly('Assembly_{}.html'.format(i),(0,0),(AG1.center_distance,0))
    AG1.Gear1.rack.SVGExportRack(5,'rack1-s{}.html'.format(str(i)))
    AG1.Gear2.rack.SVGExportRack(5,'rack2-s{}.html'.format(str(i)))
    AG1.Gear1.SVGExportRackGenere('Cremaillere_Z1-s{}.html'.format(str(i)))
    AG1.Gear2.SVGExportRackGenere('Cremaillere_Z2-s{}.html'.format(str(i)))
    AG1.SVGExportDetailDent('Creation_Dent1-s{}.html'.format(str(i)),'Z1')
    AG1.SVGExportDetailDent('Creation_Dent2-s{}.html'.format(str(i)),'Z2')
    
### Sorties data
results=[]
for GA in F1.solutions:
    results.append(GA.Dict())
