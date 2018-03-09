import math
import mechanical_components.clutches as clutches
import numpy as np
import matplotlib.pyplot as plt
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 17:29:05 2018

@author: jezequel
"""

verin = clutches.HydraulicCylinder(0.030, 0.050, 0.100, 0.005, 0, 0, 0)
emb = clutches.Clutch(0.086, 0.125, 0.0006, 4, 0.062, 875, 1.5*(1/6)*10**-4, 0, verin)

regime = np.linspace(0, 2500*math.pi/30, 100)
test = emb.DragTorque(regime, 0)
#print(test)

"""
plt.plot(regime*30/math.pi, test)
plt.show()
"""

## Export Freecad vérin
primitives = []
primitives.extend(verin.chamber_volume)
primitives.extend(verin.piston_volume)
model_verin = vm.VolumeModel(primitives)
resp=model_verin.FreeCADExport('python','cylinder','/usr/lib/freecad/lib/',['stl','fcstd'])

# Export Freecad embrayage
primitives = []

primitives.extend(emb.plate_volume)
model_emb = vm.VolumeModel(primitives)
resp=model_emb.FreeCADExport('python','plate','/usr/lib/freecad/lib/',['stl','fcstd'])

