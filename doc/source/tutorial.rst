Tutorial
========


Define and Optimize a gear mesh
-------------------------------

In this tutorial, we will define and optimize a gear mesh

.. image:: images/meshes1.png
   :height: 350px
   :alt: alternate text
   :align: center

The complete script can be found in scripts/meshes/meshes1.py

Python imports
^^^^^^^^^^^^^^

First, we import mechanical_components.optimization.meshes package and then one subpackage:
 * numpy (http://www.numpy.org)

.. code:: python

  import mechanical_components.optimization.meshes as meshes_opt
  import numpy as npy

In most scripts, the package is imported as meshes to make it shorter.

Input definition
^^^^^^^^^^^^^^^^

Mandatory parameters
********************

The minimum parameters to define one gear mesh are:
 * List define minimum and maximum center-distance
 * List of connected mesh
 * Dictionary of admissible speed

.. code:: python

  list_cd=[[0.117,0.16]]
  list_gear_set=[(5,1)]
  list_speed={5:[1000*npy.pi/30,1500*npy.pi/30],1:[4100*npy.pi/30,
                 4300*npy.pi/30]}

Optional parameters
*******************

.. code:: python

  list_rack={0:{'name':'Catalogue_A','module':[0.5*1e-3,2.6*1e-3],
                'transverse_pressure_angle_rack':[20/180*npy.pi,20/180*npy.pi],
                'coeff_gear_addendum':[1,1],'coeff_gear_dedendum':[1.25,1.25],
                'coeff_root_radius':[0.38,0.38],'coeff_circular_tooth_thickness':[0.5,0.5]}}
  list_rack_choice={5:[0],1:[0]}
  list_helix_angle={5:[0,0]}
  list_material={5:meshes_opt.hardened_alloy_steel}
  list_torque={1:186,5:'output'}
  list_cycle={1:1e12}

MeshAssemblyOptimizer definition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

  GA=meshes_opt.MeshAssemblyOptimizer(Z={},
                                 connections=list_gear_set,
                                 gear_speed=list_speed,
                                 center_distance=list_cd,
                                 rack_list=list_rack,
                                 rack_choice=list_rack_choice,
                                 helix_angle=list_helix_angle,
                                 material=list_material,
                                 torque=list_torque,
                                 cycle=list_cycle)

.. seealso::

  .. autoclass:: mechanical_components.optimization.meshes.MeshAssemblyOptimizer


Gear mesh optimization
^^^^^^^^^^^^^^^^^^^^^^

Automatic gear mesh optimize
****************************

.. code:: python

  GA.SearchOptimumCD(nb_sol=1,verbose=True)
  print('Nombre de solutions convergés:',len(GA.solutions))
  solution=GA.solutions[-1]
  solution.SVGExport('name.txt',{5:[0,0]})
  solution.FreeCADExport('Gears1')

Sequential gear mesh optimize
*****************************

  .. literalinclude:: ../../scripts/meshes/meshes1.py
