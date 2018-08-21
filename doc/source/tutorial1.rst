Define and Optimize a non ISO gear mesh
---------------------------------------

In this tutorial, we will define and optimize a non ISO gear mesh

.. image:: images/meshes1.png
  :height: 350px
  :alt: alternate text
  :align: center

The complete script can be found in scripts/meshes/meshes1.py

Python imports
^^^^^^^^^^^^^^

First, we import mechanical_components.optimization.meshes package and then one subpackage:
* numpy (http://www.numpy.org)

.. literalinclude:: ../../scripts/meshes/meshes1.py
  :lines: 1-2

In most scripts, the package is imported as meshes to make it shorter.

Input definition
^^^^^^^^^^^^^^^^

The minimum parameters to define one gear mesh are:
  * List define minimum and maximum center-distance
  * List of connected mesh
  * Dictionary of admissible speed

.. literalinclude:: ../../scripts/meshes/meshes1.py
  :lines: 4-17

MeshAssemblyOptimizer definition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../scripts/meshes/meshes1.py
  :lines: 19-28


Gear mesh optimization
^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../scripts/meshes/meshes1.py
  :lines: 30-32


Export CAD and SVG
^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../scripts/meshes/meshes1.py
  :lines: 33-35
