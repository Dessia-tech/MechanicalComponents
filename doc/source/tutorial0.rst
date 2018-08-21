Define and Optimize a simple gear mesh
--------------------------------------

In this tutorial, we will define and optimize a gear mesh

.. image:: images/meshes0.png
   :height: 350px
   :alt: alternate text
   :align: center

The complete script can be found in scripts/meshes/meshes0.py

Python imports
^^^^^^^^^^^^^^

First, we import mechanical_components.optimization.meshes package and then one subpackage:
 * numpy (http://www.numpy.org)

.. literalinclude:: ../../scripts/meshes/meshes0.py
   :lines: 1-2

In most scripts, the package is imported as meshes to make it shorter.

Input definition
^^^^^^^^^^^^^^^^

The minimum parameters to define one gear mesh are:
 * List define minimum and maximum center-distance
 * List of connected mesh
 * Dictionary of admissible speed

.. literalinclude:: ../../scripts/meshes/meshes0.py
   :lines: 4-8

MeshAssemblyOptimizer definition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../scripts/meshes/meshes0.py
   :lines: 10-12

.. seealso::

  .. autoclass:: mechanical_components.optimization.meshes.MeshAssemblyOptimizer

Gear mesh optimization
^^^^^^^^^^^^^^^^^^^^^^

Automatic gear mesh optimize
****************************

.. literalinclude:: ../../scripts/meshes/meshes0.py
   :lines: 14-17

Sequential gear mesh optimize
*****************************

.. literalinclude:: ../../scripts/meshes/meshes0.py
   :lines: 19-20

Export CAD and SVG
^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../scripts/meshes/meshes0.py
   :lines: 22-26
