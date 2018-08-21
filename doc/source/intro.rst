First steps
===========

DessIA MechanicalComponents is a software for designing mechanical components
for:

 * Meshes
 * Bearings
 * Springs
 * Clutches

It comes with an open source-based version augmented by a connection to the
DessIA software platform (https://software.dessia.tech) by API

Installation
------------

The open source package is nothing more than a python3 package to install.


Project structure
-----------------

The python install script (setup.py) at the root level installs files from the
powertransmission folder

The scripts folder gives some examples of how using this software.

To sum up::

  powertransmission/ # Package source
  doc/ # Documentation
  |-- source/ # source files of doc
  |-- build/ # build files
  scripts/ # example scripts
  |-- aeronautics/ # Aeronautics scripts
  |-- automotive/ # Automotive scripts
  setup.py # Setup file
