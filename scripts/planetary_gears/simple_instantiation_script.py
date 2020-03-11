#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 12:43:48 2020

@author: launay
"""

import matplotlib.pyplot as plt
import mechanical_components.planetary_gears as pg
import networkx as nx

sun=pg.Planetary('sun',18,'sun')
ring= pg.Planetary('ring',25,'ring')
planet_carrier= pg.PlanetCarrier('planet_carrier')
planet_1=pg.Planet('planet_1','Simple',10)
planet_2=pg.Planet('planet_2','Double',10)
planet_3=pg.Planet('planet_3','Simple',5)

planetary_gears_1= pg.PlanetaryGears('pl_1', sun, ring, [planet_1], planet_carrier)
planetary_gears_1.plot()

# planetary_gears_2= pg.PlanetaryGears('pl_2', sun, ring, [planet_1,planet_2], planet_carrier)
# planetary_gears_3= pg.PlanetaryGears('pl_3', sun, ring, [planet_1,planet_2,planet_3], planet_carrier)
# print(planetary_gears_3.ratio())
# assembly_planetary_gear=pg.AssemblyPlanetaryGears('assembly_planetary_gear', [planetary_gears_1,planetary_gears_2,planetary_gears_3], [['pl_1_sun','pl_2_sun'],['pl_1_ring','pl_2_ring'],['pl_2_sun','pl_3_sun'],['pl_2_ring','pl_3_ring']])



# assembly_planetary_gear.plot()

# planetary_gear= nx.Graph()
# planetary_gear.add_nodes_from(['Pl',"R",'S','BT1','BT2','BT3','E1','E2','Pv1','Pv2','Pv3','Pv4','PC'])
# planetary_gear.add_edges_from([('BT1','Pv1',{'color':'blue'}),('Pv1','S'),('S','E1'),('E1','Pl'),('Pl','E2'),('E2','R'),('R','Pv2'),('Pv2','BT2'),('Pl','Pv3'),('Pv3','PC'),('PC','Pv4'),('Pv4','BT3')])
# # nx.draw_networkx(planetary_gear,nlist=[['Pl','Pv3','PC','Pv4','BT3'],['BT1','Pv1','S','E1','Pl','E2','R','Pv2','BT2']], with_labels=True)
# planetary_gear.add_nodes_from(['Pl','Pv'+str(4) ],)

# planetary_gear_assembly= nx.union(planetary_gear,planetary_gear,rename=('1','2'))
# planetary_gear_assembly.add_nodes_from(['EC1','EC2'])
# planetary_gear_assembly.add_edges_from([('1PC','EC1'),('EC1','2R'),('1S','EC2'),('EC2','2S')])
# nx.draw_kamada_kawai(planetary_gear_assembly, with_labels=True)
# planetary_gear_assembly2= nx.union(planetary_gear,planetary_gear_assembly,rename=('3',''))
# planetary_gear_assembly2.add_nodes_from(['EC3','EC4'])
# planetary_gear_assembly2.add_edges_from([('1PC','EC3'),('EC3','3R'),('1S','EC4'),('EC4','3S')])
# # nx.draw_kamada_kawai(planetary_gear_assembly2, with_labels=True)

