#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 16:31:51 2020

@author: launay
"""
import mechanical_components.planetary_gears as pg

sun=pg.Planetary('sun',36,'Sun')
sun_2=pg.Planetary('sun_2',60,'Sun')
ring= pg.Planetary('ring',84,'Ring')
planet_carrier= pg.PlanetCarrier('planet_carrier')
planet_1=pg.Planet('planet_1','Simple',12)
planet_2=pg.Planet('planet_2','Double',12)
planet_3=pg.Planet('planet_3','Double',16)
connexion=[sun,planet_1,'GE'],[planet_1,planet_2,'GE'],[planet_2,ring,'GE'],[planet_2,planet_3,'D'],[planet_3,sun_2,'GI']

planetary_gears_1= pg.PlanetaryGear('pl_1', [sun,ring,sun_2], [planet_1,planet_2,planet_3], planet_carrier,connexion)



planetary_gears_1.plot_cinematic_graph()