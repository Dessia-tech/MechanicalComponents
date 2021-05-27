#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 16:31:51 2020

@author: launay
"""
import mechanical_components.planetary_gears as pg
from dessia_api_client import Client
import plot_data
sun = pg.Planetary(36, 'Sun', 'sun')
sun_2 = pg.Planetary(60, 'Sun', 'sun_2')
ring = pg.Planetary(84, 'Ring', 'ring')
planet_carrier = pg.PlanetCarrier('planet_carrier')
planet_1 = pg.Planet( 12, 'planet_1')
planet_2 = pg.Planet( 12, 'planet_2')
planet_3 = pg.Planet( 16, 'planet_3')
connection = [pg.Connection([sun, planet_1], 'GE'), 
              pg.Connection([planet_1, planet_2], 'GE'), 
              pg.Connection([planet_2, ring], 'GE'),
              pg.Connection([planet_2, planet_3], 'D'), 
              pg.Connection([planet_3, sun_2], 'GI')]

planetary_gears_1 = pg.PlanetaryGear([sun, ring, sun_2], [planet_1, planet_2, planet_3], \
                                      planet_carrier, connection, 'pl_1')


plot_data.plot_canvas(planetary_gears_1.plot_kinematic_graph()[0])