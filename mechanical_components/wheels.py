#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


"""

from dessia_common.core import DessiaObject

class Wheel(DessiaObject):
    def __init__(self, diameter, width, thickness):
        DessiaObject.__init__(self, diameter=diameter, width=width, thickness=thickness)
        

w = Wheel(0.2, 0.3, 0.002)