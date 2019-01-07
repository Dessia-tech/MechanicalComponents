#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 10:42:56 2018

@author: Pierrem
"""

import numpy as npy
import math as mt
from scipy import interpolate
import os
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math
from scipy.linalg import norm
from scipy.optimize import minimize,fsolve
from scipy.interpolate import splprep, splev
#from sympy import *
import itertools
from jinja2 import Environment, PackageLoader, select_autoescape

import mechanical_components.LibSvgD3 as LibSvg
import mechanical_components.bearings as bearings

import persistent
import pandas
from pandas.plotting import scatter_matrix

with pkg_resources.resource_stream(pkg_resources.Requirement('mechanical_components'),'mechanical_components/catalogs/rules_rlts_SKF.csv') as rules_rlts_skf:
    rules_rlts_skf=pandas.read_csv(rules_rlts_skf)
    df_rules_dict=rules_rlts_skf.to_dict()
    
dico_rules={}
for ind in rules_rlts_skf.index:
    varx=rules_rlts_skf['x'][ind]
    vary=rules_rlts_skf['y'][ind]
    typ=rules_rlts_skf['type'][ind]
    a=rules_rlts_skf['a'][ind]
    b=rules_rlts_skf['b'][ind]
    if (varx,vary,typ) not in dico_rules.keys():
        dico_rules[(varx,vary,typ)]=[a,b]