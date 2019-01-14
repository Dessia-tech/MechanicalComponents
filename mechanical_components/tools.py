#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
import numpy as npy

def StringifyDictKeys(d):
    if type(d) == list or type(d) == tuple:
        new_d = []
        for di in d:
            new_d.append(StringifyDictKeys(di))
        
    elif type(d) ==dict:
        new_d = {}
        for k,v in d.items():
            new_d[str(k)] = StringifyDictKeys(v)
    else:
        return d
    return new_d

def EnforceBuiltinTypes(variable):
    type_variable = type(variable)

    if type_variable == npy.int64:
        variable2 = int(variable)
    elif type_variable == npy.float64:
        variable2 = float(variable)
    elif type_variable == dict:
        variable2 = {}
        for k,v in variable.items():
            variable2[k] = EnforceBuiltinTypes(v)
    elif (type_variable == list) or (type_variable == tuple):
        variable2 = [EnforceBuiltinTypes(vi) for vi in variable]
    elif type_variable in [int, float]:
        variable2 = variable
    else:
        print(type_variable)
        raise NotImplementedError
            
    return variable2